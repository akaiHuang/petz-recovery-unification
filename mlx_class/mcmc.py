"""
mcmc.py — MCMC parameter fitting interface and batch parameter evaluation.

KILLER FEATURE: Evaluate C_l for multiple parameter sets simultaneously on GPU.
Even with a simple Python loop, mlx_class avoids the massive per-run setup overhead
of CLASS. Future versions will use MLX vmap for true GPU parallelism.

Components:
    1. ParameterizedSolver — compute C_l for arbitrary LCDM parameters
    2. BatchEvaluator      — evaluate C_l for multiple parameter sets
    3. MetropolisMCMC      — Metropolis-Hastings sampler
    4. MockDataGenerator   — generate mock Planck-like data for testing
    5. MCMCDiagnostics     — chain analysis and plotting

Usage:
    python -m mlx_class.mcmc [--n_steps 5000] [--mock] [--plot]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import argparse
import time
import numpy as np
import mlx.core as mx

# Default Planck 2018 fiducial parameters
FIDUCIAL_PARAMS = {
    'h':         0.6736,
    'omega_b':   0.02237,
    'omega_cdm': 0.12,
    'A_s':       2.1e-9,
    'n_s':       0.9649,
    'tau_reio':  0.0544,
}

# Flat priors (physical bounds)
DEFAULT_PRIORS = {
    'h':         (0.5, 0.9),
    'omega_b':   (0.01, 0.04),
    'omega_cdm': (0.05, 0.25),
    'A_s':       (1.0e-9, 5.0e-9),
    'n_s':       (0.8, 1.1),
    'tau_reio':  (0.01, 0.15),
}

# Default step sizes (roughly 1% of parameter range, tuned for acceptance ~25%)
DEFAULT_STEP_SIZES = {
    'h':         0.005,
    'omega_b':   0.0003,
    'omega_cdm': 0.002,
    'A_s':       0.05e-9,
    'n_s':       0.005,
    'tau_reio':  0.005,
}


# ============================================================================
# Component 1: Parameterized solver — single parameter set -> C_l
# ============================================================================

class ParameterizedSolver:
    """
    Compute C_l TT for an arbitrary set of LCDM cosmological parameters.

    This wraps the full mlx_class pipeline (Background -> IMEX -> C_l)
    with custom parameters, overriding the module-level globals.

    Parameters
    ----------
    ell_values : array
        Multipole values at which to compute C_l.
    N_k : int
        Number of k-modes for the perturbation solver.
    mode : str
        'fast' (N_k=200, no ISW) or 'standard' (N_k=500, with ISW).
    """

    def __init__(self, ell_values=None, N_k=200, mode='fast'):
        if ell_values is None:
            # Sparse ell grid for speed — sufficient for MCMC chi2
            ell_values = np.unique(np.concatenate([
                np.arange(2, 30, 2),
                np.arange(30, 100, 5),
                np.arange(100, 500, 10),
                np.arange(500, 1500, 20),
                np.arange(1500, 2501, 30),
            ])).astype(int)
        self.ell_values = np.asarray(ell_values, dtype=int)
        self.N_k = N_k
        self.mode = mode

    def compute_cl(self, params):
        """
        Compute D_l^TT for a single parameter set.

        Parameters
        ----------
        params : dict
            Must contain: h, omega_b, omega_cdm, A_s, n_s.
            Optional: tau_reio (used for reionization damping).

        Returns
        -------
        ell_out : array (N_ell,)
        Dl : array (N_ell,) in muK^2
        """
        # Unpack parameters
        h_val = params['h']
        omega_b = params['omega_b']
        omega_cdm = params['omega_cdm']
        A_s_val = params['A_s']
        n_s_val = params['n_s']
        tau_reio = params.get('tau_reio', 0.0544)

        # Derived quantities
        H0_km_s_Mpc = h_val * 100.0
        c_km_s = 299792.458
        H0_Mpc = H0_km_s_Mpc / c_km_s

        Omega_b_val = omega_b / h_val**2
        Omega_c_val = omega_cdm / h_val**2
        Omega_r_val = 2.469e-5 * (1 + 0.2271 * 3.046) / h_val**2
        Omega_m_val = Omega_b_val + Omega_c_val
        Omega_L_val = 1.0 - Omega_m_val - Omega_r_val

        # Build a custom Background with these parameters
        bg = self._build_background(
            h_val, Omega_b_val, Omega_c_val, Omega_r_val,
            Omega_m_val, Omega_L_val, H0_Mpc, omega_b
        )

        # k-grid
        k_arr = np.geomspace(5e-5, 0.35, self.N_k).astype(np.float32)

        # Solve perturbations (IMEX)
        from .perturbations_implicit import ImplicitBoltzmannSolver
        solver = ImplicitBoltzmannSolver(bg, k_arr)
        result = solver.solve()
        Theta_0, Phi, v_b = result.source_at_recombination()

        source_SW = Theta_0 + Phi

        # Upsample source for C_l integration
        from .main import upsample_source
        k_fine, source_fine = upsample_source(
            k_arr, source_SW, bg.D_A, ell_max=int(self.ell_values[-1])
        )

        # C_l integration with custom A_s, n_s
        Cl, Dl = self._compute_cl_custom(
            source_fine, k_fine, self.ell_values, bg.D_A,
            A_s_val, n_s_val
        )

        # Reionization damping: exp(-2 * tau_reio) for l > ~10
        if tau_reio > 0:
            reio_damp = np.exp(-2.0 * tau_reio)
            Dl = Dl * reio_damp

        # Early ISW correction (optional, only in standard mode)
        if self.mode == 'standard':
            from .main import early_isw_template
            T_CMB = 2.7255
            Dl_ISW = early_isw_template(self.ell_values, Dl, bg)
            Dl = Dl + Dl_ISW

        return self.ell_values, Dl

    def _build_background(self, h_val, Omega_b, Omega_c, Omega_r,
                          Omega_m, Omega_L, H0_Mpc, omega_b_phys):
        """
        Build a Background object with custom cosmological parameters.

        We temporarily override the module-level globals in background.py,
        create and solve the Background, then restore the originals.
        This is safe because mlx_class is single-threaded.
        """
        import mlx_class.background as bg_mod

        # Save originals
        saved = {}
        attrs_to_override = [
            'h', 'H0_km_s_Mpc', 'H0_Mpc', 'H0_SI',
            'omega_b', 'omega_c', 'Omega_b', 'Omega_c', 'Omega_r',
            'Omega_m', 'Omega_L', 'a_eq', 'k_eq', 'A_s', 'n_s',
        ]
        for attr in attrs_to_override:
            saved[attr] = getattr(bg_mod, attr)

        try:
            # Override
            bg_mod.h = h_val
            bg_mod.H0_km_s_Mpc = h_val * 100.0
            c_km_s = 299792.458
            bg_mod.H0_Mpc = bg_mod.H0_km_s_Mpc / c_km_s
            Mpc_SI = 3.0856775814913673e22
            bg_mod.H0_SI = bg_mod.H0_km_s_Mpc * 1e3 / Mpc_SI
            bg_mod.omega_b = omega_b_phys
            bg_mod.omega_c = Omega_c * h_val**2
            bg_mod.Omega_b = Omega_b
            bg_mod.Omega_c = Omega_c
            bg_mod.Omega_r = Omega_r
            bg_mod.Omega_m = Omega_m
            bg_mod.Omega_L = Omega_L
            bg_mod.a_eq = Omega_r / Omega_m
            bg_mod.k_eq = np.sqrt(2 * Omega_m * bg_mod.H0_Mpc**2 / bg_mod.a_eq)
            bg_mod.k_D = 0.15 * (omega_b_phys / 0.022)**0.25

            from .background import Background
            bg = Background(khronon=False)
            bg.solve()
            return bg

        finally:
            # Restore originals
            for attr, val in saved.items():
                setattr(bg_mod, attr, val)

    def _compute_cl_custom(self, source_SW, k_arr, ell_values, D_A, A_s_val, n_s_val):
        """
        C_l integration with custom primordial spectrum parameters.
        Mirrors spectra.compute_cl but with A_s, n_s as arguments.
        """
        from scipy.special import spherical_jn

        T_CMB = 2.7255
        k_pivot = 0.05

        N_ell = len(ell_values)
        N_k = len(k_arr)
        x_arr = k_arr * D_A

        # Primordial power spectrum with custom A_s, n_s
        P_R = A_s_val * (k_arr / k_pivot)**(n_s_val - 1.0)

        # Bessel table
        jl = np.zeros((N_ell, N_k), dtype=np.float32)
        for il, ell in enumerate(ell_values):
            jl[il] = spherical_jn(int(ell), x_arr)

        # Transfer: Delta_l(k) = SW(k) * j_l(k * D_A)
        Delta_l = source_SW[None, :] * jl

        # C_l = 4*pi * int d(ln k) P_R |Delta_l|^2
        integrand = (P_R[None, :] * Delta_l**2).astype(np.float32)
        lnk = np.log(k_arr)
        dlnk = np.diff(lnk).astype(np.float32)

        integrand_gpu = mx.array(integrand)
        dlnk_gpu = mx.array(dlnk)
        mid = 0.5 * (integrand_gpu[:, :-1] + integrand_gpu[:, 1:])
        Cl_gpu = 4.0 * np.pi * mx.sum(mid * dlnk_gpu[None, :], axis=1)
        mx.eval(Cl_gpu)

        Cl = np.maximum(np.array(Cl_gpu), 0.0)
        ell_f = ell_values.astype(float)
        Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (T_CMB * 1e6)**2

        return Cl, Dl


# ============================================================================
# Component 1b: Sync gauge parameterized solver
# ============================================================================

class SyncParameterizedSolver:
    """
    Compute C_l TT using the synchronous gauge solver backend.

    This wraps solver_sync.run_sync_solver() which uses scipy's Radau
    solver for each k-mode sequentially. More accurate than the IMEX solver
    (full photon+neutrino hierarchy, proper LOS integration with ISW),
    but much slower (~420s for 180 k-modes vs ~2s for IMEX).

    Parameters
    ----------
    ell_values : array or None
        Multipole values at which to compute C_l. If None, uses the
        sync solver's default ell grid.
    N_k : int
        Number of k-modes (default: 100 for speed; production: 300).
    method : str
        ODE method for scipy solve_ivp (default: 'Radau').
    rtol, atol : float
        ODE tolerances.

    Notes
    -----
    - The sync gauge solver currently uses hardcoded cosmological parameters
      from background.py. Custom parameter support requires passing parameters
      through to the Background constructor (TODO).
    - For MCMC, this is too slow (~420s per evaluation with 180 k-modes).
      GPU batched version or multiprocessing needed for actual MCMC runs.
    - For single comparisons or validation, this provides the most accurate
      C_l from mlx_class.
    """

    def __init__(self, ell_values=None, N_k=100, method='Radau',
                 rtol=1e-6, atol=1e-9):
        self.ell_values = ell_values  # None = use sync solver default
        self.N_k = N_k
        self.method = method
        self.rtol = rtol
        self.atol = atol

    def compute_cl(self, params=None):
        """
        Compute D_l^TT using the sync gauge solver.

        Parameters
        ----------
        params : dict or None
            Cosmological parameters. Currently ignored — the sync gauge
            solver uses hardcoded params from background.py.
            TODO: pass parameters through to Background constructor.

        Returns
        -------
        ell_out : array (N_ell,)
        Dl : array (N_ell,) in muK^2
        """
        from .solver_sync import run_sync_solver

        # TODO: When parameter passthrough is implemented, override
        # background.py globals here (similar to ParameterizedSolver).
        if params is not None:
            import warnings
            warnings.warn(
                "SyncParameterizedSolver currently ignores custom parameters. "
                "Using hardcoded background.py values. "
                "Parameter passthrough is a TODO.",
                UserWarning, stacklevel=2
            )

        result = run_sync_solver(
            N_k=self.N_k,
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            verbose=False,
        )

        ell_out = result['ell']
        Dl = result['Dl']

        # If custom ell_values requested, interpolate
        if self.ell_values is not None:
            from scipy.interpolate import interp1d
            f_interp = interp1d(ell_out, Dl, kind='cubic',
                                fill_value='extrapolate')
            Dl = f_interp(self.ell_values)
            Dl = np.maximum(Dl, 0.0)
            ell_out = self.ell_values

        return ell_out, Dl


# ============================================================================
# Component 2: Batch parameter evaluation
# ============================================================================

class BatchEvaluator:
    """
    Evaluate C_l for multiple parameter sets.

    The key GPU advantage: even with a Python loop, mlx_class avoids the
    heavy per-run setup overhead of CLASS (~2s overhead vs ~10s for CLASS).

    Future: MLX vmap over parameter sets for true GPU parallelism.

    Parameters
    ----------
    ell_values : array or None
        Multipole values for output. None = use default sparse grid.
    N_k : int
        Number of k-modes per evaluation.
    mode : str
        'fast' or 'standard'.
    solver_type : str
        'production' (default, IMEX conformal Newtonian gauge) or
        'sync' (synchronous gauge, more accurate but ~200x slower).
    verbose : bool
        Print progress for each evaluation.
    """

    def __init__(self, ell_values=None, N_k=200, mode='fast',
                 solver_type='production', verbose=False):
        if solver_type == 'sync':
            self.solver = SyncParameterizedSolver(
                ell_values=ell_values, N_k=N_k)
            # SyncParameterizedSolver may return its own ell grid
            if ell_values is not None:
                self.ell_values = np.asarray(ell_values, dtype=int)
            else:
                # Run once to get the ell grid (expensive but needed for API)
                # Use a minimal N_k just to get the ell layout
                from .solver_sync import run_sync_solver
                _tmp = run_sync_solver(N_k=5, verbose=False)
                self.ell_values = _tmp['ell']
        else:
            self.solver = ParameterizedSolver(
                ell_values=ell_values, N_k=N_k, mode=mode)
            self.ell_values = self.solver.ell_values
        self.solver_type = solver_type
        self.verbose = verbose
        self._cache = {}
        self._cache_hits = 0
        self._cache_misses = 0

    def _param_key(self, params):
        """Create a hashable cache key from a parameter dict."""
        return tuple(sorted(
            (k, round(v, 12)) for k, v in params.items()
        ))

    def evaluate_single(self, params):
        """
        Evaluate C_l for a single parameter set.

        Parameters
        ----------
        params : dict
            Keys: h, omega_b, omega_cdm, A_s, n_s, tau_reio.

        Returns
        -------
        Dl : array (N_ell,) — D_l in muK^2
        """
        key = self._param_key(params)
        if key in self._cache:
            self._cache_hits += 1
            return self._cache[key]

        self._cache_misses += 1

        # Suppress solver print output for batch mode
        import io
        import contextlib
        if not self.verbose:
            with contextlib.redirect_stdout(io.StringIO()):
                _, Dl = self.solver.compute_cl(params)
        else:
            _, Dl = self.solver.compute_cl(params)

        self._cache[key] = Dl
        return Dl

    def evaluate(self, params_list):
        """
        Evaluate C_l for multiple parameter sets.

        Parameters
        ----------
        params_list : list of dicts
            Each dict has: h, omega_b, omega_cdm, A_s, n_s, tau_reio.

        Returns
        -------
        results : list of arrays
            Each element is D_l (N_ell,) in muK^2.
        """
        n_total = len(params_list)
        results = []
        t0 = time.time()

        for i, params in enumerate(params_list):
            if self.verbose:
                print(f"\n--- Evaluating parameter set {i+1}/{n_total} ---")
                for k, v in params.items():
                    print(f"  {k}: {v}")

            Dl = self.evaluate_single(params)
            results.append(Dl)

            if self.verbose and (i + 1) % 10 == 0:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                eta = (n_total - i - 1) / rate
                print(f"[BatchEvaluator] {i+1}/{n_total} done, "
                      f"{rate:.1f} evals/s, ETA: {eta:.0f}s")

        elapsed = time.time() - t0
        if n_total > 0:
            print(f"[BatchEvaluator] {n_total} evaluations in {elapsed:.1f}s "
                  f"({elapsed/n_total:.2f}s per eval, "
                  f"cache hits: {self._cache_hits})")

        return results

    def clear_cache(self):
        """Clear the parameter cache."""
        self._cache.clear()
        self._cache_hits = 0
        self._cache_misses = 0


# ============================================================================
# Component 3: MCMC sampler (Metropolis-Hastings)
# ============================================================================

class MetropolisMCMC:
    """
    Metropolis-Hastings MCMC sampler for CMB parameter estimation.

    Parameters
    ----------
    data_Dl : array (N_ell,)
        Observed D_l in muK^2.
    errors_Dl : array (N_ell,)
        1-sigma error bars on D_l.
    ell_values : array (N_ell,)
        Multipole values corresponding to data_Dl.
    evaluator : BatchEvaluator
        The C_l evaluator.
    priors : dict or None
        Flat prior bounds: {param_name: (lo, hi)}.
    """

    def __init__(self, data_Dl, errors_Dl, ell_values, evaluator, priors=None):
        self.data_Dl = np.asarray(data_Dl, dtype=np.float64)
        self.errors_Dl = np.asarray(errors_Dl, dtype=np.float64)
        self.ell_values = np.asarray(ell_values)
        self.evaluator = evaluator
        self.priors = priors if priors is not None else DEFAULT_PRIORS

        # Match data ells to evaluator ells
        self._ell_mask = np.isin(evaluator.ell_values, self.ell_values)
        self._data_mask = np.isin(self.ell_values, evaluator.ell_values)

        n_matched = np.sum(self._data_mask)
        n_data = len(self.ell_values)
        print(f"[MCMC] Data: {n_data} ells, matched to evaluator: {n_matched}")

        # Precompute inverse variance on matched ells
        matched_errors = self.errors_Dl[self._data_mask]
        self._inv_var = 1.0 / (matched_errors**2)
        self._matched_data = self.data_Dl[self._data_mask]
        self._n_dof = n_matched

    def chi2(self, params):
        """
        Compute chi^2 = sum((model - data)^2 / error^2).

        Parameters
        ----------
        params : dict
            Cosmological parameters.

        Returns
        -------
        chi2_val : float
        """
        # Check priors
        for name, (lo, hi) in self.priors.items():
            if name in params:
                if params[name] < lo or params[name] > hi:
                    return np.inf

        try:
            model_Dl = self.evaluator.evaluate_single(params)
            model_matched = model_Dl[self._ell_mask]
        except Exception as e:
            print(f"[MCMC] Evaluation failed: {e}")
            return np.inf

        residual = model_matched - self._matched_data
        return float(np.sum(residual**2 * self._inv_var))

    def log_likelihood(self, params):
        """Log-likelihood = -0.5 * chi^2."""
        chi2_val = self.chi2(params)
        if np.isinf(chi2_val):
            return -np.inf
        return -0.5 * chi2_val

    def run(self, initial_params, n_steps=1000, step_sizes=None,
            burnin_fraction=0.3, thin=1, seed=42, verbose=True):
        """
        Run Metropolis-Hastings MCMC chain.

        Parameters
        ----------
        initial_params : dict
            Starting point.
        n_steps : int
            Total number of MCMC steps (including burn-in).
        step_sizes : dict or None
            Gaussian proposal widths for each parameter.
        burnin_fraction : float
            Fraction of chain to discard as burn-in.
        thin : int
            Thinning factor (keep every thin-th sample).
        seed : int
            Random seed.
        verbose : bool
            Print progress every 10% of steps.

        Returns
        -------
        result : MCMCResult
        """
        rng = np.random.RandomState(seed)
        param_names = sorted(initial_params.keys())

        if step_sizes is None:
            step_sizes = {k: DEFAULT_STEP_SIZES.get(k, abs(v) * 0.01)
                          for k, v in initial_params.items()}

        # Initialize chain storage
        chain = np.zeros((n_steps, len(param_names)))
        chi2_chain = np.zeros(n_steps)
        accepted = np.zeros(n_steps, dtype=bool)

        # Current state
        current_params = dict(initial_params)
        current_chi2 = self.chi2(current_params)
        current_logL = -0.5 * current_chi2

        chain[0] = [current_params[k] for k in param_names]
        chi2_chain[0] = current_chi2
        accepted[0] = True

        n_accept = 0
        t0 = time.time()
        report_interval = max(1, n_steps // 10)

        if verbose:
            print(f"\n{'='*60}")
            print(f"MCMC: Metropolis-Hastings, {n_steps} steps")
            print(f"Parameters: {param_names}")
            print(f"Initial chi^2 = {current_chi2:.1f} "
                  f"(chi^2/dof = {current_chi2/self._n_dof:.2f})")
            print(f"{'='*60}")

        for step in range(1, n_steps):
            # Propose new parameters (Gaussian random walk)
            proposed = {}
            for name in param_names:
                proposed[name] = current_params[name] + \
                    rng.normal(0, step_sizes[name])

            # Evaluate
            proposed_chi2 = self.chi2(proposed)
            proposed_logL = -0.5 * proposed_chi2

            # Metropolis acceptance
            log_alpha = proposed_logL - current_logL
            if log_alpha > 0 or rng.uniform() < np.exp(log_alpha):
                current_params = proposed
                current_chi2 = proposed_chi2
                current_logL = proposed_logL
                accepted[step] = True
                n_accept += 1

            chain[step] = [current_params[k] for k in param_names]
            chi2_chain[step] = current_chi2

            if verbose and (step + 1) % report_interval == 0:
                elapsed = time.time() - t0
                rate = (step + 1) / elapsed
                accept_rate = n_accept / (step + 1)
                print(f"  Step {step+1}/{n_steps}: "
                      f"chi^2={current_chi2:.1f}, "
                      f"accept={accept_rate:.1%}, "
                      f"{rate:.1f} steps/s")

        elapsed = time.time() - t0
        accept_rate = n_accept / n_steps

        if verbose:
            print(f"\n{'='*60}")
            print(f"MCMC complete: {elapsed:.1f}s "
                  f"({n_steps/elapsed:.1f} steps/s)")
            print(f"Acceptance rate: {accept_rate:.1%}")
            print(f"Cache: {self.evaluator._cache_hits} hits, "
                  f"{self.evaluator._cache_misses} misses")
            print(f"{'='*60}")

        # Build result
        n_burnin = int(n_steps * burnin_fraction)
        return MCMCResult(
            chain=chain,
            chi2_chain=chi2_chain,
            accepted=accepted,
            param_names=param_names,
            n_burnin=n_burnin,
            thin=thin,
            accept_rate=accept_rate,
            elapsed=elapsed,
            true_params=initial_params,  # for mock data tests
        )


class MCMCResult:
    """Container for MCMC chain results with analysis methods."""

    def __init__(self, chain, chi2_chain, accepted, param_names,
                 n_burnin, thin, accept_rate, elapsed, true_params=None):
        self.full_chain = chain
        self.full_chi2 = chi2_chain
        self.accepted = accepted
        self.param_names = param_names
        self.n_burnin = n_burnin
        self.thin = thin
        self.accept_rate = accept_rate
        self.elapsed = elapsed
        self.true_params = true_params

        # Post-burn-in, thinned chain
        self.chain = chain[n_burnin::thin]
        self.chi2 = chi2_chain[n_burnin::thin]
        self.n_samples = len(self.chain)

    def summary(self):
        """Print parameter summary: mean, std, 68% CI."""
        print(f"\n{'='*70}")
        print(f"MCMC Summary ({self.n_samples} post-burn-in samples)")
        print(f"{'='*70}")
        print(f"{'Parameter':<12} {'Mean':>12} {'Std':>12} "
              f"{'68% CI':>24} {'True':>12}")
        print(f"{'-'*70}")

        for i, name in enumerate(self.param_names):
            samples = self.chain[:, i]
            mean = np.mean(samples)
            std = np.std(samples)
            lo = np.percentile(samples, 16)
            hi = np.percentile(samples, 84)

            true_str = ""
            if self.true_params and name in self.true_params:
                true_val = self.true_params[name]
                bias_sigma = (mean - true_val) / std if std > 0 else 0
                true_str = f"{true_val:.6g} ({bias_sigma:+.1f}s)"

            # Format based on parameter scale
            if name == 'A_s':
                print(f"{name:<12} {mean:.4e} {std:.4e} "
                      f"[{lo:.4e}, {hi:.4e}] {true_str:>12}")
            else:
                print(f"{name:<12} {mean:.6f} {std:.6f} "
                      f"[{lo:.6f}, {hi:.6f}] {true_str:>12}")

        best_idx = np.argmin(self.chi2)
        print(f"\nBest-fit chi^2 = {self.chi2[best_idx]:.1f}")
        print(f"Best-fit params: ", end="")
        for i, name in enumerate(self.param_names):
            val = self.chain[best_idx, i]
            if name == 'A_s':
                print(f"{name}={val:.4e}", end="  ")
            else:
                print(f"{name}={val:.6f}", end="  ")
        print()

    def best_fit(self):
        """Return best-fit parameter dict."""
        best_idx = np.argmin(self.chi2)
        return {name: self.chain[best_idx, i]
                for i, name in enumerate(self.param_names)}

    def convergence_check(self):
        """
        Basic convergence diagnostics.

        Returns True if chain appears converged (Geweke test + trace stability).
        """
        print(f"\n--- Convergence Diagnostics ---")

        converged = True
        for i, name in enumerate(self.param_names):
            samples = self.chain[:, i]
            n = len(samples)

            # Geweke test: compare first 10% and last 50%
            first = samples[:max(1, n // 10)]
            last = samples[max(1, n // 2):]
            z_score = (np.mean(first) - np.mean(last)) / \
                np.sqrt(np.var(first) / len(first) + np.var(last) / len(last) + 1e-30)

            status = "OK" if abs(z_score) < 2.0 else "WARN"
            if abs(z_score) >= 2.0:
                converged = False
            print(f"  {name:<12} Geweke z = {z_score:+.2f}  [{status}]")

        # Chi^2 stability: is the chain exploring the same region?
        chi2_first_half = np.mean(self.chi2[:len(self.chi2)//2])
        chi2_second_half = np.mean(self.chi2[len(self.chi2)//2:])
        chi2_ratio = chi2_second_half / (chi2_first_half + 1e-30)
        chi2_status = "OK" if 0.8 < chi2_ratio < 1.2 else "WARN"
        print(f"  chi^2 ratio (2nd/1st half): {chi2_ratio:.3f}  [{chi2_status}]")

        if converged:
            print("  CONVERGED (all Geweke |z| < 2)")
        else:
            print("  NOT YET CONVERGED — increase n_steps or tune step sizes")

        return converged


# ============================================================================
# Component 4: Mock data generator
# ============================================================================

class MockDataGenerator:
    """
    Generate mock Planck-like C_l data for testing MCMC.

    Creates mock data by running mlx_class at known ("true") parameters
    and adding realistic noise based on Planck-like error bars.
    """

    def __init__(self, evaluator):
        self.evaluator = evaluator

    def generate(self, true_params=None, noise_level=1.0, seed=12345):
        """
        Generate mock D_l data with Planck-like noise.

        Parameters
        ----------
        true_params : dict or None
            True cosmological parameters. None = use fiducial.
        noise_level : float
            Scale factor for noise (1.0 = Planck-like).
        seed : int
            Random seed for reproducibility.

        Returns
        -------
        ell_values : array
        data_Dl : array (muK^2)
        errors_Dl : array (muK^2)
        true_params : dict
        """
        if true_params is None:
            true_params = dict(FIDUCIAL_PARAMS)

        print(f"[MockData] Generating mock data at true parameters:")
        for k, v in true_params.items():
            if k == 'A_s':
                print(f"  {k} = {v:.4e}")
            else:
                print(f"  {k} = {v}")

        # Compute true C_l
        ell_values, Dl_true = self.evaluator.solver.compute_cl(true_params)

        # Planck-like errors: cosmic variance + instrument noise
        # sigma_Dl ~ Dl * sqrt(2/(2l+1)) * f_sky^{-0.5} + noise
        f_sky = 0.70  # Planck effective sky fraction
        ell_f = ell_values.astype(float)

        # Cosmic variance
        sigma_cv = Dl_true * np.sqrt(2.0 / (2.0 * ell_f + 1.0)) / np.sqrt(f_sky)

        # Instrument noise (Planck-like: ~40 uK-arcmin for 143 GHz)
        theta_beam = 7.0 / 60.0 * np.pi / 180.0  # 7 arcmin beam in radians
        sigma_noise_uK_arcmin = 40.0  # muK-arcmin
        N_l = (sigma_noise_uK_arcmin * np.pi / 180.0 / 60.0)**2 * \
            np.exp(ell_f * (ell_f + 1) * theta_beam**2 / (8.0 * np.log(2)))
        sigma_noise = ell_f * (ell_f + 1) / (2 * np.pi) * N_l

        # Total error
        errors_Dl = noise_level * np.sqrt(sigma_cv**2 + sigma_noise**2)

        # Minimum error floor (prevent division by zero at low ell)
        errors_Dl = np.maximum(errors_Dl, 1.0)

        # Add noise to true signal
        rng = np.random.RandomState(seed)
        data_Dl = Dl_true + rng.normal(0, errors_Dl)
        data_Dl = np.maximum(data_Dl, 0.0)  # D_l >= 0

        print(f"[MockData] Generated {len(ell_values)} data points")
        print(f"[MockData] D_l range: [{np.min(Dl_true):.0f}, {np.max(Dl_true):.0f}] muK^2")
        print(f"[MockData] Noise range: [{np.min(errors_Dl):.0f}, {np.max(errors_Dl):.0f}] muK^2")
        print(f"[MockData] SNR range: [{np.min(Dl_true/(errors_Dl+1e-10)):.1f}, "
              f"{np.max(Dl_true/(errors_Dl+1e-10)):.1f}]")

        return ell_values, data_Dl, errors_Dl, true_params


# ============================================================================
# Component 5: Diagnostics and plotting
# ============================================================================

class MCMCDiagnostics:
    """Plotting and diagnostics for MCMC chains."""

    @staticmethod
    def plot_chains(result, save_path=None):
        """Plot MCMC trace plots for all parameters."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            n_params = len(result.param_names)
            fig, axes = plt.subplots(n_params + 1, 1,
                                      figsize=(12, 3 * (n_params + 1)),
                                      sharex=True)

            steps = np.arange(len(result.full_chain))

            for i, name in enumerate(result.param_names):
                ax = axes[i]
                ax.plot(steps, result.full_chain[:, i], 'b-', lw=0.3, alpha=0.7)
                ax.axvline(result.n_burnin, color='red', ls='--', alpha=0.5,
                           label='burn-in')

                if result.true_params and name in result.true_params:
                    ax.axhline(result.true_params[name], color='green',
                               ls='-', lw=1.5, label='true')

                mean_val = np.mean(result.chain[:, i])
                ax.axhline(mean_val, color='orange', ls='--', lw=1, label='mean')

                ax.set_ylabel(name, fontsize=11)
                if i == 0:
                    ax.legend(fontsize=9, loc='upper right')

            # Chi^2 trace
            ax = axes[-1]
            ax.plot(steps, result.full_chi2, 'k-', lw=0.3, alpha=0.7)
            ax.axvline(result.n_burnin, color='red', ls='--', alpha=0.5)
            ax.set_ylabel(r'$\chi^2$', fontsize=11)
            ax.set_xlabel('Step', fontsize=11)

            plt.suptitle(f'MCMC Trace (accept={result.accept_rate:.1%}, '
                         f'{result.elapsed:.0f}s)', fontsize=13)
            plt.tight_layout()

            if save_path is None:
                save_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    'mcmc_trace.png'
                )
            plt.savefig(save_path, dpi=150)
            print(f"Trace plot: {save_path}")
            plt.close()

        except ImportError:
            print("matplotlib not available for trace plots")

    @staticmethod
    def plot_corner(result, save_path=None):
        """Plot 2D marginalized posteriors (corner plot)."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            n_params = len(result.param_names)
            fig, axes = plt.subplots(n_params, n_params,
                                      figsize=(3 * n_params, 3 * n_params))

            for i in range(n_params):
                for j in range(n_params):
                    ax = axes[i, j]
                    if j > i:
                        ax.axis('off')
                        continue

                    if i == j:
                        # 1D histogram
                        ax.hist(result.chain[:, i], bins=40,
                                color='steelblue', alpha=0.7,
                                density=True)
                        if result.true_params:
                            name = result.param_names[i]
                            if name in result.true_params:
                                ax.axvline(result.true_params[name],
                                           color='red', ls='--', lw=1.5)
                    else:
                        # 2D scatter
                        ax.scatter(result.chain[:, j], result.chain[:, i],
                                   s=1, alpha=0.2, c='steelblue')
                        if result.true_params:
                            name_x = result.param_names[j]
                            name_y = result.param_names[i]
                            if name_x in result.true_params and name_y in result.true_params:
                                ax.axvline(result.true_params[name_x],
                                           color='red', ls='--', lw=0.8, alpha=0.5)
                                ax.axhline(result.true_params[name_y],
                                           color='red', ls='--', lw=0.8, alpha=0.5)

                    if i == n_params - 1:
                        ax.set_xlabel(result.param_names[j], fontsize=9)
                    if j == 0:
                        ax.set_ylabel(result.param_names[i], fontsize=9)
                    ax.tick_params(labelsize=7)

            plt.suptitle('MCMC Posterior', fontsize=13)
            plt.tight_layout()

            if save_path is None:
                save_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    'mcmc_corner.png'
                )
            plt.savefig(save_path, dpi=120)
            print(f"Corner plot: {save_path}")
            plt.close()

        except ImportError:
            print("matplotlib not available for corner plot")

    @staticmethod
    def plot_bestfit(result, evaluator, data_Dl, errors_Dl, ell_data,
                     save_path=None):
        """Plot best-fit C_l vs data."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import io
            import contextlib

            # Merge best-fit with fiducial to get full parameter set
            best_params = dict(FIDUCIAL_PARAMS)
            best_params.update(result.best_fit())
            with contextlib.redirect_stdout(io.StringIO()):
                _, Dl_bestfit = evaluator.solver.compute_cl(best_params)

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8),
                                            gridspec_kw={'height_ratios': [3, 1]})

            # Top: D_l
            ax1.errorbar(ell_data, data_Dl, yerr=errors_Dl,
                         fmt='.', color='gray', alpha=0.5, ms=2,
                         label='Mock data')
            ax1.plot(evaluator.ell_values, Dl_bestfit, 'r-', lw=1.5,
                     label='Best-fit')

            if result.true_params:
                true_full = dict(FIDUCIAL_PARAMS)
                true_full.update(result.true_params)
                with contextlib.redirect_stdout(io.StringIO()):
                    _, Dl_true = evaluator.solver.compute_cl(true_full)
                ax1.plot(evaluator.ell_values, Dl_true, 'b--', lw=1,
                         alpha=0.7, label='True')

            ax1.set_ylabel(r'$D_\ell$ [$\mu K^2$]', fontsize=12)
            ax1.legend(fontsize=11)
            ax1.set_xlim(2, int(evaluator.ell_values[-1]))
            ax1.set_title('MCMC Best-fit vs Data', fontsize=13)

            # Bottom: residuals
            # Interpolate best-fit to data ells
            from scipy.interpolate import interp1d
            f_bf = interp1d(evaluator.ell_values, Dl_bestfit,
                            kind='linear', fill_value='extrapolate')
            Dl_bf_at_data = f_bf(ell_data)
            residual = (data_Dl - Dl_bf_at_data) / (errors_Dl + 1e-10)

            ax2.scatter(ell_data, residual, s=2, alpha=0.5, c='steelblue')
            ax2.axhline(0, color='k', ls='-', lw=0.5)
            ax2.axhline(1, color='gray', ls='--', lw=0.5)
            ax2.axhline(-1, color='gray', ls='--', lw=0.5)
            ax2.set_ylabel(r'Residual [$\sigma$]', fontsize=11)
            ax2.set_xlabel(r'$\ell$', fontsize=12)
            ax2.set_xlim(2, int(evaluator.ell_values[-1]))
            ax2.set_ylim(-5, 5)

            plt.tight_layout()

            if save_path is None:
                save_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    'mcmc_bestfit.png'
                )
            plt.savefig(save_path, dpi=150)
            print(f"Best-fit plot: {save_path}")
            plt.close()

        except ImportError:
            print("matplotlib not available for best-fit plot")


# ============================================================================
# Main: demo / test
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='mlx_class MCMC parameter fitting')
    parser.add_argument('--n_steps', type=int, default=200,
                        help='Number of MCMC steps (default: 200 for demo)')
    parser.add_argument('--n_params', type=int, default=3,
                        help='Number of parameters to vary (3=h,omega_b,omega_cdm; '
                             '6=all LCDM)')
    parser.add_argument('--mock', action='store_true',
                        help='Generate mock data and run MCMC recovery test')
    parser.add_argument('--batch_test', action='store_true',
                        help='Run batch evaluation speed test')
    parser.add_argument('--plot', action='store_true',
                        help='Generate diagnostic plots')
    parser.add_argument('--N_k', type=int, default=200,
                        help='k-modes per evaluation (default: 200)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    print("=" * 60)
    print("mlx_class MCMC Parameter Fitting")
    print(f"GPU: {mx.default_device()}")
    print("=" * 60)

    # ---- Batch evaluation speed test ----
    if args.batch_test:
        print("\n--- Batch Evaluation Speed Test ---")
        evaluator = BatchEvaluator(N_k=args.N_k, mode='fast', verbose=False)

        # Create 5 parameter sets with small variations
        params_list = []
        rng = np.random.RandomState(args.seed)
        for i in range(5):
            p = dict(FIDUCIAL_PARAMS)
            p['h'] += rng.normal(0, 0.005)
            p['omega_b'] += rng.normal(0, 0.0003)
            p['omega_cdm'] += rng.normal(0, 0.002)
            params_list.append(p)

        results = evaluator.evaluate(params_list)
        print(f"\nResults: {len(results)} spectra computed")
        for i, Dl in enumerate(results):
            print(f"  Set {i}: max(D_l) = {np.max(Dl):.0f} muK^2")
        return

    # ---- Mock data MCMC test ----
    if args.mock:
        print("\n--- Mock Data MCMC Recovery Test ---")
        t_total = time.time()

        # Step 1: Setup evaluator
        evaluator = BatchEvaluator(N_k=args.N_k, mode='fast', verbose=False)

        # Step 2: Generate mock data
        print("\n--- Generating mock data ---")
        mock_gen = MockDataGenerator(evaluator)
        ell_data, data_Dl, errors_Dl, true_params = mock_gen.generate(
            noise_level=1.0, seed=args.seed
        )

        # Step 3: Select parameters to vary
        if args.n_params == 3:
            vary_names = ['h', 'omega_b', 'omega_cdm']
        elif args.n_params == 4:
            vary_names = ['h', 'omega_b', 'omega_cdm', 'n_s']
        else:
            vary_names = ['h', 'omega_b', 'omega_cdm', 'A_s', 'n_s', 'tau_reio']

        # Initial parameters: offset from true values to test convergence
        rng = np.random.RandomState(args.seed + 1)
        initial_params = dict(FIDUCIAL_PARAMS)  # start at fiducial
        for name in vary_names:
            # Offset by ~2 sigma from true
            offset = 2.0 * DEFAULT_STEP_SIZES[name] * rng.choice([-1, 1])
            initial_params[name] = true_params[name] + offset

        # Fix non-varied parameters at true values
        for name in FIDUCIAL_PARAMS:
            if name not in vary_names:
                initial_params[name] = true_params[name]

        print(f"\nVarying: {vary_names}")
        print(f"Starting point offsets from true:")
        for name in vary_names:
            diff = initial_params[name] - true_params[name]
            sigma = DEFAULT_STEP_SIZES[name]
            print(f"  {name}: {diff/sigma:+.1f} sigma")

        # Step 4: Run MCMC
        print(f"\n--- Running MCMC ({args.n_steps} steps) ---")
        mcmc = MetropolisMCMC(
            data_Dl=data_Dl,
            errors_Dl=errors_Dl,
            ell_values=ell_data,
            evaluator=evaluator,
        )

        # Only include varied parameters in the MCMC
        initial_varied = {k: initial_params[k] for k in vary_names}
        step_sizes_varied = {k: DEFAULT_STEP_SIZES[k] for k in vary_names}

        # Wrap evaluator to fix non-varied parameters
        class FixedParamEvaluator:
            """Wrapper that fixes some parameters at known values."""
            def __init__(self, base_evaluator, fixed_params):
                self.base = base_evaluator
                self.fixed = dict(fixed_params)
                self.ell_values = base_evaluator.ell_values
                self._cache_hits = 0
                self._cache_misses = 0

            def evaluate_single(self, varied_params):
                full_params = dict(self.fixed)
                full_params.update(varied_params)
                result = self.base.evaluate_single(full_params)
                self._cache_hits = self.base._cache_hits
                self._cache_misses = self.base._cache_misses
                return result

        fixed_params = {k: true_params[k] for k in FIDUCIAL_PARAMS
                        if k not in vary_names}
        wrapped_evaluator = FixedParamEvaluator(evaluator, fixed_params)

        mcmc_wrapped = MetropolisMCMC(
            data_Dl=data_Dl,
            errors_Dl=errors_Dl,
            ell_values=ell_data,
            evaluator=wrapped_evaluator,
        )

        result = mcmc_wrapped.run(
            initial_params=initial_varied,
            n_steps=args.n_steps,
            step_sizes=step_sizes_varied,
            seed=args.seed,
        )

        # Store true params for diagnostics
        result.true_params = {k: true_params[k] for k in vary_names}

        # Step 5: Results
        result.summary()
        result.convergence_check()

        # Step 6: Plots
        if args.plot:
            MCMCDiagnostics.plot_chains(result)
            if len(vary_names) <= 6:
                MCMCDiagnostics.plot_corner(result)
            MCMCDiagnostics.plot_bestfit(
                result, evaluator, data_Dl, errors_Dl, ell_data
            )

        # Save chain
        out_dir = os.path.dirname(os.path.abspath(__file__))
        chain_path = os.path.join(out_dir, 'mcmc_chain.npz')
        np.savez(chain_path,
                 chain=result.chain,
                 chi2=result.chi2,
                 param_names=result.param_names,
                 true_params=np.array([true_params[k] for k in result.param_names]))
        print(f"\nChain saved: {chain_path}")

        t_total_elapsed = time.time() - t_total
        print(f"\nTotal time: {t_total_elapsed:.0f}s")
        print(f"Equivalent CLASS time: ~{args.n_steps * 10:.0f}s "
              f"(speedup: ~{args.n_steps * 10 / t_total_elapsed:.0f}x)")
        return

    # Default: just run a single evaluation to verify the pipeline works
    print("\n--- Single evaluation test ---")
    solver = ParameterizedSolver(N_k=args.N_k, mode='fast')
    t0 = time.time()
    ell_out, Dl = solver.compute_cl(FIDUCIAL_PARAMS)
    t_eval = time.time() - t0
    print(f"\nSingle evaluation: {t_eval:.2f}s")
    print(f"ell range: [{ell_out[0]}, {ell_out[-1]}]")
    print(f"D_l range: [{np.min(Dl):.0f}, {np.max(Dl):.0f}] muK^2")
    print(f"\nTo run MCMC: python -m mlx_class.mcmc --mock --n_steps 200 --plot")


if __name__ == '__main__':
    main()
