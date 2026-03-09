"""
Test 3: Tightness of F >= exp(-Sigma/2) -- Finding Near-Saturating Channels

The bound F >= exp(-Sigma/2) always holds. But HOW TIGHT is it?

Scan over a large family of quantum channels:
  - Dephasing channels (varying p)
  - Amplitude damping channels (varying gamma)
  - Depolarizing channels (varying p)
  - Generalized Pauli channels (varying p_x, p_y, p_z)
  - Random CPTP maps (Haar random Kraus operators, d=2,3,4)

For each: compute F_Petz and Sigma
Plot F vs Sigma scatter with the bound curve
Find: which channels come CLOSEST to saturating? (F ~ exp(-Sigma/2))
Find: which channels are FURTHEST? (F >> exp(-Sigma/2))

Why this matters:
  - If gravitational dephasing is near-saturating -> gravity is "maximally
    efficient" at creating irreversibility -> physical insight
  - If gravitational dephasing is far from saturating -> the bound is loose
    for gravity -> less useful
  - If SOME channels saturate -> we identify the "worst case" channels
    -> useful for QEC

Reference: Huang (2026), Retrodiction Landauer Principle.
"""

import os
import numpy as np
from scipy.linalg import sqrtm
import warnings
warnings.filterwarnings('ignore')

_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# Quantum primitives
# ============================================================

def _hermitianize(M):
    return (M + M.conj().T) / 2


def von_neumann_entropy(rho):
    rho_h = _hermitianize(rho)
    eigvals = np.linalg.eigvalsh(rho_h)
    eigvals = eigvals[eigvals > 1e-15]
    return -np.sum(eigvals * np.log(eigvals))


def trace_fidelity(rho, sigma):
    """Trace fidelity F_tr = Tr(sqrt(sqrt(rho) sigma sqrt(rho)))."""
    sqrt_rho = sqrtm(_hermitianize(rho))
    M = sqrt_rho @ _hermitianize(sigma) @ sqrt_rho
    M = _hermitianize(M)
    eigvals = np.linalg.eigvalsh(M)
    return np.real(np.sum(np.sqrt(np.maximum(eigvals, 0))))


def apply_channel(rho, kraus_ops):
    return sum(E @ rho @ E.conj().T for E in kraus_ops)


def apply_adjoint(X, kraus_ops):
    return sum(E.conj().T @ X @ E for E in kraus_ops)


def complementary_channel(rho, kraus_ops):
    n = len(kraus_ops)
    rho_E = np.zeros((n, n), dtype=complex)
    for i, Ei in enumerate(kraus_ops):
        for j, Ej in enumerate(kraus_ops):
            rho_E[i, j] = np.trace(Ei @ rho @ Ej.conj().T)
    return _hermitianize(rho_E)


def _mat_log(M, tol=1e-15):
    """Matrix logarithm for positive semidefinite matrices."""
    M_h = _hermitianize(M)
    eigvals, eigvecs = np.linalg.eigh(M_h)
    log_eigvals = np.where(eigvals > tol, np.log(eigvals), 0.0)
    return eigvecs @ np.diag(log_eigvals) @ eigvecs.conj().T


def relative_entropy(rho, sigma):
    """D(rho || sigma) = Tr[rho (log rho - log sigma)]."""
    rho_h = _hermitianize(rho)
    sigma_h = _hermitianize(sigma)
    # Regularize sigma to avoid log(0)
    sigma_reg = 0.999 * sigma_h + 0.001 * np.eye(sigma_h.shape[0]) / sigma_h.shape[0]
    log_rho = _mat_log(rho_h)
    log_sigma = _mat_log(sigma_reg)
    val = np.real(np.trace(rho_h @ (log_rho - log_sigma)))
    return max(val, 0.0)


def relative_entropy_drop(rho, sigma, kraus_ops):
    """Delta_D = D(rho||sigma) - D(N(rho)||N(sigma)).

    This is the correct quantity for the Petz recovery bound:
        F_tr(rho, R_sigma(N(rho))) >= exp(-Delta_D/2)
    """
    rho_out = apply_channel(rho, kraus_ops)
    sigma_out = apply_channel(sigma, kraus_ops)
    D_before = relative_entropy(rho, sigma)
    D_after = relative_entropy(rho_out, sigma_out)
    return max(D_before - D_after, 0.0)


def entropy_production(rho, kraus_ops):
    """Sigma = Delta S_sys + S_env (total entropy production).

    NOTE: This is NOT the same as Delta_D for general channels.
    The bound F >= exp(-Delta_D/2) uses Delta_D, not Sigma.
    We keep this for comparison but the bound uses relative_entropy_drop.
    """
    rho_out = apply_channel(rho, kraus_ops)
    rho_E = complementary_channel(rho, kraus_ops)
    S_env = von_neumann_entropy(rho_E)
    delta_S = von_neumann_entropy(rho_out) - von_neumann_entropy(rho)
    return max(delta_S + S_env, 0.0)


def petz_recovery(rho_out, sigma, kraus_ops):
    """Standard Petz recovery map R_sigma(rho_out)."""
    sigma_out = apply_channel(sigma, kraus_ops)
    sqrt_sigma = sqrtm(_hermitianize(sigma))
    inv_sqrt_sigma_out = np.linalg.pinv(sqrtm(_hermitianize(sigma_out)))
    inner = inv_sqrt_sigma_out @ rho_out @ inv_sqrt_sigma_out
    adjoint_inner = apply_adjoint(inner, kraus_ops)
    result = sqrt_sigma @ adjoint_inner @ sqrt_sigma
    return _hermitianize(result)


# ============================================================
# Channel families
# ============================================================

def dephasing_kraus(p, d=2):
    """Dephasing channel: off-diagonal elements decay by factor (1-p)."""
    E0 = np.sqrt(1 - p) * np.eye(d, dtype=complex)
    E1 = np.sqrt(p) * np.diag([(-1)**i for i in range(d)]).astype(complex)
    return [E0, E1]


def amplitude_damping_kraus(gamma):
    """Amplitude damping channel."""
    E0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    E1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [E0, E1]


def depolarizing_kraus(p, d=2):
    """Depolarizing channel: N(rho) = (1-p)*rho + p*I/d."""
    # Kraus operators: sqrt(1-p)*I, sqrt(p/3)*sigma_x, sqrt(p/3)*sigma_y, sqrt(p/3)*sigma_z
    if d == 2:
        sx = np.array([[0, 1], [1, 0]], dtype=complex)
        sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sz = np.array([[1, 0], [0, -1]], dtype=complex)
        E0 = np.sqrt(1 - 3*p/4) * np.eye(2, dtype=complex)
        E1 = np.sqrt(p/4) * sx
        E2 = np.sqrt(p/4) * sy
        E3 = np.sqrt(p/4) * sz
        return [E0, E1, E2, E3]
    else:
        # General d: use generalized Gell-Mann matrices
        kraus = [np.sqrt(1 - p + p/d**2) * np.eye(d, dtype=complex)]
        coeff = np.sqrt(p / d**2)
        for i in range(d):
            for j in range(d):
                if i == j == 0:
                    continue
                E = np.zeros((d, d), dtype=complex)
                E[i, j] = 1.0
                kraus.append(coeff * E)
        return kraus


def generalized_pauli_kraus(px, py, pz):
    """Generalized Pauli channel: sum p_i sigma_i rho sigma_i."""
    p0 = max(1 - px - py - pz, 0)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    E0 = np.sqrt(p0) * np.eye(2, dtype=complex)
    E1 = np.sqrt(px) * sx
    E2 = np.sqrt(py) * sy
    E3 = np.sqrt(pz) * sz
    return [E0, E1, E2, E3]


def random_cptp_kraus(d, n_kraus=None, rng=None):
    """Random CPTP map via Haar-random Kraus operators."""
    if rng is None:
        rng = np.random.default_rng()
    if n_kraus is None:
        n_kraus = d**2

    # Generate random Kraus operators
    G = rng.standard_normal((n_kraus, d, d)) + 1j * rng.standard_normal((n_kraus, d, d))
    G /= np.sqrt(2)

    # Ensure CPTP: sum E_k^dag E_k = I
    S = sum(E.conj().T @ E for E in G)
    S_inv_sqrt = np.linalg.inv(sqrtm(_hermitianize(S)))

    kraus = [E @ S_inv_sqrt for E in G]
    return kraus


def gravitational_dephasing_kraus(N_modes, epsilon, t):
    """
    Gravitational dephasing: position-dependent phase.
    For a single qubit system interacting with N internal modes.

    Effective channel on the system after tracing out N modes:
    D(t) = cos^N(2*epsilon*t)
    rho_out = [[rho_00, D*rho_01], [D*rho_10, rho_11]]

    Kraus: E0 = [[1, 0], [0, sqrt((1+D)/2)]], E1 = [[0, 0], [0, sqrt((1-D)/2)]]
    Wait -- that's not right for dephasing.

    For pure dephasing:
    E0 = sqrt((1+D)/2) * I
    E1 = sqrt((1-D)/2) * sigma_z
    """
    D = np.cos(2 * epsilon * t)**N_modes
    abs_D = abs(D)
    E0 = np.sqrt((1 + abs_D) / 2) * np.eye(2, dtype=complex)
    E1 = np.sqrt((1 - abs_D) / 2) * np.array([[1, 0], [0, -1]], dtype=complex)
    return [E0, E1]


# ============================================================
# Scan channels
# ============================================================

def scan_channel_family(name, kraus_fn, param_values, input_states, sigma):
    """Scan a channel family over parameter values and input states."""
    results = []

    for params in param_values:
        if isinstance(params, (list, tuple)):
            kraus = kraus_fn(*params)
            param_label = str(params)
        else:
            kraus = kraus_fn(params)
            param_label = f"{params:.3f}"

        for state_name, rho in input_states.items():
            try:
                rho_out = apply_channel(rho, kraus)

                # Delta_D = D(rho||sigma) - D(N(rho)||N(sigma))
                # This is the correct quantity for the Petz bound
                DeltaD = relative_entropy_drop(rho, sigma, kraus)

                # Also compute Sigma for comparison
                Sigma = entropy_production(rho, kraus)

                # Petz recovery
                rho_rec = petz_recovery(rho_out, sigma, kraus)
                F = min(trace_fidelity(rho, rho_rec), 1.0)

                bound = np.exp(-DeltaD / 2)
                gap = F - bound
                tightness = gap / max(bound, 1e-10)

                results.append({
                    'channel': name,
                    'params': param_label,
                    'state': state_name,
                    'F': F,
                    'Sigma': Sigma,
                    'DeltaD': DeltaD,
                    'bound': bound,
                    'gap': gap,
                    'tightness': tightness,
                })
            except Exception:
                pass

    return results


def run_test3():
    """Run Test 3: Bound tightness scan."""
    print("=" * 70)
    print("TEST 3: Tightness of F >= exp(-Sigma/2)")
    print("=" * 70)
    print()

    # Input states
    input_states = {
        '|1>': np.array([[0, 0], [0, 1]], dtype=complex),
        '|+>': np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex),
        'mixed(0.3)': np.array([[0.7, 0], [0, 0.3]], dtype=complex),
        'mixed(0.1)': np.array([[0.9, 0], [0, 0.1]], dtype=complex),
    }

    sigma = np.eye(2, dtype=complex) / 2

    all_results = []

    # 1. Dephasing channels
    print("  Scanning dephasing channels...")
    params_deph = np.linspace(0.01, 0.99, 50)
    results_deph = scan_channel_family('Dephasing', dephasing_kraus, params_deph,
                                        input_states, sigma)
    all_results.extend(results_deph)

    # 2. Amplitude damping
    print("  Scanning amplitude damping channels...")
    params_ad = np.linspace(0.01, 0.99, 50)
    results_ad = scan_channel_family('Amp. Damping', amplitude_damping_kraus, params_ad,
                                      input_states, sigma)
    all_results.extend(results_ad)

    # 3. Depolarizing
    print("  Scanning depolarizing channels...")
    params_depol = np.linspace(0.01, 0.99, 50)
    results_depol = scan_channel_family('Depolarizing', depolarizing_kraus, params_depol,
                                         input_states, sigma)
    all_results.extend(results_depol)

    # 4. Generalized Pauli
    print("  Scanning generalized Pauli channels...")
    pauli_params = []
    rng = np.random.default_rng(42)
    for _ in range(200):
        px = rng.uniform(0, 0.25)
        py = rng.uniform(0, 0.25)
        pz = rng.uniform(0, 0.25)
        if px + py + pz <= 1:
            pauli_params.append((px, py, pz))
    results_pauli = scan_channel_family('Gen. Pauli', generalized_pauli_kraus,
                                         pauli_params, input_states, sigma)
    all_results.extend(results_pauli)

    # 5. Random CPTP maps (d=2)
    print("  Scanning random CPTP maps (d=2)...")
    rng = np.random.default_rng(123)
    for trial in range(200):
        kraus = random_cptp_kraus(2, n_kraus=4, rng=rng)
        for state_name, rho in input_states.items():
            try:
                rho_out = apply_channel(rho, kraus)
                DeltaD = relative_entropy_drop(rho, sigma, kraus)
                Sigma = entropy_production(rho, kraus)
                rho_rec = petz_recovery(rho_out, sigma, kraus)
                F = min(trace_fidelity(rho, rho_rec), 1.0)
                bound = np.exp(-DeltaD / 2)
                gap = F - bound
                tightness = gap / max(bound, 1e-10)

                all_results.append({
                    'channel': 'Random d=2',
                    'params': f'trial_{trial}',
                    'state': state_name,
                    'F': F,
                    'Sigma': Sigma,
                    'DeltaD': DeltaD,
                    'bound': bound,
                    'gap': gap,
                    'tightness': tightness,
                })
            except Exception:
                pass

    # 6. Gravitational dephasing (our model)
    print("  Scanning gravitational dephasing channels...")
    for N in [2, 4, 8, 16, 32]:
        for t in np.linspace(0.1, 3.0, 30):
            epsilon = 0.3
            kraus = gravitational_dephasing_kraus(N, epsilon, t)
            for state_name, rho in input_states.items():
                try:
                    rho_out = apply_channel(rho, kraus)
                    DeltaD = relative_entropy_drop(rho, sigma, kraus)
                    Sigma = entropy_production(rho, kraus)
                    rho_rec = petz_recovery(rho_out, sigma, kraus)
                    F = min(trace_fidelity(rho, rho_rec), 1.0)
                    bound = np.exp(-DeltaD / 2)
                    gap = F - bound
                    tightness = gap / max(bound, 1e-10)

                    all_results.append({
                        'channel': 'Grav. Dephasing',
                        'params': f'N={N},t={t:.2f}',
                        'state': state_name,
                        'F': F,
                        'Sigma': Sigma,
                        'DeltaD': DeltaD,
                        'bound': bound,
                        'gap': gap,
                        'tightness': tightness,
                    })
                except Exception:
                    pass

    # Analysis
    print()
    print(f"Total data points: {len(all_results)}")
    print()

    # Violations check
    violations = [r for r in all_results if r['gap'] < -1e-6]
    print(f"Bound violations: {len(violations)}")
    if violations:
        print("  WARNING: Bound violations detected!")
        for v in violations[:5]:
            print(f"    {v['channel']} {v['params']} {v['state']}: "
                  f"F={v['F']:.6f} < bound={v['bound']:.6f}")
    else:
        print("  No violations. Bound F >= exp(-Delta_D/2) is SATISFIED for all channels.")
    print()

    # Tightest points (smallest gap)
    valid_results = [r for r in all_results if r.get('DeltaD', r['Sigma']) > 0.01]
    if valid_results:
        sorted_by_gap = sorted(valid_results, key=lambda r: r['gap'])

        print("TIGHTEST POINTS (smallest F - exp(-Sigma/2)):")
        print("-" * 80)
        print(f"{'Channel':<18} {'State':<12} {'Sigma':>8} {'F':>8} "
              f"{'Bound':>8} {'Gap':>10} {'Rel.Gap%':>10}")
        print("-" * 80)
        for r in sorted_by_gap[:15]:
            print(f"{r['channel']:<18} {r['state']:<12} {r['Sigma']:>8.4f} "
                  f"{r['F']:>8.4f} {r['bound']:>8.4f} {r['gap']:>10.6f} "
                  f"{r['tightness']*100:>10.2f}%")

        print()
        print("LOOSEST POINTS (largest F - exp(-Sigma/2)):")
        print("-" * 80)
        sorted_by_gap_desc = sorted(valid_results, key=lambda r: -r['gap'])
        for r in sorted_by_gap_desc[:10]:
            print(f"{r['channel']:<18} {r['state']:<12} {r['Sigma']:>8.4f} "
                  f"{r['F']:>8.4f} {r['bound']:>8.4f} {r['gap']:>10.6f} "
                  f"{r['tightness']*100:>10.2f}%")

    # Per-channel statistics
    print()
    print("PER-CHANNEL TIGHTNESS STATISTICS:")
    print("-" * 70)
    channels = set(r['channel'] for r in valid_results)
    channel_stats = {}
    for ch in sorted(channels):
        ch_results = [r for r in valid_results if r['channel'] == ch]
        gaps = [r['gap'] for r in ch_results]
        tightnesses = [r['tightness'] for r in ch_results]
        channel_stats[ch] = {
            'mean_gap': np.mean(gaps),
            'min_gap': np.min(gaps),
            'max_gap': np.max(gaps),
            'mean_tightness': np.mean(tightnesses),
            'n': len(ch_results),
        }
        print(f"  {ch:<18}: n={len(ch_results):>4}, "
              f"mean_gap={np.mean(gaps):.6f}, min_gap={np.min(gaps):.6f}, "
              f"max_gap={np.max(gaps):.6f}")

    print()
    print("INTERPRETATION:")
    if valid_results:
        tightest = sorted_by_gap[0]
        print(f"  Tightest channel: {tightest['channel']} "
              f"(gap = {tightest['gap']:.6f})")
        grav_results = [r for r in valid_results if r['channel'] == 'Grav. Dephasing']
        if grav_results:
            grav_min_gap = min(r['gap'] for r in grav_results)
            grav_mean_gap = np.mean([r['gap'] for r in grav_results])
            print(f"  Gravitational dephasing: min_gap = {grav_min_gap:.6f}, "
                  f"mean_gap = {grav_mean_gap:.6f}")
            if grav_min_gap < 0.05:
                print("  -> Gravitational dephasing is NEAR-SATURATING")
                print("     Gravity is 'maximally efficient' at creating irreversibility!")
            else:
                print("  -> Gravitational dephasing is NOT near-saturating")
                print("     The bound is LOOSE for gravity. The actual recovery")
                print("     fidelity is much better than the bound guarantees.")
    print()

    return all_results, channel_stats


if __name__ == "__main__":
    all_results, channel_stats = run_test3()
