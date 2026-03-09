"""
Supplement verification: F >= exp(-Sigma/2) for quantum channels
studied by Singh et al. (2025) NMR experiments.

Computes EXACT analytical values for:
  (1) Amplitude Damping (AD) channel
  (2) Phase Damping (PD) channel

For each channel, verifies the JRSWW bound:
    F(rho, R_Petz(N(rho))) >= exp(-Sigma/2)

where Sigma = D(rho||sigma) - D(N(rho)||N(sigma)) is the entropy production
(information loss through the channel relative to the reference state).

For the dephasing channel, also compares:
    tau_Petz vs tau_optimal vs tau_JRSWW

Reference: Huang (2026), "Petz Recovery Map as Retrodiction Functor".
Singh et al. (2025), NMR implementation of Petz recovery.
"""

import os
import csv
import numpy as np
from scipy.linalg import logm, sqrtm

_DIR = os.path.dirname(os.path.abspath(__file__))

# =============================================================================
# Quantum information utilities
# =============================================================================

def ensure_hermitian(M):
    """Force a matrix to be Hermitian."""
    return (M + M.conj().T) / 2


def ensure_density_matrix(rho):
    """Project onto valid density matrix (Hermitian, PSD, trace 1)."""
    rho = ensure_hermitian(rho)
    eigvals, eigvecs = np.linalg.eigh(rho)
    eigvals = np.maximum(eigvals, 0)
    rho = eigvecs @ np.diag(eigvals) @ eigvecs.conj().T
    tr = np.real(np.trace(rho))
    if tr > 1e-15:
        rho = rho / tr
    return rho


def matrix_sqrt_psd(M):
    """Compute PSD matrix square root via eigendecomposition."""
    M_h = ensure_hermitian(M)
    eigvals, eigvecs = np.linalg.eigh(M_h)
    eigvals = np.maximum(eigvals, 0)
    return eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.conj().T


def matrix_inv_sqrt(M, eps=1e-12):
    """Compute M^{-1/2} with regularization on support."""
    M_h = ensure_hermitian(M)
    eigvals, eigvecs = np.linalg.eigh(M_h)
    inv_sqrt_vals = np.zeros_like(eigvals)
    mask = eigvals > eps
    inv_sqrt_vals[mask] = 1.0 / np.sqrt(eigvals[mask])
    return eigvecs @ np.diag(inv_sqrt_vals) @ eigvecs.conj().T


def uhlmann_fidelity(rho, sigma):
    """
    Compute Uhlmann fidelity F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2.

    For pure state rho = |psi><psi|, this simplifies to F = <psi|sigma|psi>.
    """
    rho = ensure_density_matrix(rho)
    sigma = ensure_density_matrix(sigma)
    sqrt_rho = matrix_sqrt_psd(rho)
    product = sqrt_rho @ sigma @ sqrt_rho
    product = ensure_hermitian(product)
    eigvals = np.linalg.eigvalsh(product)
    eigvals = np.maximum(eigvals, 0)
    F = np.real(np.sum(np.sqrt(eigvals)))**2
    return min(F, 1.0)


def umegaki_relative_entropy(rho, sigma, eps=1e-15):
    """
    Compute D(rho||sigma) = Tr[rho (log rho - log sigma)].

    Returns np.inf if supp(rho) is not contained in supp(sigma).
    """
    rho = ensure_density_matrix(rho)
    sigma = ensure_density_matrix(sigma)

    # Get eigendecompositions
    eig_rho, vec_rho = np.linalg.eigh(rho)
    eig_sigma, vec_sigma = np.linalg.eigh(sigma)

    # Check support condition: if rho has support where sigma doesn't, D = inf
    for i, lr in enumerate(eig_rho):
        if lr > eps:
            # Check if this eigenvector of rho has overlap with null space of sigma
            v = vec_rho[:, i]
            sigma_v = sigma @ v
            overlap = np.real(np.vdot(v, sigma_v))
            if overlap < eps:
                return np.inf

    # Compute via eigendecompositions for numerical stability
    # D = sum_i lambda_i (log lambda_i - <v_i|log sigma|v_i>)
    log_sigma = np.zeros_like(sigma, dtype=complex)
    for i, ls in enumerate(eig_sigma):
        if ls > eps:
            log_sigma += np.log(ls) * np.outer(vec_sigma[:, i], vec_sigma[:, i].conj())
        else:
            # Assign large negative value for zero eigenvalues
            log_sigma += (-1e15) * np.outer(vec_sigma[:, i], vec_sigma[:, i].conj())

    D = 0.0
    for i, lr in enumerate(eig_rho):
        if lr > eps:
            v = vec_rho[:, i]
            log_sigma_elem = np.real(np.vdot(v, log_sigma @ v))
            D += lr * (np.log(lr) - log_sigma_elem)

    return np.real(D)


def apply_channel(kraus_ops, rho):
    """Apply quantum channel N(rho) = sum_k K_k rho K_k^dag."""
    rho = ensure_density_matrix(rho)
    out = np.zeros_like(rho, dtype=complex)
    for K in kraus_ops:
        out += K @ rho @ K.conj().T
    return ensure_density_matrix(out)


# =============================================================================
# Channel definitions
# =============================================================================

def ad_kraus(p):
    """Amplitude damping Kraus operators."""
    K0 = np.array([[1, 0], [0, np.sqrt(1 - p)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(p)], [0, 0]], dtype=complex)
    return [K0, K1]


def pd_kraus(p):
    """Phase damping (dephasing) Kraus operators."""
    K0 = np.sqrt(1 - p / 2) * np.eye(2, dtype=complex)
    K1 = np.sqrt(p / 2) * np.array([[1, 0], [0, -1]], dtype=complex)
    return [K0, K1]


# =============================================================================
# Petz recovery map constructions
# =============================================================================

def petz_recovery_ad(rho, p, eps_ref):
    """
    Apply Petz recovery map for amplitude damping channel.

    Reference state: sigma = (1-eps)|0><0| + eps|1><1|

    Petz recovery Kraus operators for AD:
        M0 = [[sqrt((1-eps)/(1-(1-p)*eps)), 0], [0, 1]]
        M1 = [[0, 0], [sqrt(p*eps/(1-(1-p)*eps)), 0]]
    """
    # Channel output state
    N_sigma_00 = 1 - (1 - p) * eps_ref  # (1-eps) + p*eps = 1 - eps + p*eps

    if N_sigma_00 < 1e-15 or N_sigma_00 > 1 - 1e-15:
        # Degenerate case
        return rho

    M0 = np.array([
        [np.sqrt((1 - eps_ref) / N_sigma_00), 0],
        [0, 1]
    ], dtype=complex)

    M1 = np.array([
        [0, 0],
        [np.sqrt(p * eps_ref / N_sigma_00), 0]
    ], dtype=complex)

    return apply_channel([M0, M1], rho)


def petz_recovery_general(rho_in, kraus_ops, sigma):
    """
    General Petz recovery map:
        R_{sigma,N}(Y) = sigma^{1/2} N^dag(N(sigma)^{-1/2} Y N(sigma)^{-1/2}) sigma^{1/2}

    Applied to Y = N(rho_in).
    """
    sigma = ensure_density_matrix(sigma)
    N_rho = apply_channel(kraus_ops, rho_in)
    N_sigma = apply_channel(kraus_ops, sigma)

    sqrt_sigma = matrix_sqrt_psd(sigma)
    inv_sqrt_N_sigma = matrix_inv_sqrt(N_sigma)

    # Compute N^dag(inv_sqrt_N_sigma @ N_rho @ inv_sqrt_N_sigma)
    inner = inv_sqrt_N_sigma @ N_rho @ inv_sqrt_N_sigma

    # N^dag(Y) = sum_k K_k^dag Y K_k
    adjoint_result = np.zeros((2, 2), dtype=complex)
    for K in kraus_ops:
        adjoint_result += K.conj().T @ inner @ K

    recovered = sqrt_sigma @ adjoint_result @ sqrt_sigma
    return ensure_density_matrix(recovered)


# =============================================================================
# Input states
# =============================================================================

ket0 = np.array([1, 0], dtype=complex)
ket1 = np.array([0, 1], dtype=complex)
ket_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
ket_minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
# Singh state: 0.9268|0> + 0.3754i|1>  (normalized)
ket_singh_raw = np.array([0.9268, 0.3754j], dtype=complex)
ket_singh = ket_singh_raw / np.linalg.norm(ket_singh_raw)

INPUT_STATES = {
    "|0>": ket0,
    "|1>": ket1,
    "|+>": ket_plus,
    "Singh": ket_singh,
}


def ket_to_dm(ket):
    """Convert ket vector to density matrix."""
    return np.outer(ket, ket.conj())


# =============================================================================
# Main verification
# =============================================================================

def verify_channel(channel_name, kraus_fn, sigma_fn, eps_values, p_values,
                   petz_fn=None):
    """
    Verify F >= exp(-Sigma/2) for a given channel across all parameters.

    Parameters
    ----------
    channel_name : str
    kraus_fn : callable(p) -> list of Kraus operators
    sigma_fn : callable(eps) -> density matrix (reference state)
    eps_values : list of float
    p_values : array of float
    petz_fn : callable(rho, p, eps) -> recovered state, or None for general

    Returns
    -------
    results : list of dicts with all computed quantities
    """
    results = []
    violations = 0

    for eps in eps_values:
        sigma = sigma_fn(eps)

        for state_name, ket in INPUT_STATES.items():
            rho = ket_to_dm(ket)

            # Check support condition before looping
            D_rho_sigma = umegaki_relative_entropy(rho, sigma)
            if np.isinf(D_rho_sigma):
                print(f"  [{channel_name}] eps={eps}, state={state_name}: "
                      f"D(rho||sigma)=inf, SKIPPING (support condition)")
                continue

            for p in p_values:
                kraus = kraus_fn(p)

                # Channel output
                N_rho = apply_channel(kraus, rho)
                N_sigma = apply_channel(kraus, sigma)

                # Relative entropies
                D_before = umegaki_relative_entropy(rho, sigma)
                D_after = umegaki_relative_entropy(N_rho, N_sigma)

                if np.isinf(D_before) or np.isinf(D_after):
                    continue

                Sigma = D_before - D_after  # entropy production / info loss

                # Sigma should be >= 0 by data processing inequality
                if Sigma < -1e-10:
                    print(f"  WARNING: DPI violation! Sigma={Sigma:.6e} "
                          f"at p={p}, eps={eps}, state={state_name}")

                Sigma = max(Sigma, 0)  # clip numerical noise

                # Petz recovery
                if petz_fn is not None:
                    recovered = petz_fn(N_rho, p, eps)
                else:
                    recovered = petz_recovery_general(rho, kraus, sigma)

                # Fidelity
                F = uhlmann_fidelity(rho, recovered)

                # Bound
                bound = np.exp(-Sigma / 2)

                # tau values
                tau = 1 - F
                tau_JRSWW = 1 - bound

                # Check bound
                holds = F >= bound - 1e-10  # small tolerance for numerics
                if not holds:
                    violations += 1
                    print(f"  *** VIOLATION *** at p={p:.2f}, eps={eps}, "
                          f"state={state_name}: F={F:.8f} < bound={bound:.8f}")

                row = {
                    'channel': channel_name,
                    'state': state_name,
                    'eps': eps,
                    'p': p,
                    'D_before': D_before,
                    'D_after': D_after,
                    'Sigma': Sigma,
                    'F_petz': F,
                    'bound': bound,
                    'tau_petz': tau,
                    'tau_JRSWW': tau_JRSWW,
                    'gap': F - bound,
                    'holds': holds,
                }
                results.append(row)

    return results, violations


def dephasing_special_analysis(eps_values, p_values):
    """
    Special analysis for dephasing channel comparing three tau values:
        tau_Petz:    exact Petz recovery for |+> input
        tau_optimal: identity map recovery (best for dephasing)
        tau_JRSWW:   the JRSWW bound 1 - exp(-Sigma/2)

    For dephasing with |+> input and sigma = (1-eps)|+><+| + eps|-><-|:
        - N(|+><+|) has fidelity with |+> equal to (1+1-p)/2 = 1-p/2
          Actually, let's compute exactly.
    """
    results = []

    print("\n" + "=" * 80)
    print("DEPHASING CHANNEL SPECIAL ANALYSIS: tau_Petz vs tau_optimal vs tau_JRSWW")
    print("=" * 80)

    for eps in eps_values:
        sigma = (1 - eps) * ket_to_dm(ket_plus) + eps * ket_to_dm(ket_minus)
        rho_plus = ket_to_dm(ket_plus)

        print(f"\n--- eps = {eps} ---")
        print(f"{'p':>6s}  {'tau_Petz':>12s}  {'tau_optimal':>12s}  "
              f"{'tau_JRSWW':>12s}  {'F_Petz':>12s}  {'Sigma':>12s}")
        print("-" * 75)

        for p in p_values:
            kraus = pd_kraus(p)

            # Channel output of |+>
            N_rho = apply_channel(kraus, rho_plus)
            N_sigma = apply_channel(kraus, sigma)

            # --- tau_optimal: identity recovery (just use N(rho) directly) ---
            # For dephasing, N(|+><+|) = [[0.5, 0.5(1-p)], [0.5(1-p), 0.5]]
            # F(|+>, N(|+>)) = <+|N(|+><+|)|+> = 0.5 + 0.5(1-p) = 1 - p/2
            F_identity = uhlmann_fidelity(rho_plus, N_rho)
            tau_optimal = 1 - F_identity

            # --- tau_Petz: Petz recovery ---
            recovered = petz_recovery_general(rho_plus, kraus, sigma)
            F_petz = uhlmann_fidelity(rho_plus, recovered)
            tau_petz = 1 - F_petz

            # --- Sigma and tau_JRSWW ---
            D_before = umegaki_relative_entropy(rho_plus, sigma)
            D_after = umegaki_relative_entropy(N_rho, N_sigma)

            if np.isinf(D_before) or np.isinf(D_after):
                continue

            Sigma = max(D_before - D_after, 0)
            tau_JRSWW = 1 - np.exp(-Sigma / 2)

            print(f"{p:6.2f}  {tau_petz:12.8f}  {tau_optimal:12.8f}  "
                  f"{tau_JRSWW:12.8f}  {F_petz:12.8f}  {Sigma:12.8f}")

            results.append({
                'eps': eps,
                'p': p,
                'tau_petz': tau_petz,
                'tau_optimal': tau_optimal,
                'tau_JRSWW': tau_JRSWW,
                'F_petz': F_petz,
                'Sigma': Sigma,
            })

    return results


def print_results_table(results, channel_name):
    """Print formatted results table for a channel."""
    print(f"\n{'=' * 110}")
    print(f"CHANNEL: {channel_name}")
    print(f"{'=' * 110}")

    current_eps = None
    current_state = None

    for r in results:
        if r['channel'] != channel_name:
            continue

        if r['eps'] != current_eps or r['state'] != current_state:
            current_eps = r['eps']
            current_state = r['state']
            print(f"\n  eps={r['eps']}, state={r['state']}")
            print(f"  {'p':>6s}  {'F_petz':>12s}  {'bound':>12s}  "
                  f"{'gap':>12s}  {'Sigma':>12s}  {'tau_petz':>12s}  "
                  f"{'tau_JRSWW':>12s}  {'holds':>6s}")
            print(f"  {'-'*90}")

        status = "YES" if r['holds'] else "*** NO ***"
        print(f"  {r['p']:6.2f}  {r['F_petz']:12.8f}  {r['bound']:12.8f}  "
              f"{r['gap']:12.8f}  {r['Sigma']:12.8f}  {r['tau_petz']:12.8f}  "
              f"{r['tau_JRSWW']:12.8f}  {status:>6s}")


def save_csv(results, filename):
    """Save results to CSV file."""
    if not results:
        return
    filepath = os.path.join(_DIR, filename)
    keys = results[0].keys()
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved: {filepath}")


def main():
    print("=" * 110)
    print("SUPPLEMENT VERIFICATION: F >= exp(-Sigma/2) for Singh et al. (2025) channels")
    print("=" * 110)

    eps_values = [0.2, 0.5, 0.8]
    p_values = np.arange(0.0, 1.05, 0.05)
    # Clean up floating point: round to avoid 0.050000000000000044 etc.
    p_values = np.round(p_values, 4)

    # --- Amplitude Damping ---
    print("\n\n>>> AMPLITUDE DAMPING CHANNEL <<<")
    print("    K0 = [[1,0],[0,sqrt(1-p)]], K1 = [[0,sqrt(p)],[0,0]]")
    print("    sigma = (1-eps)|0><0| + eps|1><1|")

    def ad_sigma(eps):
        return (1 - eps) * ket_to_dm(ket0) + eps * ket_to_dm(ket1)

    def ad_petz(N_rho, p, eps):
        """Petz recovery for AD using explicit Kraus operators."""
        N_sigma_00 = 1 - (1 - p) * eps
        if N_sigma_00 < 1e-15:
            return N_rho

        M0 = np.array([
            [np.sqrt((1 - eps) / N_sigma_00), 0],
            [0, 1]
        ], dtype=complex)

        M1 = np.array([
            [0, 0],
            [np.sqrt(p * eps / N_sigma_00), 0]
        ], dtype=complex)

        return apply_channel([M0, M1], N_rho)

    ad_results, ad_violations = verify_channel(
        "Amplitude Damping", ad_kraus, ad_sigma, eps_values, p_values,
        petz_fn=ad_petz
    )

    # Also verify with general Petz for cross-check (just for a subset)
    print("\n  Cross-checking AD explicit Kraus vs general Petz (subset)...")
    for eps in [0.5]:
        sigma = ad_sigma(eps)
        for state_name, ket in INPUT_STATES.items():
            rho = ket_to_dm(ket)
            D_check = umegaki_relative_entropy(rho, sigma)
            if np.isinf(D_check):
                continue
            for p in [0.1, 0.5, 0.9]:
                kraus = ad_kraus(p)
                N_rho = apply_channel(kraus, rho)
                rec_explicit = ad_petz(N_rho, p, eps)
                rec_general = petz_recovery_general(rho, kraus, sigma)
                F_explicit = uhlmann_fidelity(rho, rec_explicit)
                F_general = uhlmann_fidelity(rho, rec_general)
                diff = abs(F_explicit - F_general)
                if diff > 1e-6:
                    print(f"    MISMATCH at p={p}, eps={eps}, state={state_name}: "
                          f"F_explicit={F_explicit:.8f}, F_general={F_general:.8f}, "
                          f"diff={diff:.2e}")
                else:
                    print(f"    OK: p={p}, eps={eps}, state={state_name}: "
                          f"F_explicit={F_explicit:.8f}, F_general={F_general:.8f}")

    print_results_table(ad_results, "Amplitude Damping")

    # --- Phase Damping ---
    print("\n\n>>> PHASE DAMPING (DEPHASING) CHANNEL <<<")
    print("    K0 = sqrt(1-p/2)*I, K1 = sqrt(p/2)*sigma_z")
    print("    sigma = (1-eps)|+><+| + eps|-><-|")

    def pd_sigma(eps):
        return (1 - eps) * ket_to_dm(ket_plus) + eps * ket_to_dm(ket_minus)

    pd_results, pd_violations = verify_channel(
        "Phase Damping", pd_kraus, pd_sigma, eps_values, p_values
    )

    print_results_table(pd_results, "Phase Damping")

    # --- Dephasing special analysis ---
    deph_special = dephasing_special_analysis(eps_values, p_values)

    # --- Summary ---
    all_results = ad_results + pd_results
    total = len(all_results)
    total_violations = ad_violations + pd_violations

    print("\n" + "=" * 110)
    print("SUMMARY")
    print("=" * 110)
    print(f"  Total data points computed:      {total}")
    print(f"  Amplitude Damping points:        {len(ad_results)}")
    print(f"  Phase Damping points:            {len(pd_results)}")
    print(f"  Total violations of F>=exp(-S/2): {total_violations}")

    if total_violations == 0:
        print("\n  *** ALL POINTS SATISFY F >= exp(-Sigma/2) ***")
        print("  The JRSWW bound is verified for all tested configurations.")
    else:
        print(f"\n  *** WARNING: {total_violations} VIOLATIONS FOUND ***")

    # Statistics
    if all_results:
        gaps = [r['gap'] for r in all_results]
        taus = [r['tau_petz'] for r in all_results]
        sigmas = [r['Sigma'] for r in all_results]

        print(f"\n  Gap statistics (F - bound):")
        print(f"    Min gap:  {min(gaps):.8f}")
        print(f"    Max gap:  {max(gaps):.8f}")
        print(f"    Mean gap: {np.mean(gaps):.8f}")

        print(f"\n  tau_Petz statistics:")
        print(f"    Min:  {min(taus):.8f}")
        print(f"    Max:  {max(taus):.8f}")
        print(f"    Mean: {np.mean(taus):.8f}")

        print(f"\n  Sigma (entropy production) statistics:")
        print(f"    Min:  {min(sigmas):.8f}")
        print(f"    Max:  {max(sigmas):.8f}")
        print(f"    Mean: {np.mean(sigmas):.8f}")

    # Identify tightest bound cases (smallest gap)
    if all_results:
        sorted_by_gap = sorted(all_results, key=lambda r: r['gap'])
        print(f"\n  Top 5 tightest bound cases (smallest gap F - bound):")
        for i, r in enumerate(sorted_by_gap[:5]):
            print(f"    {i+1}. {r['channel']}, state={r['state']}, "
                  f"eps={r['eps']}, p={r['p']:.2f}: "
                  f"gap={r['gap']:.10f}, F={r['F_petz']:.8f}, "
                  f"bound={r['bound']:.8f}")

    # --- Save CSVs ---
    save_csv(ad_results, "supplement_verification_AD.csv")
    save_csv(pd_results, "supplement_verification_PD.csv")
    save_csv(deph_special, "supplement_verification_dephasing_special.csv")

    # Save combined
    save_csv(all_results, "supplement_verification_all.csv")

    print("\n" + "=" * 110)
    print("VERIFICATION COMPLETE")
    print("=" * 110)

    return total_violations


if __name__ == "__main__":
    n_violations = main()
    if n_violations > 0:
        exit(1)
