"""
IFF Saturation Search: Find the exact necessary and sufficient condition for
    F²(ρ, R_Petz ∘ N(ρ)) = exp(−ΔD)

where:
    - R_Petz is the Petz recovery map with reference σ = I/2
    - F is the fidelity
    - ΔD = D(ρ‖σ) - D(N(ρ)‖N(σ))
    - D is the quantum relative entropy

For σ = I/2 (maximally mixed state), many formulas simplify.

Author: Sheng-Kai Huang
Date: 2026-03-08
"""

import numpy as np
from scipy.linalg import sqrtm, logm, expm
from scipy.stats import unitary_group
import warnings
import time
from collections import defaultdict

warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def random_pure_state(d=2):
    """Random pure state |ψ⟩ as a d-dim vector."""
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return psi

def pure_state_dm(psi):
    """Pure state density matrix |ψ⟩⟨ψ|."""
    return np.outer(psi, psi.conj())

def random_cptp_stinespring(d=2, d_env=None):
    """
    Generate random CPTP map via Stinespring dilation.
    V is a (d*d_env, d) isometry, N(rho) = Tr_E(V rho V^dag).
    Kraus operators: K_k = <k_E| V.
    """
    if d_env is None:
        d_env = np.random.choice([2, 3, 4])

    # Random unitary on d*d_env, take first d columns
    U = np.random.randn(d * d_env, d * d_env) + 1j * np.random.randn(d * d_env, d * d_env)
    Q, R = np.linalg.qr(U)
    V = Q[:, :d]  # isometry: V^dag V = I_d

    kraus = []
    for k in range(d_env):
        K = V[k * d:(k + 1) * d, :]  # shape (d, d)
        kraus.append(K)

    return kraus

def apply_channel(kraus, rho):
    """Apply quantum channel N(ρ) = Σ_k K_k ρ K_k†."""
    d = rho.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K @ rho @ K.conj().T
    return result

def apply_adjoint(kraus, X):
    """Apply adjoint channel N†(X) = Σ_k K_k† X K_k."""
    d = X.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K.conj().T @ X @ K
    return result

def matrix_log(A):
    """Compute matrix logarithm, handling near-singular cases."""
    return logm(A)

def matrix_sqrt(A):
    """Compute matrix square root."""
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 0)
    return eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.conj().T

def matrix_inv_sqrt(A):
    """Compute A^{-1/2} for positive definite A."""
    eigvals, eigvecs = np.linalg.eigh(A)
    mask = eigvals > 1e-14
    inv_sqrt_vals = np.zeros_like(eigvals)
    inv_sqrt_vals[mask] = 1.0 / np.sqrt(eigvals[mask])
    return eigvecs @ np.diag(inv_sqrt_vals) @ eigvecs.conj().T

def fidelity_squared(rho, sigma):
    """
    Compute F²(ρ, σ) = (Tr √(√ρ σ √ρ))².
    For pure ρ = |ψ⟩⟨ψ|: F²(ρ, σ) = ⟨ψ|σ|ψ⟩.
    """
    sqrt_rho = matrix_sqrt(rho)
    inner = sqrt_rho @ sigma @ sqrt_rho
    # Eigenvalues of inner should be non-negative
    eigvals = np.linalg.eigvalsh(inner)
    eigvals = np.maximum(eigvals, 0)
    F = np.sum(np.sqrt(eigvals))
    return F**2

def fidelity_squared_pure(psi, sigma):
    """F²(|ψ⟩⟨ψ|, σ) = ⟨ψ|σ|ψ⟩ (real for Hermitian σ)."""
    return np.real(psi.conj() @ sigma @ psi)

def quantum_relative_entropy(rho, sigma):
    """
    D(ρ‖σ) = Tr(ρ(ln ρ - ln σ)).
    For pure ρ = |ψ⟩⟨ψ|: D(ρ‖σ) = -⟨ψ|ln σ|ψ⟩ (since S(ρ)=0, Tr(ρ ln ρ) = 0).

    Requires supp(ρ) ⊆ supp(σ), otherwise D = +∞.
    """
    # Check support condition
    eigvals_s = np.linalg.eigvalsh(sigma)
    eigvals_r = np.linalg.eigvalsh(rho)

    # Compute in eigenbasis of sigma
    log_sigma = logm(sigma)
    log_rho = logm(rho + 1e-300 * np.eye(rho.shape[0]))  # regularize

    D = np.real(np.trace(rho @ (log_rho - log_sigma)))
    return D

def relative_entropy_pure(psi, sigma):
    """D(|ψ⟩⟨ψ| ‖ σ) = -⟨ψ|ln σ|ψ⟩."""
    log_sigma = logm(sigma)
    return -np.real(psi.conj() @ log_sigma @ psi)

def petz_recovery_map(kraus, sigma, X):
    """
    Petz recovery map: R_{σ,N}(X) = σ^{1/2} N†(N(σ)^{-1/2} X N(σ)^{-1/2}) σ^{1/2}

    For σ = I/2:
      σ^{1/2} = I/√2
      N(σ) = N(I/2) = (1/2) Σ_k K_k K_k†
      R(X) = (1/2) N†(N(I/2)^{-1/2} X N(I/2)^{-1/2})
    """
    d = X.shape[0]

    # Compute N(σ)
    N_sigma = apply_channel(kraus, sigma)

    # Compute N(σ)^{-1/2}
    N_sigma_inv_sqrt = matrix_inv_sqrt(N_sigma)

    # Compute N(σ)^{-1/2} X N(σ)^{-1/2}
    sandwiched = N_sigma_inv_sqrt @ X @ N_sigma_inv_sqrt

    # Apply adjoint
    adj_result = apply_adjoint(kraus, sandwiched)

    # Sandwich with σ^{1/2}
    sigma_sqrt = matrix_sqrt(sigma)
    result = sigma_sqrt @ adj_result @ sigma_sqrt

    return result

def compute_saturation_data(psi, kraus, sigma):
    """
    Compute all relevant quantities for a (ρ=|ψ⟩⟨ψ|, N, σ) triple.
    Returns dict with F², exp(-ΔD), ratio, and diagnostic quantities.
    """
    d = len(psi)
    rho = pure_state_dm(psi)

    # Apply channel
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)

    # Check if N(σ) is invertible
    eigvals_Ns = np.linalg.eigvalsh(N_sigma)
    if np.min(eigvals_Ns) < 1e-12:
        return None  # Skip non-invertible cases

    # Compute Petz recovery
    recovered = petz_recovery_map(kraus, sigma, N_rho)

    # F²(ρ, R∘N(ρ))
    F2 = fidelity_squared_pure(psi, recovered)

    # ΔD = D(ρ‖σ) - D(N(ρ)‖N(σ))
    D_before = relative_entropy_pure(psi, sigma)
    D_after = quantum_relative_entropy(N_rho, N_sigma)
    Delta_D = D_before - D_after

    if Delta_D < -1e-10:
        return None  # DPI violation, numerical issue

    exp_neg_Delta_D = np.exp(-Delta_D)

    # Ratio
    if exp_neg_Delta_D > 1e-15:
        ratio = F2 / exp_neg_Delta_D
    else:
        ratio = float('inf')

    return {
        'F2': F2,
        'exp_neg_DeltaD': exp_neg_Delta_D,
        'ratio': ratio,
        'Delta_D': Delta_D,
        'N_rho': N_rho,
        'recovered': recovered,
        'rho': rho,
        'N_sigma': N_sigma,
    }

# ============================================================
# PHASE 1: SYSTEMATIC QUBIT SURVEY
# ============================================================

def phase1_survey(n_channels=10000, n_states=100, tol=1e-10):
    """Sweep over random CPTP qubit maps and pure input states."""
    print("=" * 70)
    print("PHASE 1: SYSTEMATIC QUBIT SURVEY")
    print(f"  {n_channels} channels × {n_states} states = {n_channels * n_states} pairs")
    print("=" * 70)

    d = 2
    sigma = np.eye(d, dtype=complex) / d  # I/2

    saturating_pairs = []
    near_saturating = []  # ratio < 1.001
    total_tested = 0
    total_valid = 0

    t0 = time.time()

    for i_ch in range(n_channels):
        if (i_ch + 1) % 1000 == 0:
            elapsed = time.time() - t0
            print(f"  Channel {i_ch+1}/{n_channels}  "
                  f"({elapsed:.1f}s elapsed, "
                  f"{len(saturating_pairs)} saturating found)")

        # Generate random CPTP map via Stinespring
        try:
            kraus = random_cptp_stinespring(d)
        except:
            continue

        # Pre-compute N(σ) and check invertibility
        N_sigma = apply_channel(kraus, sigma)
        if np.min(np.linalg.eigvalsh(N_sigma)) < 1e-12:
            continue

        # Compute Bloch parameters for this channel
        # N(ρ) = (I + (t + Λn)·σ)/2 where ρ = (I + n·σ)/2
        pauli = [np.array([[0,1],[1,0]], dtype=complex),
                 np.array([[0,-1j],[1j,0]], dtype=complex),
                 np.array([[1,0],[0,-1]], dtype=complex)]

        # Translation vector t: N(I/2) = (I + t·σ)/2
        t_vec = np.array([0.5 * np.real(np.trace(p @ N_sigma)) for p in pauli])

        # Contraction matrix Λ: apply N to σ_i/2 states
        Lambda = np.zeros((3, 3))
        for j in range(3):
            rho_j_plus = (np.eye(d) + pauli[j]) / 2
            N_rho_j = apply_channel(kraus, rho_j_plus)
            bloch_out = np.array([0.5 * np.real(np.trace(p @ N_rho_j)) for p in pauli])
            Lambda[:, j] = bloch_out - t_vec

        for i_st in range(n_states):
            total_tested += 1
            psi = random_pure_state(d)

            data = compute_saturation_data(psi, kraus, sigma)
            if data is None:
                continue
            total_valid += 1

            if abs(data['ratio'] - 1.0) < tol:
                # SATURATING!
                # Compute Bloch vector of input
                n_vec = np.array([np.real(np.trace(p @ data['rho'])) for p in pauli])

                saturating_pairs.append({
                    'channel_idx': i_ch,
                    'state_idx': i_st,
                    'psi': psi.copy(),
                    'kraus': kraus,
                    'F2': data['F2'],
                    'exp_neg_DeltaD': data['exp_neg_DeltaD'],
                    'ratio': data['ratio'],
                    'Delta_D': data['Delta_D'],
                    'N_rho': data['N_rho'],
                    'recovered': data['recovered'],
                    'rho': data['rho'],
                    'N_sigma': N_sigma,
                    't_vec': t_vec,
                    'Lambda': Lambda,
                    'n_vec': n_vec,
                    'kraus_rank': len(kraus),
                })
            elif data['ratio'] < 1.001:
                near_saturating.append({
                    'ratio': data['ratio'],
                    'Delta_D': data['Delta_D'],
                })

    elapsed = time.time() - t0
    print(f"\nPhase 1 complete in {elapsed:.1f}s")
    print(f"  Total tested: {total_tested}")
    print(f"  Total valid: {total_valid}")
    print(f"  Saturating (|ratio-1| < {tol}): {len(saturating_pairs)}")
    print(f"  Near-saturating (ratio < 1.001): {len(near_saturating)}")

    return saturating_pairs, near_saturating


# ============================================================
# PHASE 2: ANALYZE SATURATING PAIRS
# ============================================================

def phase2_analyze(saturating_pairs):
    """Deep analysis of each saturating pair."""
    print("\n" + "=" * 70)
    print("PHASE 2: ANALYZE SATURATING PAIRS")
    print("=" * 70)

    if len(saturating_pairs) == 0:
        print("  No saturating pairs found!")
        return {}

    d = 2
    sigma = np.eye(d, dtype=complex) / d
    pauli = [np.array([[0,1],[1,0]], dtype=complex),
             np.array([[0,-1j],[1j,0]], dtype=complex),
             np.array([[1,0],[0,-1]], dtype=complex)]

    properties = defaultdict(list)

    for idx, pair in enumerate(saturating_pairs):
        psi = pair['psi']
        kraus = pair['kraus']
        rho = pair['rho']
        N_rho = pair['N_rho']
        N_sigma = pair['N_sigma']
        recovered = pair['recovered']
        t_vec = pair['t_vec']
        Lambda = pair['Lambda']
        n_vec = pair['n_vec']

        # 1. Is ρ a fixed point? ||N(ρ) - ρ|| < ε
        is_fixed = np.linalg.norm(N_rho - rho) < 1e-8
        properties['is_fixed_point'].append(is_fixed)

        # 2. Is N unital? ||N(I/2) - I/2|| < ε
        is_unital = np.linalg.norm(t_vec) < 1e-8
        properties['is_unital'].append(is_unital)

        # 3. Does [ρ, N(ρ)] = 0?
        comm_rho_Nrho = rho @ N_rho - N_rho @ rho
        commutes = np.linalg.norm(comm_rho_Nrho) < 1e-8
        properties['rho_Nrho_commute'].append(commutes)

        # 4. Does [ρ, recovered] = 0?
        comm_rho_rec = rho @ recovered - recovered @ rho
        commutes_rec = np.linalg.norm(comm_rho_rec) < 1e-8
        properties['rho_recovered_commute'].append(commutes_rec)

        # 5. Eigenvalues of N(ρ)
        eigvals_Nrho = np.linalg.eigvalsh(N_rho)
        properties['eigvals_Nrho'].append(eigvals_Nrho)

        # 6. Kraus rank
        kraus_rank = len(kraus)
        properties['kraus_rank'].append(kraus_rank)

        # 7. (removed choi spectrum, using kraus_rank above instead)

        # 8. Is ρ an eigenstate of N†?
        # N†(ρ) = c·ρ for some scalar c?
        adj_rho = apply_adjoint(kraus, rho)
        # Check if adj_rho ∝ ρ
        # For pure ρ: adj_rho = c·ρ iff adj_rho = Tr(adj_rho · ρ) · ρ + ...
        c = np.real(np.trace(adj_rho @ rho))
        residual_adj = adj_rho - c * rho
        is_adj_eigenstate = np.linalg.norm(residual_adj) < 1e-8
        properties['is_adjoint_eigenstate'].append(is_adj_eigenstate)
        properties['adjoint_eigenvalue'].append(c)

        # 9. ΔD value
        properties['Delta_D'].append(pair['Delta_D'])

        # 10. Is ΔD = 0? (i.e., no information loss)
        is_zero_loss = pair['Delta_D'] < 1e-10
        properties['is_zero_DeltaD'].append(is_zero_loss)

        # 11. Bloch vector analysis
        Ln = Lambda @ n_vec  # Should be the Bloch vector shift
        n_out = t_vec + Ln  # Output Bloch vector
        properties['bloch_in'].append(n_vec)
        properties['bloch_out'].append(n_out)

        # 12. Are input and output Bloch vectors parallel?
        if np.linalg.norm(n_vec) > 1e-10 and np.linalg.norm(n_out) > 1e-10:
            cos_angle = np.dot(n_vec, n_out) / (np.linalg.norm(n_vec) * np.linalg.norm(n_out))
            is_parallel = abs(abs(cos_angle) - 1) < 1e-8
        else:
            cos_angle = 0
            is_parallel = False
        properties['bloch_parallel'].append(is_parallel)
        properties['bloch_cos_angle'].append(cos_angle)

        # 13. Is N(ρ) = α·ρ + (1-α)·σ for some α?
        # This means N(ρ) - σ = α(ρ - σ), i.e., output is on the line ρ-σ
        diff_out = N_rho - sigma
        diff_in = rho - sigma
        if np.linalg.norm(diff_in) > 1e-12:
            # Check if diff_out = α * diff_in
            # Use Frobenius inner product
            alpha = np.real(np.trace(diff_out.conj().T @ diff_in)) / np.real(np.trace(diff_in.conj().T @ diff_in))
            residual_line = diff_out - alpha * diff_in
            on_line = np.linalg.norm(residual_line) < 1e-8
        else:
            alpha = 0
            on_line = False
        properties['N_rho_on_rho_sigma_line'].append(on_line)
        properties['alpha_contraction'].append(alpha)

        # 14. Singular values of Λ
        sv_Lambda = np.linalg.svd(Lambda, compute_uv=False)
        properties['Lambda_singular_values'].append(sv_Lambda)

        # 15. Is Λ proportional to identity? (isotropic channel)
        Lambda_iso_residual = Lambda - (np.trace(Lambda)/3) * np.eye(3)
        is_isotropic = np.linalg.norm(Lambda_iso_residual) < 1e-8
        properties['is_isotropic'].append(is_isotropic)

        # 16. Is Λ diagonal in some basis related to n_vec?
        # Check if n_vec is an eigenvector of Λ^T Λ
        LtL = Lambda.T @ Lambda
        LtL_n = LtL @ n_vec
        if np.linalg.norm(n_vec) > 1e-10:
            proj = np.dot(LtL_n, n_vec) / np.dot(n_vec, n_vec)
            residual_ev = LtL_n - proj * n_vec
            n_is_eigvec_LtL = np.linalg.norm(residual_ev) < 1e-8
        else:
            n_is_eigvec_LtL = False
        properties['n_eigvec_of_LtL'].append(n_is_eigvec_LtL)

        # 17. Check: does Λ·n̂ ∥ n̂? (i.e., n̂ is eigenvector of Λ)
        Ln_normalized = Lambda @ n_vec
        if np.linalg.norm(Ln_normalized) > 1e-10 and np.linalg.norm(n_vec) > 1e-10:
            cos_Ln_n = np.dot(Ln_normalized, n_vec) / (np.linalg.norm(Ln_normalized) * np.linalg.norm(n_vec))
            n_eigvec_Lambda = abs(abs(cos_Ln_n) - 1) < 1e-8
        else:
            n_eigvec_Lambda = False
            cos_Ln_n = 0
        properties['n_eigvec_of_Lambda'].append(n_eigvec_Lambda)

        # 18. If n̂ is eigvec of Λ, what's the eigenvalue?
        if n_eigvec_Lambda and np.linalg.norm(n_vec) > 1e-10:
            lambda_n = np.dot(Lambda @ n_vec, n_vec) / np.dot(n_vec, n_vec)
        else:
            lambda_n = None
        properties['lambda_n_eigenvalue'].append(lambda_n)

        # 19. Check t ⊥ n̂ (translation perpendicular to input Bloch vector)
        if np.linalg.norm(t_vec) > 1e-10 and np.linalg.norm(n_vec) > 1e-10:
            t_dot_n = abs(np.dot(t_vec, n_vec)) / (np.linalg.norm(t_vec) * np.linalg.norm(n_vec))
            t_perp_n = t_dot_n < 1e-8
            t_parallel_n = abs(t_dot_n - 1) < 1e-8
        else:
            t_perp_n = True  # t=0 is vacuously perp
            t_parallel_n = np.linalg.norm(t_vec) < 1e-10
            t_dot_n = 0
        properties['t_perp_n'].append(t_perp_n)
        properties['t_parallel_n'].append(t_parallel_n)
        properties['t_dot_n_normalized'].append(t_dot_n)

        # 20. Check if the channel is "classical" on ρ:
        # Does there exist a basis containing |ψ⟩ such that N is classical?
        # Equivalent: does [K_i, |ψ⟩⟨ψ|] = 0 for all i?
        all_kraus_commute = True
        for K in kraus:
            comm = K @ rho - rho @ K
            if np.linalg.norm(comm) > 1e-8:
                all_kraus_commute = False
                break
        properties['all_kraus_commute_with_rho'].append(all_kraus_commute)

        # 21. Check: N(ρ) proportional to identity? (completely mixed output)
        is_Nrho_mixed = np.linalg.norm(N_rho - sigma) < 1e-8
        properties['N_rho_is_mixed'].append(is_Nrho_mixed)

        # 22. F² value itself
        properties['F2_value'].append(pair['F2'])

        # 23. Check reversibility: R∘N(ρ) = ρ exactly?
        is_perfectly_recovered = np.linalg.norm(recovered - rho) < 1e-8
        properties['perfectly_recovered'].append(is_perfectly_recovered)

        # 24. Check: is Λ symmetric?
        Lambda_sym_residual = Lambda - Lambda.T
        is_Lambda_sym = np.linalg.norm(Lambda_sym_residual) < 1e-8
        properties['Lambda_symmetric'].append(is_Lambda_sym)

        # 25. Compute: Λ^T Λ eigenvalues
        LtL_eigvals = np.linalg.eigvalsh(Lambda.T @ Lambda)
        properties['LtL_eigenvalues'].append(LtL_eigvals)

        # 26. NEW: Check if n is eigenvector of Λ^T (not just Λ)
        Lt_n = Lambda.T @ n_vec
        if np.linalg.norm(Lt_n) > 1e-10 and np.linalg.norm(n_vec) > 1e-10:
            cos_Ltn_n = np.dot(Lt_n, n_vec) / (np.linalg.norm(Lt_n) * np.linalg.norm(n_vec))
            n_eigvec_LambdaT = abs(abs(cos_Ltn_n) - 1) < 1e-8
        else:
            n_eigvec_LambdaT = False
        properties['n_eigvec_of_LambdaT'].append(n_eigvec_LambdaT)

        # 27. Check combined condition: t ∥ Λ·n (translation along output direction)
        if np.linalg.norm(t_vec) > 1e-10 and np.linalg.norm(Ln_normalized) > 1e-10:
            cos_t_Ln = abs(np.dot(t_vec, Ln_normalized)) / (np.linalg.norm(t_vec) * np.linalg.norm(Ln_normalized))
            t_parallel_Ln = abs(cos_t_Ln - 1) < 1e-8
        else:
            t_parallel_Ln = np.linalg.norm(t_vec) < 1e-10
        properties['t_parallel_Ln'].append(t_parallel_Ln)

        # 28. Combined: n is eigvec of Λ AND (t=0 OR t ∥ n)
        combined_cond1 = n_eigvec_Lambda and (np.linalg.norm(t_vec) < 1e-8 or t_parallel_n)
        properties['cond_n_eigvec_and_t_parallel'].append(combined_cond1)

    # Print summary
    print(f"\nAnalyzed {len(saturating_pairs)} saturating pairs:")
    print("-" * 50)

    bool_props = [
        'is_fixed_point', 'is_unital', 'rho_Nrho_commute', 'rho_recovered_commute',
        'is_adjoint_eigenstate', 'is_zero_DeltaD', 'bloch_parallel',
        'N_rho_on_rho_sigma_line', 'is_isotropic', 'n_eigvec_of_LtL',
        'n_eigvec_of_Lambda', 'n_eigvec_of_LambdaT',
        't_perp_n', 't_parallel_n',
        'all_kraus_commute_with_rho', 'N_rho_is_mixed',
        'perfectly_recovered', 'Lambda_symmetric',
        't_parallel_Ln', 'cond_n_eigvec_and_t_parallel',
    ]

    for prop in bool_props:
        if prop in properties:
            vals = properties[prop]
            n_true = sum(vals)
            pct = 100 * n_true / len(vals) if vals else 0
            marker = " <<<" if pct == 100 else (" !!!" if pct == 0 else "")
            print(f"  {prop:45s}: {n_true:5d}/{len(vals):5d}  ({pct:6.2f}%){marker}")

    print("\n--- Continuous properties ---")
    for prop in ['Delta_D', 'F2_value', 'alpha_contraction']:
        if prop in properties:
            vals = [v for v in properties[prop] if v is not None]
            if vals:
                print(f"  {prop}: min={min(vals):.6f}, max={max(vals):.6f}, "
                      f"mean={np.mean(vals):.6f}, std={np.std(vals):.6f}")

    return properties


# ============================================================
# PHASE 3: ALGEBRAIC BLOCH SPHERE ANALYSIS
# ============================================================

def phase3_algebraic(saturating_pairs, properties):
    """
    For σ = I/2, everything has clean Bloch sphere expressions.

    ρ = (I + n̂·σ⃗)/2, |n̂| = 1 (pure state)
    N(ρ) = (I + m⃗·σ⃗)/2 where m⃗ = t⃗ + Λ·n̂
    N(I/2) = (I + t⃗·σ⃗)/2

    For σ = I/2:
      D(ρ‖σ) = ln 2 - S(ρ) = ln 2  (for pure ρ)
      D(N(ρ)‖N(σ)) = relative entropy of two qubit states in Bloch form

    F²(ρ, R∘N(ρ)): need Petz map in Bloch form
    """
    print("\n" + "=" * 70)
    print("PHASE 3: ALGEBRAIC BLOCH SPHERE ANALYSIS")
    print("=" * 70)

    if len(saturating_pairs) == 0:
        print("No saturating pairs to analyze algebraically.")
        return

    # Categorize saturating pairs
    categories = defaultdict(list)

    for idx, pair in enumerate(saturating_pairs):
        t_norm = np.linalg.norm(pair['t_vec'])
        n_vec = pair['n_vec']
        Lambda = pair['Lambda']

        # Classify channel type
        if t_norm < 1e-8:
            # Unital channel
            if properties['is_isotropic'][idx]:
                cat = 'unital_isotropic'
            elif properties['n_eigvec_of_Lambda'][idx]:
                cat = 'unital_n_is_eigvec'
            else:
                cat = 'unital_other'
        else:
            if properties['N_rho_on_rho_sigma_line'][idx]:
                cat = 'nonunital_on_line'
            elif properties['n_eigvec_of_Lambda'][idx]:
                cat = 'nonunital_n_eigvec'
            else:
                cat = 'nonunital_other'

        categories[cat].append(idx)

    print("\nCategories of saturating pairs:")
    for cat, indices in sorted(categories.items()):
        print(f"  {cat}: {len(indices)} pairs")
        if len(indices) <= 5:
            for idx in indices:
                p = saturating_pairs[idx]
                print(f"    F²={p['F2']:.8f}, ΔD={p['Delta_D']:.8f}, "
                      f"|t|={np.linalg.norm(p['t_vec']):.6f}")

    # Deep dive: for each category, try to find the algebraic pattern
    print("\n--- Deep algebraic analysis ---")

    # For unital isotropic channels: Λ = λI, t=0
    # N(ρ) = (I + λn·σ)/2, N(I/2) = I/2
    # This is the depolarizing channel with parameter λ
    if 'unital_isotropic' in categories:
        print("\n[Unital Isotropic (Depolarizing-like)]")
        for idx in categories['unital_isotropic'][:3]:
            p = saturating_pairs[idx]
            sv = properties['Lambda_singular_values'][idx]
            print(f"  Λ singular values: {sv}")
            print(f"  F² = {p['F2']:.10f}")
            print(f"  ΔD = {p['Delta_D']:.10f}")
            # For depolarizing: N(ρ) = λρ + (1-λ)I/2
            # D(ρ‖I/2) = ln 2, D(N(ρ)‖N(I/2)) = D(N(ρ)‖I/2)
            # For pure input: N(ρ) has eigenvalues (1+λ)/2, (1-λ)/2
            # D(N(ρ)‖I/2) = ln 2 - S(N(ρ)) = ln 2 - H((1+λ)/2)
            # ΔD = H((1+λ)/2) = -(1+λ)/2 ln(1+λ)/2 - (1-λ)/2 ln(1-λ)/2

    # For unital channels where n is eigenvector of Λ
    if 'unital_n_is_eigvec' in categories:
        print("\n[Unital, n̂ eigenvector of Λ]")
        for idx in categories['unital_n_is_eigvec'][:5]:
            p = saturating_pairs[idx]
            n = p['n_vec']
            Ln = p['Lambda'] @ n
            lambda_n = np.dot(Ln, n) / np.dot(n, n)
            sv = properties['Lambda_singular_values'][idx]
            print(f"  Λ eigenvalue along n̂: {lambda_n:.8f}")
            print(f"  Λ singular values: {sv}")
            print(f"  F² = {p['F2']:.10f}, ΔD = {p['Delta_D']:.10f}")

    # Critical test: For ALL saturating pairs, check the combined condition
    print("\n--- CRITICAL: Testing candidate IFF conditions ---")

    # Condition A: n̂ is eigenvector of Λ AND t ∥ n̂ (or t=0)
    # This means the output Bloch vector m = t + Λn is parallel to n
    cond_A_count = sum(1 for i in range(len(saturating_pairs))
                       if properties['cond_n_eigvec_and_t_parallel'][i])

    # Condition B: output Bloch vector parallel to input
    cond_B_count = sum(1 for i in range(len(saturating_pairs))
                       if properties['bloch_parallel'][i])

    # Condition C: N(ρ) on the ρ-σ line (N(ρ) = αρ + (1-α)σ)
    cond_C_count = sum(1 for i in range(len(saturating_pairs))
                       if properties['N_rho_on_rho_sigma_line'][i])

    # Condition D: ΔD = 0 (trivial saturation, perfect recovery)
    cond_D_count = sum(1 for i in range(len(saturating_pairs))
                       if properties['is_zero_DeltaD'][i])

    # Condition E: [ρ, N(ρ)] = 0
    cond_E_count = sum(1 for i in range(len(saturating_pairs))
                       if properties['rho_Nrho_commute'][i])

    n_total = len(saturating_pairs)
    print(f"  Cond A (n̂ eigvec Λ, t∥n̂):    {cond_A_count}/{n_total}")
    print(f"  Cond B (m⃗ ∥ n̂):               {cond_B_count}/{n_total}")
    print(f"  Cond C (N(ρ)=αρ+(1-α)σ):       {cond_C_count}/{n_total}")
    print(f"  Cond D (ΔD=0):                  {cond_D_count}/{n_total}")
    print(f"  Cond E ([ρ,N(ρ)]=0):            {cond_E_count}/{n_total}")

    # Check: C ↔ B equivalence
    C_implies_B = all(properties['bloch_parallel'][i]
                      for i in range(n_total)
                      if properties['N_rho_on_rho_sigma_line'][i])
    B_implies_C = all(properties['N_rho_on_rho_sigma_line'][i]
                      for i in range(n_total)
                      if properties['bloch_parallel'][i])
    print(f"\n  C → B (on-line → parallel): {C_implies_B}")
    print(f"  B → C (parallel → on-line): {B_implies_C}")

    # For any non-100% condition, find counterexamples
    if cond_B_count < n_total:
        print(f"\n  Counterexamples to Cond B (m⃗ not ∥ n̂):")
        for i in range(min(5, n_total)):
            if not properties['bloch_parallel'][i]:
                p = saturating_pairs[i]
                print(f"    n̂={p['n_vec']}, m⃗={properties['bloch_out'][i]}, "
                      f"cos={properties['bloch_cos_angle'][i]:.8f}")

    if cond_C_count < n_total:
        print(f"\n  Counterexamples to Cond C (not on ρ-σ line):")
        count = 0
        for i in range(n_total):
            if not properties['N_rho_on_rho_sigma_line'][i] and count < 5:
                p = saturating_pairs[i]
                print(f"    n̂={p['n_vec']}, |t|={np.linalg.norm(p['t_vec']):.6f}, "
                      f"α={properties['alpha_contraction'][i]:.6f}")
                count += 1

    return categories


# ============================================================
# PHASE 4: VALIDATION
# ============================================================

def make_depolarizing_kraus(p, d=2):
    """Depolarizing channel: N(ρ) = (1-p)ρ + p·I/d. Parameter p ∈ [0, 1]."""
    # Kraus: sqrt(1-p)·I, sqrt(p/3)·σ_x, sqrt(p/3)·σ_y, sqrt(p/3)·σ_z
    # Actually for qubit: N(ρ) = (1-p)ρ + p·Tr(ρ)·I/2
    # = (1 - 4p/3)ρ + (4p/3)·I/2  if we use Pauli Kraus
    # Let's use: N(ρ) = (1-λ)·I/2 + λ·ρ where λ = 1-p
    pauli = [np.eye(2, dtype=complex),
             np.array([[0,1],[1,0]], dtype=complex),
             np.array([[0,-1j],[1j,0]], dtype=complex),
             np.array([[1,0],[0,-1]], dtype=complex)]

    lam = 1 - p  # contraction parameter
    # Kraus operators: sqrt((1+3λ/1)/4)·I and sqrt((1-λ)/4)·σ_i
    # Actually: depolarizing N(ρ) = (1-q)ρ + q·I/2
    # Kraus: K_0 = sqrt(1 - 3q/4)·I, K_i = sqrt(q/4)·σ_i
    q = p
    if q < 0 or q > 1:
        raise ValueError(f"p must be in [0,1], got {p}")

    K0 = np.sqrt(1 - 3*q/4) * pauli[0]
    K1 = np.sqrt(q/4) * pauli[1]
    K2 = np.sqrt(q/4) * pauli[2]
    K3 = np.sqrt(q/4) * pauli[3]
    return [K0, K1, K2, K3]

def make_amplitude_damping_kraus(gamma):
    """Amplitude damping channel."""
    K0 = np.array([[1, 0], [0, np.sqrt(1-gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]

def make_dephasing_kraus(p):
    """Dephasing channel: N(ρ) = (1-p)ρ + p·Z·ρ·Z."""
    K0 = np.sqrt(1-p) * np.eye(2, dtype=complex)
    K1 = np.sqrt(p) * np.array([[1,0],[0,-1]], dtype=complex)
    return [K0, K1]

def make_generalized_amplitude_damping_kraus(p, gamma):
    """Generalized amplitude damping."""
    K0 = np.sqrt(p) * np.array([[1, 0], [0, np.sqrt(1-gamma)]], dtype=complex)
    K1 = np.sqrt(p) * np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    K2 = np.sqrt(1-p) * np.array([[np.sqrt(1-gamma), 0], [0, 1]], dtype=complex)
    K3 = np.sqrt(1-p) * np.array([[0, 0], [np.sqrt(gamma), 0]], dtype=complex)
    return [K0, K1, K2, K3]

def make_erasure_kraus(p):
    """Erasure channel to qutrit: with prob p, erase to |2⟩.
    Actually for qubit-to-qubit analysis, erasure maps to a flag state.
    We model as: N(ρ) = (1-p)ρ + p·|0⟩⟨0| (replacement to |0⟩ with prob p).
    This is actually a "replacement channel" version.
    """
    # Actually let's do proper erasure-like: N(ρ) = (1-p)ρ + p·σ
    # where σ = I/2. This is just depolarizing!
    # For a proper erasure channel we'd need to go to higher dim.
    # Let's use replacement channel: N(ρ) = (1-p)ρ + p·|0⟩⟨0|
    # Kraus: K0 = sqrt(1-p)·I, then we need to find K1 such that K1 ρ K1† = p·|0⟩⟨0|
    # K1 = sqrt(p)·|0⟩⟨0|, K2 = sqrt(p)·|0⟩⟨1|
    K0 = np.sqrt(1-p) * np.eye(2, dtype=complex)
    K1 = np.sqrt(p) * np.array([[1,0],[0,0]], dtype=complex)
    K2 = np.sqrt(p) * np.array([[0,1],[0,0]], dtype=complex)
    return [K0, K1, K2]

def make_completely_depolarizing_kraus():
    """Completely depolarizing: N(ρ) = I/2 for all ρ."""
    pauli = [np.eye(2, dtype=complex),
             np.array([[0,1],[1,0]], dtype=complex),
             np.array([[0,-1j],[1j,0]], dtype=complex),
             np.array([[1,0],[0,-1]], dtype=complex)]
    return [0.5 * p for p in pauli]

def test_channel_saturation(kraus, channel_name, n_states=1000):
    """Test a specific channel for saturation across many states."""
    d = 2
    sigma = np.eye(d, dtype=complex) / d

    N_sigma = apply_channel(kraus, sigma)
    if np.min(np.linalg.eigvalsh(N_sigma)) < 1e-12:
        print(f"  {channel_name}: N(σ) not invertible, skipping")
        return

    pauli = [np.array([[0,1],[1,0]], dtype=complex),
             np.array([[0,-1j],[1j,0]], dtype=complex),
             np.array([[1,0],[0,-1]], dtype=complex)]

    t_vec = np.array([0.5 * np.real(np.trace(p @ N_sigma)) for p in pauli])
    Lambda = np.zeros((3, 3))
    for j in range(3):
        rho_j = (np.eye(d) + pauli[j]) / 2
        N_rho_j = apply_channel(kraus, rho_j)
        bloch_out = np.array([0.5 * np.real(np.trace(p @ N_rho_j)) for p in pauli])
        Lambda[:, j] = bloch_out - t_vec

    n_sat = 0
    min_ratio = float('inf')
    max_ratio = 0
    ratios = []
    sat_details = []

    for _ in range(n_states):
        psi = random_pure_state(d)
        data = compute_saturation_data(psi, kraus, sigma)
        if data is None:
            continue

        ratios.append(data['ratio'])
        if data['ratio'] < min_ratio:
            min_ratio = data['ratio']
        if data['ratio'] > max_ratio:
            max_ratio = data['ratio']

        if abs(data['ratio'] - 1.0) < 1e-8:
            n_sat += 1
            n_vec = np.array([np.real(np.trace(p @ data['rho'])) for p in pauli])

            # Check if n is eigvec of Lambda AND (t=0 or t∥n)
            Ln = Lambda @ n_vec
            if np.linalg.norm(Ln) > 1e-10 and np.linalg.norm(n_vec) > 1e-10:
                cos_Ln_n = np.dot(Ln, n_vec) / (np.linalg.norm(Ln) * np.linalg.norm(n_vec))
                n_eigvec = abs(abs(cos_Ln_n) - 1) < 1e-6
            else:
                n_eigvec = True  # Ln=0 means lambda=0, any n works

            on_line = False
            diff_out = data['N_rho'] - sigma
            diff_in = data['rho'] - sigma
            if np.linalg.norm(diff_in) > 1e-12:
                alpha = np.real(np.trace(diff_out.conj().T @ diff_in)) / np.real(np.trace(diff_in.conj().T @ diff_in))
                residual = diff_out - alpha * diff_in
                on_line = np.linalg.norm(residual) < 1e-8

            sat_details.append({
                'n_eigvec_Lambda': n_eigvec,
                'on_rho_sigma_line': on_line,
            })

    if ratios:
        all_sat = (n_sat == len(ratios))
        cond_check = ""
        if sat_details:
            all_eigvec = all(s['n_eigvec_Lambda'] for s in sat_details)
            all_on_line = all(s['on_rho_sigma_line'] for s in sat_details)
            cond_check = f" | n̂ eigvec(Λ): {all_eigvec} | on ρ-σ line: {all_on_line}"

        print(f"  {channel_name:35s}: {n_sat:4d}/{len(ratios):4d} saturating | "
              f"ratio range [{min_ratio:.10f}, {max_ratio:.10f}]{cond_check}")
        if not all_sat and n_sat > 0:
            print(f"    → PARTIAL SATURATION (interesting!)")

    return n_sat, len(ratios), sat_details


def phase4_validate(saturating_pairs, properties):
    """Validate conjectured condition on new data and special cases."""
    print("\n" + "=" * 70)
    print("PHASE 4: VALIDATION")
    print("=" * 70)

    d = 2
    sigma = np.eye(d, dtype=complex) / d

    # ---- Part A: Known channel families ----
    print("\n--- Part A: Known channel families ---")

    # Depolarizing
    for p in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        kraus = make_depolarizing_kraus(p)
        test_channel_saturation(kraus, f"Depolarizing(p={p})")

    print()

    # Amplitude damping
    for gamma in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        kraus = make_amplitude_damping_kraus(gamma)
        test_channel_saturation(kraus, f"AmplDamp(γ={gamma})")

    print()

    # Dephasing
    for p in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        kraus = make_dephasing_kraus(p)
        test_channel_saturation(kraus, f"Dephasing(p={p})")

    print()

    # Generalized amplitude damping
    for p, gamma in [(0.5, 0.3), (0.3, 0.5), (0.5, 0.5), (0.7, 0.3), (0.5, 1.0)]:
        kraus = make_generalized_amplitude_damping_kraus(p, gamma)
        test_channel_saturation(kraus, f"GAD(p={p},γ={gamma})")

    print()

    # Completely depolarizing
    kraus = make_completely_depolarizing_kraus()
    test_channel_saturation(kraus, "CompletelyDepolarizing")

    print()

    # Replacement/erasure-like
    for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
        kraus = make_erasure_kraus(p)
        test_channel_saturation(kraus, f"Replacement(p={p})")

    # ---- Part B: Targeted tests of conjectured condition ----
    print("\n--- Part B: Testing conjectured IFF condition ---")
    print("Conjecture: F² = exp(-ΔD) iff N(ρ) lies on the ρ-σ line")
    print("(equivalently: output Bloch vector m⃗ ∥ input n̂)")

    # Test: construct channels that deliberately satisfy/violate the condition

    # Test 1: Unital channel with Λ = diag(a, b, c)
    # n̂ along z-axis should always saturate (eigenvector of Λ)
    # n̂ along x-axis should saturate only if it's an eigenvector
    print("\n  Test: Diagonal Λ channels with specific input states")
    for a, b, c in [(0.5, 0.3, 0.7), (0.5, 0.5, 0.7), (0.3, 0.3, 0.3)]:
        # Construct unital channel with Λ = diag(a, b, c)
        # Kraus operators for this
        # For Pauli channel: N(ρ) = p0 ρ + p1 XρX + p2 YρY + p3 ZρZ
        # Λ = diag(p0+p1-p2-p3, p0-p1+p2-p3, p0-p1-p2+p3)
        # with p0+p1+p2+p3 = 1
        # a = 1-2(p2+p3), b = 1-2(p1+p3), c = 1-2(p1+p2)
        # → p1 = (1-b-c+a)/4... let me just verify
        # Actually: a = p0+p1-p2-p3, b = p0-p1+p2-p3, c = p0-p1-p2+p3
        # p0 = (1+a+b+c)/4, p1 = (1+a-b-c)/4, p2 = (1-a+b-c)/4, p3 = (1-a-b+c)/4
        p0 = (1+a+b+c)/4
        p1 = (1+a-b-c)/4
        p2 = (1-a+b-c)/4
        p3 = (1-a-b+c)/4

        if min(p0,p1,p2,p3) < -1e-10:
            continue

        pauli_ops = [np.eye(2, dtype=complex),
                     np.array([[0,1],[1,0]], dtype=complex),
                     np.array([[0,-1j],[1j,0]], dtype=complex),
                     np.array([[1,0],[0,-1]], dtype=complex)]

        kraus = [np.sqrt(max(pi, 0)) * P for pi, P in zip([p0,p1,p2,p3], pauli_ops)]
        kraus = [K for K in kraus if np.linalg.norm(K) > 1e-12]

        # Test eigenvector states
        for name, psi in [("z+", np.array([1,0], dtype=complex)),
                          ("z-", np.array([0,1], dtype=complex)),
                          ("x+", np.array([1,1], dtype=complex)/np.sqrt(2)),
                          ("y+", np.array([1,1j], dtype=complex)/np.sqrt(2)),
                          ("random", random_pure_state())]:
            data = compute_saturation_data(psi, kraus, sigma)
            if data is not None:
                sat = "SAT" if abs(data['ratio'] - 1) < 1e-8 else "---"
                rho = data['rho']
                n = np.array([np.real(np.trace(P @ rho)) for P in pauli_ops[1:]])
                Ln = np.diag([a,b,c]) @ n
                if np.linalg.norm(Ln) > 1e-10 and np.linalg.norm(n) > 1e-10:
                    cos = np.dot(Ln, n) / (np.linalg.norm(Ln) * np.linalg.norm(n))
                    is_eig = abs(abs(cos) - 1) < 1e-6
                else:
                    is_eig = True
                print(f"    Λ=diag({a},{b},{c}), |{name}⟩: {sat} ratio={data['ratio']:.10f} "
                      f"n̂_eigvec={is_eig}")

    # ---- Part C: 10000 new random channels ----
    print("\n--- Part C: 10000 new random channels validation ---")

    pauli_mats = [np.array([[0,1],[1,0]], dtype=complex),
                  np.array([[0,-1j],[1j,0]], dtype=complex),
                  np.array([[1,0],[0,-1]], dtype=complex)]

    n_new_channels = 10000
    n_states_per = 50
    n_conjecture_correct = 0
    n_conjecture_total = 0
    n_false_positive = 0  # Conjecture says sat but not
    n_false_negative = 0  # Conjecture says not sat but it is

    false_neg_examples = []
    false_pos_examples = []

    np.random.seed(12345)  # Different seed

    for i_ch in range(n_new_channels):
        if (i_ch + 1) % 2000 == 0:
            print(f"  Validated {i_ch+1}/{n_new_channels}...")

        try:
            kraus = random_cptp_stinespring(d)
        except:
            continue

        N_sigma = apply_channel(kraus, sigma)
        if np.min(np.linalg.eigvalsh(N_sigma)) < 1e-12:
            continue

        t_vec = np.array([0.5 * np.real(np.trace(p @ N_sigma)) for p in pauli_mats])
        Lambda = np.zeros((3, 3))
        for j in range(3):
            rho_j = (np.eye(d) + pauli_mats[j]) / 2
            N_rho_j = apply_channel(kraus, rho_j)
            bloch_out = np.array([0.5 * np.real(np.trace(p @ N_rho_j)) for p in pauli_mats])
            Lambda[:, j] = bloch_out - t_vec

        for i_st in range(n_states_per):
            psi = random_pure_state(d)
            data = compute_saturation_data(psi, kraus, sigma)
            if data is None:
                continue

            rho = data['rho']
            n_vec = np.array([np.real(np.trace(p @ rho)) for p in pauli_mats])

            # Check conjecture: N(ρ) on ρ-σ line
            diff_out = data['N_rho'] - sigma
            diff_in = rho - sigma
            if np.linalg.norm(diff_in) > 1e-12:
                alpha = np.real(np.trace(diff_out.conj().T @ diff_in)) / np.real(np.trace(diff_in.conj().T @ diff_in))
                residual = diff_out - alpha * diff_in
                conjecture_predicts_sat = np.linalg.norm(residual) < 1e-7
            else:
                conjecture_predicts_sat = True

            actual_sat = abs(data['ratio'] - 1.0) < 1e-7

            n_conjecture_total += 1
            if conjecture_predicts_sat == actual_sat:
                n_conjecture_correct += 1
            elif conjecture_predicts_sat and not actual_sat:
                n_false_positive += 1
                if len(false_pos_examples) < 5:
                    false_pos_examples.append({
                        'ratio': data['ratio'], 'residual': np.linalg.norm(residual),
                        'n_vec': n_vec, 't_vec': t_vec
                    })
            else:
                n_false_negative += 1
                if len(false_neg_examples) < 5:
                    false_neg_examples.append({
                        'ratio': data['ratio'], 'residual': np.linalg.norm(residual),
                        'n_vec': n_vec, 't_vec': t_vec
                    })

    accuracy = n_conjecture_correct / n_conjecture_total if n_conjecture_total > 0 else 0
    print(f"\n  Conjecture validation results:")
    print(f"    Total pairs tested:  {n_conjecture_total}")
    print(f"    Correct predictions: {n_conjecture_correct} ({100*accuracy:.4f}%)")
    print(f"    False positives:     {n_false_positive}")
    print(f"    False negatives:     {n_false_negative}")

    if false_pos_examples:
        print(f"\n  False positive examples (conjecture says sat, but not):")
        for ex in false_pos_examples:
            print(f"    ratio={ex['ratio']:.10f}, residual={ex['residual']:.2e}")

    if false_neg_examples:
        print(f"\n  False negative examples (conjecture says not sat, but is):")
        for ex in false_neg_examples:
            print(f"    ratio={ex['ratio']:.10f}, residual={ex['residual']:.2e}")

    return accuracy


# ============================================================
# PHASE 2.5: DEEPER ALGEBRAIC INVESTIGATION
# ============================================================

def phase2_5_deeper_analysis(saturating_pairs, properties):
    """
    If the simple conjecture fails, dig deeper.
    Compute more exotic quantities.
    """
    print("\n" + "=" * 70)
    print("PHASE 2.5: DEEPER ALGEBRAIC INVESTIGATION")
    print("=" * 70)

    if len(saturating_pairs) == 0:
        print("No saturating pairs.")
        return

    d = 2
    sigma = np.eye(d, dtype=complex) / d
    pauli = [np.array([[0,1],[1,0]], dtype=complex),
             np.array([[0,-1j],[1j,0]], dtype=complex),
             np.array([[1,0],[0,-1]], dtype=complex)]

    print("\nDetailed examination of first 20 saturating pairs:")

    for idx in range(min(20, len(saturating_pairs))):
        p = saturating_pairs[idx]
        psi = p['psi']
        kraus = p['kraus']
        rho = p['rho']
        N_rho = p['N_rho']
        N_sigma = p['N_sigma']
        n_vec = p['n_vec']
        t_vec = p['t_vec']
        Lambda = p['Lambda']

        m_vec = t_vec + Lambda @ n_vec  # output Bloch vector

        # 1. Decompose m along n and perp to n
        m_parallel = np.dot(m_vec, n_vec) * n_vec / np.dot(n_vec, n_vec)
        m_perp = m_vec - m_parallel

        # 2. Eigendecomposition of Λ
        eigvals_L, eigvecs_L = np.linalg.eig(Lambda)

        # 3. Check: is m_perp = 0?
        m_perp_norm = np.linalg.norm(m_perp)

        # 4. Check: is t_perp_to_n = 0? (component of t perpendicular to n)
        t_parallel = np.dot(t_vec, n_vec) * n_vec / np.dot(n_vec, n_vec)
        t_perp = t_vec - t_parallel
        t_perp_norm = np.linalg.norm(t_perp)

        # 5. Λ·n decomposition
        Ln = Lambda @ n_vec
        Ln_parallel = np.dot(Ln, n_vec) * n_vec / np.dot(n_vec, n_vec)
        Ln_perp = Ln - Ln_parallel
        Ln_perp_norm = np.linalg.norm(Ln_perp)

        if idx < 10:
            print(f"\n  Pair {idx}: F²={p['F2']:.10f}, ΔD={p['Delta_D']:.10f}")
            print(f"    n̂ = [{n_vec[0]:.6f}, {n_vec[1]:.6f}, {n_vec[2]:.6f}]")
            print(f"    t⃗ = [{t_vec[0]:.6f}, {t_vec[1]:.6f}, {t_vec[2]:.6f}]  |t|={np.linalg.norm(t_vec):.6f}")
            print(f"    m⃗ = [{m_vec[0]:.6f}, {m_vec[1]:.6f}, {m_vec[2]:.6f}]  |m|={np.linalg.norm(m_vec):.6f}")
            print(f"    |m⊥| = {m_perp_norm:.2e}, |t⊥| = {t_perp_norm:.2e}, |Λn⊥| = {Ln_perp_norm:.2e}")
            print(f"    Λ eigenvalues: {eigvals_L}")
            print(f"    Λ·n̂ = [{Ln[0]:.6f}, {Ln[1]:.6f}, {Ln[2]:.6f}]")

    # Global check: for ALL saturating pairs, is m_perp = 0?
    all_m_parallel = True
    max_m_perp = 0
    for idx in range(len(saturating_pairs)):
        p = saturating_pairs[idx]
        n_vec = p['n_vec']
        m_vec = p['t_vec'] + p['Lambda'] @ n_vec
        m_parallel = np.dot(m_vec, n_vec) * n_vec / np.dot(n_vec, n_vec)
        m_perp = m_vec - m_parallel
        m_perp_norm = np.linalg.norm(m_perp)
        if m_perp_norm > max_m_perp:
            max_m_perp = m_perp_norm
        if m_perp_norm > 1e-7:
            all_m_parallel = False

    print(f"\n  GLOBAL CHECK: m⃗ ∥ n̂ for ALL saturating pairs: {all_m_parallel}")
    print(f"  Max |m⊥|: {max_m_perp:.2e}")

    # Also check converse: for non-saturating but m∥n pairs
    # (This would be checked in Phase 4)


# ============================================================
# BONUS: Analytical verification for depolarizing channel
# ============================================================

def verify_depolarizing_analytical():
    """
    For depolarizing channel N(ρ) = (1-p)ρ + p·I/2:
    - Λ = (1-p)·I, t = 0
    - For ANY pure input |ψ⟩: N(ρ) = (1-p)ρ + p·I/2 lies on ρ-σ line
    - So saturation should hold for ALL states

    Verify analytically:
    - D(ρ‖I/2) = ln 2 for pure ρ
    - N(ρ) has eigenvalues (1±(1-p))/2 = ((2-p)/2, p/2)
    - N(I/2) = I/2
    - D(N(ρ)‖I/2) = ln 2 - H((2-p)/2) where H is binary entropy
    - ΔD = H((2-p)/2)

    Petz map with σ=I/2, N unital (N(I/2)=I/2):
    R(X) = (1/2)·N†(2X) = N†(X) = N(X) since depolarizing is self-adjoint
    So R∘N(ρ) = N(N(ρ)) = (1-p)²ρ + (1-(1-p)²)I/2

    F²(ρ, R∘N(ρ)) = ⟨ψ|R∘N(ρ)|ψ⟩ = (1-p)² + (1-(1-p)²)/2 = (1+(1-p)²)/2

    exp(-ΔD) = exp(-H((2-p)/2)) ... need to compute
    """
    print("\n" + "=" * 70)
    print("ANALYTICAL VERIFICATION: DEPOLARIZING CHANNEL")
    print("=" * 70)

    for p in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        lam = 1 - p  # contraction

        # F² = (1 + λ²) / 2
        F2_analytical = (1 + lam**2) / 2

        # N(ρ) eigenvalues: (1+λ)/2, (1-λ)/2
        e1, e2 = (1+lam)/2, (1-lam)/2

        # D(N(ρ)‖I/2) = ln 2 - S(N(ρ))
        if e1 > 0 and e2 > 0:
            S_Nrho = -e1 * np.log(e1) - e2 * np.log(e2)
        elif e1 > 0:
            S_Nrho = -e1 * np.log(e1)
        else:
            S_Nrho = 0

        D_after = np.log(2) - S_Nrho
        Delta_D = np.log(2) - D_after  # = S(N(ρ))

        exp_neg_DD = np.exp(-Delta_D)
        ratio = F2_analytical / exp_neg_DD if exp_neg_DD > 0 else float('inf')

        # Numerical check
        kraus = make_depolarizing_kraus(p)
        psi = np.array([1, 0], dtype=complex)
        sigma = np.eye(2, dtype=complex) / 2
        data = compute_saturation_data(psi, kraus, sigma)

        if data:
            print(f"  p={p:.1f}: F²_anal={F2_analytical:.8f}, F²_num={data['F2']:.8f}, "
                  f"exp(-ΔD)_anal={exp_neg_DD:.8f}, exp(-ΔD)_num={data['exp_neg_DeltaD']:.8f}, "
                  f"ratio_anal={ratio:.8f}, ratio_num={data['ratio']:.8f}")
        else:
            print(f"  p={p:.1f}: numerical computation failed")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("IFF SATURATION SEARCH")
    print("Finding the exact condition for F²(ρ, R_Petz∘N(ρ)) = exp(-ΔD)")
    print("Reference state: σ = I/2 (maximally mixed)")
    print("Input states: pure (ρ = |ψ⟩⟨ψ|)")
    print("=" * 70)

    t_start = time.time()

    # Phase 0: Analytical check
    verify_depolarizing_analytical()

    # Phase 1: Survey
    saturating_pairs, near_saturating = phase1_survey(
        n_channels=10000, n_states=100, tol=1e-8)

    # Phase 2: Analyze
    properties = phase2_analyze(saturating_pairs)

    # Phase 2.5: Deeper analysis
    phase2_5_deeper_analysis(saturating_pairs, properties)

    # Phase 3: Algebraic
    categories = phase3_algebraic(saturating_pairs, properties)

    # Phase 4: Validate
    accuracy = phase4_validate(saturating_pairs, properties)

    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"Total time: {time.time() - t_start:.1f}s")
    print(f"Saturating pairs found: {len(saturating_pairs)}")
    if properties:
        n = len(saturating_pairs)
        print(f"\nProperties that hold for ALL saturating pairs:")
        bool_props = [
            'is_fixed_point', 'is_unital', 'rho_Nrho_commute', 'rho_recovered_commute',
            'is_adjoint_eigenstate', 'is_zero_DeltaD', 'bloch_parallel',
            'N_rho_on_rho_sigma_line', 'is_isotropic', 'n_eigvec_of_LtL',
            'n_eigvec_of_Lambda', 'n_eigvec_of_LambdaT',
            't_perp_n', 't_parallel_n',
            'all_kraus_commute_with_rho', 'N_rho_is_mixed',
            'perfectly_recovered', 'Lambda_symmetric',
            't_parallel_Ln', 'cond_n_eigvec_and_t_parallel',
        ]
        for prop in bool_props:
            if prop in properties:
                count = sum(properties[prop])
                if count == n:
                    print(f"  ✓ {prop}")

    if accuracy is not None:
        print(f"\nConjecture validation accuracy: {100*accuracy:.4f}%")


if __name__ == "__main__":
    main()
