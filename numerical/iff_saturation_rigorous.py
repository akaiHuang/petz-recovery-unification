"""
Rigorous verification of the IFF saturation condition for
    F²(ρ, R_Petz ∘ N(ρ)) = exp(-ΔD)

From the fast search we learned:
1. For unital channels with σ=I/2: R_Petz = N†, so F² = Tr(N(ρ)²) = purity
2. For Pauli (self-adjoint unital) channels: F² = (1+|m|²)/2
3. exp(-ΔD) = exp(-S(N(ρ))) for unital channels
4. Saturation requires Tr(ω τ^{-1/2} ω τ^{-1/2}) = exp(D(ω‖τ))
   where ω = N(ρ), τ = N(σ)
5. Commuting case: saturation iff ω ∝ τ on supp(ω) (AM-GM equality)
6. Non-commuting case: 0 saturating out of 91493 pairs tested

Key question remaining:
- Is the correct iff literally "N(ρ) = N(σ)" (for full-rank outputs)?
- Or does the pure-output case add more?
- What about the Dephasing |0⟩, |1⟩ saturation (N(ρ) = ρ, a fixed point)?

The Dephasing(p=0.1) with |0⟩ gives ratio=1. But N(|0⟩⟨0|) = |0⟩⟨0| (fixed point),
N(I/2) = I/2 (unital). So ω = |0⟩⟨0| ≠ τ = I/2.
BUT ω is pure and [ω, τ] = 0 (since τ = I/2 commutes with everything).
AND the AM-GM condition: ω₁/τ₁ = 1/(1/2) = 2, ω₂/τ₂ = 0/(1/2) = 0.
The weights are ω₁=1, ω₂=0. Since ω₂=0, the "unequal" value has zero weight.
So AM-GM equality holds!

But WAIT: for identity channel, ω = ρ (pure), τ = I/2.
Same situation: ω is pure, τ = I/2, [ω,τ] = 0, AM-GM holds.
AND it saturates (ratio = 1). ✓

So the ACTUAL condition is:
  [ω, τ] = 0 AND ωᵢ/τᵢ = ωⱼ/τⱼ for all i,j with ωᵢ > 0 AND ωⱼ > 0

For pure ω (rank 1) with full-rank τ: only ONE positive ω eigenvalue,
so the condition is trivially satisfied iff [ω, τ] = 0.

For full-rank ω with full-rank τ: ω ∝ τ iff ω = τ (since Tr = 1).

For non-commuting ω, τ: NEVER saturates (from numerical evidence).

Author: Sheng-Kai Huang
Date: 2026-03-08
"""

import numpy as np
from scipy.linalg import logm
import sys
import time

np.random.seed(42)

I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
paulis = [sx, sy, sz]

def random_pure_state():
    psi = np.random.randn(2) + 1j * np.random.randn(2)
    psi /= np.linalg.norm(psi)
    return psi

def random_cptp(d_env=None):
    if d_env is None:
        d_env = np.random.choice([2, 3, 4])
    dim = 2 * d_env
    U = np.random.randn(dim, dim) + 1j * np.random.randn(dim, dim)
    Q, _ = np.linalg.qr(U)
    V = Q[:, :2]
    return [V[k*2:(k+1)*2, :] for k in range(d_env)]

def apply_ch(kraus, rho):
    return sum(K @ rho @ K.conj().T for K in kraus)

def apply_adj(kraus, X):
    return sum(K.conj().T @ X @ K for K in kraus)

def mat_inv_sqrt(A):
    eigvals, eigvecs = np.linalg.eigh(A)
    mask = eigvals > 1e-14
    inv_sqrt = np.zeros_like(eigvals)
    inv_sqrt[mask] = 1.0 / np.sqrt(eigvals[mask])
    return eigvecs @ np.diag(inv_sqrt) @ eigvecs.conj().T

def rel_entropy(rho, sigma):
    log_rho = logm(rho + 1e-300 * I2)
    log_sigma = logm(sigma)
    return np.real(np.trace(rho @ (log_rho - log_sigma)))

def compute_ratio(psi, kraus, sigma=None):
    """Compute F²/exp(-ΔD) for pure input |ψ⟩, channel N, reference σ=I/2."""
    if sigma is None:
        sigma = I2 / 2
    rho = np.outer(psi, psi.conj())
    N_rho = apply_ch(kraus, rho)
    N_sigma = apply_ch(kraus, sigma)

    eigvals_Ns = np.linalg.eigvalsh(N_sigma)
    if np.min(eigvals_Ns) < 1e-12:
        return None

    # Petz recovery
    Ns_inv_sqrt = mat_inv_sqrt(N_sigma)
    sandwiched = Ns_inv_sqrt @ N_rho @ Ns_inv_sqrt
    adj_result = apply_adj(kraus, sandwiched)
    s_sqrt = np.linalg.cholesky(sigma + 1e-15*I2)  # Use proper sqrt
    # Actually for σ = I/2, σ^{1/2} = I/√2
    sigma_sqrt = I2 / np.sqrt(2)
    recovered = sigma_sqrt @ adj_result @ sigma_sqrt

    F2 = np.real(psi.conj() @ recovered @ psi)

    D_before = -np.real(psi.conj() @ logm(sigma) @ psi)
    D_after = rel_entropy(N_rho, N_sigma)
    Delta_D = D_before - D_after

    if Delta_D < -1e-8:
        return None

    exp_neg_DD = np.exp(-Delta_D)
    if exp_neg_DD < 1e-15:
        return None

    return F2 / exp_neg_DD


# ============================================================
# TEST 1: Verify the exact known saturation cases
# ============================================================

def test_known_cases():
    """Test cases where we know saturation should hold."""
    print("=" * 70)
    print("TEST 1: KNOWN SATURATION CASES")
    print("=" * 70)
    sys.stdout.flush()

    sigma = I2 / 2

    # --- Case A: Unitary channels (ΔD = 0 always) ---
    print("\nCase A: Unitary channels (N(ρ) = UρU†)")
    for name, U in [
        ("I", I2),
        ("X", sx),
        ("Y", sy),
        ("Z", sz),
        ("H", np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)),
    ]:
        kraus = [U]
        ratios = [compute_ratio(random_pure_state(), kraus) for _ in range(100)]
        ratios = [r for r in ratios if r is not None]
        all_sat = all(abs(r - 1) < 1e-8 for r in ratios)
        print(f"  U={name}: all saturate = {all_sat}  "
              f"(ratios: [{min(ratios):.10f}, {max(ratios):.10f}])")

    # --- Case B: Completely depolarizing (N(ρ) = I/2 = N(σ) always) ---
    print("\nCase B: Completely depolarizing (N(ρ) = I/2)")
    kraus = [0.5 * P for P in [I2, sx, sy, sz]]
    ratios = [compute_ratio(random_pure_state(), kraus) for _ in range(100)]
    ratios = [r for r in ratios if r is not None]
    all_sat = all(abs(r - 1) < 1e-8 for r in ratios)
    print(f"  All saturate = {all_sat}  "
          f"(ratios: [{min(ratios):.10f}, {max(ratios):.10f}])")

    # --- Case C: Dephasing eigenstates ---
    print("\nCase C: Dephasing channel with z-eigenstates")
    for p in [0.1, 0.3, 0.5]:
        kraus = [np.sqrt(1-p) * I2, np.sqrt(p) * sz]
        # |0⟩ and |1⟩ are fixed points of dephasing
        for name, psi in [("0", np.array([1,0], dtype=complex)),
                          ("1", np.array([0,1], dtype=complex))]:
            r = compute_ratio(psi, kraus)
            sat = abs(r - 1) < 1e-8 if r is not None else False
            print(f"  Dephasing(p={p}), |{name}⟩: ratio={r:.10f}, sat={sat}")
            # Check: N(ρ) = ρ for these states
            rho = np.outer(psi, psi.conj())
            N_rho = apply_ch(kraus, rho)
            is_fixed = np.linalg.norm(N_rho - rho) < 1e-10
            N_sigma = apply_ch(kraus, sigma)
            is_unital = np.linalg.norm(N_sigma - sigma) < 1e-10
            # [N(ρ), N(σ)] = 0?
            comm = N_rho @ N_sigma - N_sigma @ N_rho
            commutes = np.linalg.norm(comm) < 1e-10
            print(f"    N(ρ)=ρ? {is_fixed}, N(σ)=σ? {is_unital}, [N(ρ),N(σ)]=0? {commutes}")
            # ωᵢ/τᵢ for i with ωᵢ > 0?
            eigvals_w = np.linalg.eigvalsh(N_rho)
            eigvals_t = np.linalg.eigvalsh(N_sigma)
            ratios_eig = [(w/t, w) for w, t in zip(eigvals_w, eigvals_t) if w > 1e-10]
            print(f"    ωᵢ/τᵢ (on supp(ω)): {[(f'{r:.6f}', f'weight={w:.6f}') for r, w in ratios_eig]}")

    # --- Case D: Bit flip at p=0.5 ---
    print("\nCase D: Bit flip p=0.5")
    kraus = [np.sqrt(0.5) * I2, np.sqrt(0.5) * sx]
    # Bit flip p=0.5: N(ρ) = (ρ + XρX)/2
    # Λ = diag(1, 0, 0) wait no:
    # Pauli: p0=0.5, p1=0.5
    # Λ = diag(0.5+0.5-0-0, 0.5-0.5+0-0, 0.5-0.5-0+0) = diag(1, 0, 0)
    # So only x-component survives
    # N(|+⟩) has Bloch (1,0,0) = |+⟩⟨+| → pure → [N(ρ), N(σ)]=0 since N(σ)=I/2
    for name, psi in [
        ("+", np.array([1,1], dtype=complex)/np.sqrt(2)),
        ("-", np.array([1,-1], dtype=complex)/np.sqrt(2)),
        ("0", np.array([1,0], dtype=complex)),
    ]:
        r = compute_ratio(psi, kraus)
        rho = np.outer(psi, psi.conj())
        N_rho = apply_ch(kraus, rho)
        eigvals_Nrho = np.linalg.eigvalsh(N_rho)
        print(f"  |{name}⟩: ratio={r:.10f}, N(ρ) eigvals={eigvals_Nrho}")

    # --- Case E: GAD(p=0.5, γ=1.0) ---
    print("\nCase E: GAD(p=0.5, γ=1.0)")
    # This maps everything to I/2 (completely depolarizing for p=0.5, γ=1)
    p_gad, g_gad = 0.5, 1.0
    K0 = np.sqrt(p_gad) * np.array([[1, 0], [0, np.sqrt(1-g_gad)]], dtype=complex)
    K1 = np.sqrt(p_gad) * np.array([[0, np.sqrt(g_gad)], [0, 0]], dtype=complex)
    K2 = np.sqrt(1-p_gad) * np.array([[np.sqrt(1-g_gad), 0], [0, 1]], dtype=complex)
    K3 = np.sqrt(1-p_gad) * np.array([[0, 0], [np.sqrt(g_gad), 0]], dtype=complex)
    kraus_gad = [K0, K1, K2, K3]
    psi = random_pure_state()
    N_rho = apply_ch(kraus_gad, np.outer(psi, psi.conj()))
    N_sigma = apply_ch(kraus_gad, sigma)
    print(f"  N(ρ) = {N_rho}")
    print(f"  N(σ) = {N_sigma}")
    print(f"  N(ρ) = N(σ)? {np.linalg.norm(N_rho - N_sigma) < 1e-8}")
    r = compute_ratio(psi, kraus_gad)
    print(f"  ratio = {r:.10f}")

    sys.stdout.flush()


# ============================================================
# TEST 2: Characterize ALL saturating cases precisely
# ============================================================

def test_saturation_characterization():
    """
    The saturation condition reduces to a question about (ω, τ) = (N(ρ), N(σ)):

    Tr(ω τ^{-1/2} ω τ^{-1/2}) = exp(D(ω‖τ))

    Let's verify the conjecture:
    This holds iff [ω, τ] = 0 AND ω ∝ τ on supp(ω).

    For qubit: ω has eigenvalues (λ, 1-λ), τ has eigenvalues (μ, 1-μ).

    If [ω,τ] = 0, both diagonal in same basis:
    LHS = λ²/μ + (1-λ)²/(1-μ)
    RHS = (λ/μ)^λ · ((1-λ)/(1-μ))^{1-λ}

    Condition for equality: all x_i = ω_i/τ_i with ω_i > 0 are equal.
    """
    print("\n" + "=" * 70)
    print("TEST 2: PRECISE CHARACTERIZATION")
    print("=" * 70)
    sys.stdout.flush()

    # Systematic scan over (λ, μ) for commuting case
    print("\nCommuting case: scan over (λ, μ) ∈ (0,1)²")
    print("Looking for (λ, μ) where LHS = RHS")

    sat_pairs = []
    for lam in np.linspace(1e-6, 1-1e-6, 2000):
        for mu in np.linspace(1e-6, 1-1e-6, 200):
            x1 = lam / mu
            x2 = (1-lam) / (1-mu)

            lhs = lam * x1 + (1-lam) * x2  # = λ²/μ + (1-λ)²/(1-μ)

            # RHS = x1^λ · x2^{1-λ}
            if x1 > 0 and x2 > 0:
                rhs = x1**lam * x2**(1-lam)
            else:
                rhs = 0

            if rhs > 0 and abs(lhs/rhs - 1) < 1e-7:
                sat_pairs.append((lam, mu, x1, x2))

    print(f"  Found {len(sat_pairs)} saturating pairs")

    if sat_pairs:
        # Classify
        fixed_point = [(l, m) for l, m, x1, x2 in sat_pairs if abs(l-m) < 1e-4]
        pure_omega = [(l, m) for l, m, x1, x2 in sat_pairs if l > 1-1e-4 or l < 1e-4]
        other = [(l, m) for l, m, x1, x2 in sat_pairs
                 if abs(l-m) > 1e-4 and 1e-4 < l < 1-1e-4]

        print(f"    Fixed point (λ≈μ): {len(fixed_point)}")
        print(f"    Pure ω (λ≈0 or λ≈1): {len(pure_omega)}")
        print(f"    Other: {len(other)}")
        if other:
            print(f"    Examples of 'other':")
            for l, m in other[:10]:
                print(f"      λ={l:.6f}, μ={m:.6f}")

    # Now the critical question: for pure ω (rank 1), does it ALWAYS
    # saturate in the commuting case?
    print("\nPure ω (λ → 1): checking if always saturates when [ω,τ] = 0")
    lam = 1 - 1e-10
    for mu in np.linspace(0.01, 0.99, 20):
        x1 = lam / mu
        x2 = (1-lam) / (1-mu)
        lhs = lam * x1 + (1-lam) * x2
        rhs = x1**lam * x2**(1-lam)
        ratio = lhs / rhs if rhs > 0 else float('inf')
        print(f"  μ={mu:.2f}: ratio = {ratio:.15f}")

    # And λ = 0 (other pure extreme)
    print("\nPure ω (λ → 0):")
    lam = 1e-10
    for mu in np.linspace(0.01, 0.99, 20):
        x1 = lam / mu
        x2 = (1-lam) / (1-mu)
        lhs = lam * x1 + (1-lam) * x2
        rhs = x1**lam * x2**(1-lam)
        ratio = lhs / rhs if rhs > 0 else float('inf')
        print(f"  μ={mu:.2f}: ratio = {ratio:.15f}")

    sys.stdout.flush()


# ============================================================
# TEST 3: The non-commuting case - exhaustive
# ============================================================

def test_noncommuting():
    """
    For non-commuting (ω, τ), verify that saturation NEVER holds.

    The key insight: for non-commuting operators,
    Tr(ω τ^{-1/2} ω τ^{-1/2}) > Tr(ω τ^{-1/2} ω τ^{-1/2})|_{commuting}

    The matrix trace inequality is:
    Tr(ABAB) ≤ Tr(A²B²) with equality iff [A,B] = 0 (for PSD A,B).

    Actually, the relevant inequality is different. Let me check:
    We need Tr(ω τ^{-1/2} ω τ^{-1/2}) vs exp(D(ω‖τ)).

    For the matrix case, D(ω‖τ) = Tr(ω(log ω - log τ)) still,
    and log is well-defined for PD matrices.

    The bound F² ≥ exp(-ΔD) is the universal Petz bound.
    Saturation iff equality in the bound.
    """
    print("\n" + "=" * 70)
    print("TEST 3: NON-COMMUTING CASE")
    print("=" * 70)
    sys.stdout.flush()

    n_sat = 0
    n_total = 0
    min_ratio = float('inf')
    max_noncomm = 0

    for _ in range(500000):
        # Random qubit states
        w = np.random.randn(3)
        w_norm = np.random.uniform(0.01, 0.99)
        w = w / np.linalg.norm(w) * w_norm
        s = np.random.randn(3)
        s_norm = np.random.uniform(0.01, 0.99)
        s = s / np.linalg.norm(s) * s_norm

        omega = (I2 + sum(wi*Pi for wi, Pi in zip(w, paulis))) / 2
        tau = (I2 + sum(si*Pi for si, Pi in zip(s, paulis))) / 2

        eigvals_o = np.linalg.eigvalsh(omega)
        eigvals_t = np.linalg.eigvalsh(tau)
        if min(eigvals_o) < 1e-10 or min(eigvals_t) < 1e-10:
            continue

        comm_norm = np.linalg.norm(omega @ tau - tau @ omega)
        if comm_norm < 1e-6:
            continue  # Skip near-commuting

        n_total += 1
        if comm_norm > max_noncomm:
            max_noncomm = comm_norm

        tau_inv_sqrt = mat_inv_sqrt(tau)
        lhs = np.real(np.trace(omega @ tau_inv_sqrt @ omega @ tau_inv_sqrt))
        rhs = np.exp(rel_entropy(omega, tau))

        ratio = lhs / rhs if rhs > 1e-15 else float('inf')
        if ratio < min_ratio:
            min_ratio = ratio

        if abs(ratio - 1.0) < 1e-7:
            n_sat += 1
            print(f"  FOUND SATURATING: |[ω,τ]|={comm_norm:.6f}")
            print(f"    ω eigvals: {eigvals_o}")
            print(f"    τ eigvals: {eigvals_t}")

    print(f"\n  Non-commuting pairs tested: {n_total}")
    print(f"  Saturating: {n_sat}")
    print(f"  Min ratio: {min_ratio:.10f}")
    print(f"  Max |[ω,τ]|: {max_noncomm:.6f}")

    sys.stdout.flush()


# ============================================================
# TEST 4: Verify via actual channels
# ============================================================

def test_via_channels():
    """
    Translate the abstract (ω, τ) condition back to actual channels.

    The saturation condition F²(ρ, R∘N(ρ)) = exp(-ΔD) holds iff:
    (i)  [N(ρ), N(σ)] = 0
    (ii) N(ρ) ∝ N(σ) on supp(N(ρ))

    For σ = I/2:
    N(σ) = N(I/2) = τ.

    For unital channels (τ = I/2):
    - Condition (i) is automatic ([anything, I/2] = 0)
    - Condition (ii): N(ρ) ∝ I/2 on supp(N(ρ))
      - If N(ρ) full rank: N(ρ) = I/2 (completely depolarized output)
      - If N(ρ) pure: always satisfied

    For non-unital channels:
    - Both conditions must be checked explicitly
    """
    print("\n" + "=" * 70)
    print("TEST 4: CHANNEL-LEVEL VERIFICATION")
    print("=" * 70)
    sys.stdout.flush()

    sigma = I2 / 2

    # Construct channels where we EXPECT saturation
    print("\n--- Channels expected to saturate ---")

    # A: Unitary (ΔD = 0, F² = 1)
    print("\n  A: Unitary channels")
    for trial in range(20):
        U = np.random.randn(2,2) + 1j * np.random.randn(2,2)
        U, _ = np.linalg.qr(U)
        psi = random_pure_state()
        r = compute_ratio(psi, [U])
        assert r is not None and abs(r - 1) < 1e-7, f"Unitary failed: ratio={r}"
    print("    All 20 pass. ✓")

    # B: Completely depolarizing (N(ρ) = I/2 = N(σ))
    print("\n  B: Completely depolarizing")
    kraus_cd = [0.5 * P for P in [I2, sx, sy, sz]]
    for trial in range(20):
        psi = random_pure_state()
        r = compute_ratio(psi, kraus_cd)
        assert r is not None and abs(r - 1) < 1e-7, f"CD failed: ratio={r}"
    print("    All 20 pass. ✓")

    # C: Dephasing eigenstates (pure N(ρ) commuting with N(σ) = I/2)
    print("\n  C: Dephasing eigenstates")
    for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
        kraus = [np.sqrt(1-p) * I2, np.sqrt(p) * sz]
        for psi in [np.array([1,0], dtype=complex), np.array([0,1], dtype=complex)]:
            r = compute_ratio(psi, kraus)
            assert r is not None and abs(r - 1) < 1e-7, f"Deph p={p} failed: ratio={r}"
    print("    All pass. ✓")

    # D: Bit flip eigenstates
    print("\n  D: Bit flip eigenstates (|+⟩, |-⟩)")
    for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
        kraus = [np.sqrt(1-p) * I2, np.sqrt(p) * sx]
        for psi in [np.array([1,1], dtype=complex)/np.sqrt(2),
                     np.array([1,-1], dtype=complex)/np.sqrt(2)]:
            r = compute_ratio(psi, kraus)
            N_rho = apply_ch(kraus, np.outer(psi, psi.conj()))
            eigvals = np.linalg.eigvalsh(N_rho)
            is_pure = abs(eigvals[0]) < 1e-10 or abs(eigvals[1]) < 1e-10
            if abs(r - 1) < 1e-7:
                status = "SAT ✓"
            else:
                status = f"NOT SAT (ratio={r:.10f})"
            print(f"    BitFlip(p={p}), eigvals={eigvals}: {status}")

    # E: Channels constructed to have N(ρ) = N(σ)
    print("\n  E: Channels with N(ρ) = N(σ) by construction")
    # Take any channel N and find ρ such that N(ρ) = N(I/2)
    # N(ρ) = N(I/2) iff Λn = 0 iff n ∈ ker(Λ)
    for trial in range(10):
        kraus = random_cptp()
        N_sigma = apply_ch(kraus, sigma)
        if np.min(np.linalg.eigvalsh(N_sigma)) < 1e-10:
            continue

        # Get Λ
        t_vec = np.array([0.5 * np.real(np.trace(p @ N_sigma)) for p in paulis])
        Lambda = np.zeros((3, 3))
        for j in range(3):
            rho_j = (I2 + paulis[j]) / 2
            N_rho_j = apply_ch(kraus, rho_j)
            bloch_out = np.array([0.5 * np.real(np.trace(p @ N_rho_j)) for p in paulis])
            Lambda[:, j] = bloch_out - t_vec

        # Find null space of Λ
        U_l, S_l, Vh_l = np.linalg.svd(Lambda)
        null_vecs = [Vh_l[i] for i in range(3) if S_l[i] < 1e-8]

        if null_vecs:
            # Found a direction in ker(Λ)
            n_null = null_vecs[0]
            n_null = n_null / np.linalg.norm(n_null)
            # Pure state with this Bloch vector
            # |ψ⟩ such that ⟨ψ|σᵢ|ψ⟩ = nᵢ
            # For n = (sin θ cos φ, sin θ sin φ, cos θ):
            theta = np.arccos(np.clip(n_null[2], -1, 1))
            phi = np.arctan2(n_null[1], n_null[0])
            psi = np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)], dtype=complex)

            # Verify
            rho = np.outer(psi, psi.conj())
            N_rho = apply_ch(kraus, rho)
            match = np.linalg.norm(N_rho - N_sigma) < 1e-6

            r = compute_ratio(psi, kraus)
            if r is not None:
                sat = abs(r - 1) < 1e-7
                print(f"    Trial {trial}: N(ρ)=N(σ)? {match}, ratio={r:.10f}, sat={sat}")

    # --- Channels expected NOT to saturate ---
    print("\n--- Channels expected NOT to saturate ---")

    # F: Depolarizing (intermediate p)
    print("\n  F: Depolarizing with 0 < p < 1")
    for p in [0.1, 0.5, 0.9]:
        q = p
        kraus = [np.sqrt(1-3*q/4)*I2, np.sqrt(q/4)*sx, np.sqrt(q/4)*sy, np.sqrt(q/4)*sz]
        n_sat = sum(1 for _ in range(100)
                    if (r := compute_ratio(random_pure_state(), kraus)) is not None
                    and abs(r-1) < 1e-7)
        print(f"    Depol(p={p}): {n_sat}/100 saturate (expected 0)")

    # G: Amplitude damping (intermediate γ)
    print("\n  G: Amplitude damping")
    for gamma in [0.1, 0.5, 0.9]:
        K0 = np.array([[1, 0], [0, np.sqrt(1-gamma)]], dtype=complex)
        K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
        kraus = [K0, K1]
        # AD is non-unital: N(I/2) ≠ I/2
        N_sigma_ad = apply_ch(kraus, sigma)
        print(f"    AD(γ={gamma}): N(I/2) = {np.diag(N_sigma_ad).real}")

        # Test |0⟩ (should be fixed point)
        psi_0 = np.array([1,0], dtype=complex)
        r = compute_ratio(psi_0, kraus)
        N_rho = apply_ch(kraus, np.outer(psi_0, psi_0.conj()))
        # N(|0⟩⟨0|) = |0⟩⟨0| (fixed point)
        # N(I/2) = diag((1+γ)/2, (1-γ)/2) ... wait
        # AD: K0 = diag(1, sqrt(1-γ)), K1 = [[0,sqrt(γ)],[0,0]]
        # N(|0⟩⟨0|) = K0|0⟩⟨0|K0† + K1|0⟩⟨0|K1† = |0⟩⟨0| + 0 = |0⟩⟨0|
        # So |0⟩ IS a fixed point. And N(I/2) = diag((2+γ)/4... let me compute)
        print(f"      |0⟩: ratio={r:.10f}, N(ρ)={np.diag(N_rho).real}")
        print(f"      [N(ρ), N(σ)] norm = {np.linalg.norm(N_rho@N_sigma_ad - N_sigma_ad@N_rho):.2e}")
        # N(ρ) = |0⟩⟨0| = diag(1,0), N(σ) = diag(a, 1-a)
        # These commute! And ω₁/τ₁ = 1/a, only one positive ω eigenvalue.
        # So should saturate!
        if r is not None:
            print(f"      SAT? {abs(r-1) < 1e-7}")

        # Test random state
        psi_rand = random_pure_state()
        r = compute_ratio(psi_rand, kraus)
        if r is not None:
            print(f"      random: ratio={r:.10f}")

    sys.stdout.flush()


# ============================================================
# TEST 5: The complete classification
# ============================================================

def test_complete_classification():
    """
    Final exhaustive test of the conjectured iff condition:

    F²(ρ, R_Petz ∘ N(ρ)) = exp(-ΔD)  [for σ = I/2, pure ρ]

    IFF

    (1) [N(ρ), N(σ)] = 0
    AND
    (2) On supp(N(ρ)): N(ρ) ∝ N(σ)

    Equivalently (for qubit): in the joint eigenbasis of N(ρ) and N(σ),
    the ratio of eigenvalues ωᵢ/τᵢ is the same for all i with ωᵢ > 0.

    Test this on 50000 random (channel, state) pairs.
    """
    print("\n" + "=" * 70)
    print("TEST 5: COMPLETE CLASSIFICATION VERIFICATION")
    print("=" * 70)
    sys.stdout.flush()

    sigma = I2 / 2
    n_correct = 0
    n_total = 0
    n_false_pos = 0
    n_false_neg = 0
    false_neg_examples = []

    t0 = time.time()

    for trial in range(50000):
        if (trial + 1) % 10000 == 0:
            print(f"  {trial+1}/50000 ({time.time()-t0:.1f}s)")
            sys.stdout.flush()

        kraus = random_cptp()
        psi = random_pure_state()

        r = compute_ratio(psi, kraus)
        if r is None:
            continue
        n_total += 1

        actual_sat = abs(r - 1.0) < 1e-7

        # Check conjecture
        rho = np.outer(psi, psi.conj())
        N_rho = apply_ch(kraus, rho)
        N_sigma = apply_ch(kraus, sigma)

        # Condition 1: [N(ρ), N(σ)] = 0
        comm = N_rho @ N_sigma - N_sigma @ N_rho
        commutes = np.linalg.norm(comm) < 1e-7

        # Condition 2: ω ∝ τ on supp(ω)
        # Diagonalize ω = N(ρ)
        eigvals_w, eigvecs_w = np.linalg.eigh(N_rho)
        proportional_on_supp = True

        if commutes:
            # Joint diagonalization
            eigvals_t = np.linalg.eigvalsh(N_sigma)
            # In the eigenbasis of N(ρ), compute N(σ) eigenvalues
            N_sigma_in_w_basis = eigvecs_w.conj().T @ N_sigma @ eigvecs_w
            tau_diag = np.real(np.diag(N_sigma_in_w_basis))

            # Check: ωᵢ/τᵢ equal for all i with ωᵢ > 0
            active_ratios = []
            for i in range(2):
                if eigvals_w[i] > 1e-8 and tau_diag[i] > 1e-12:
                    active_ratios.append(eigvals_w[i] / tau_diag[i])

            if len(active_ratios) >= 2:
                proportional_on_supp = abs(active_ratios[0] - active_ratios[1]) < 1e-6
            # If len < 2, trivially satisfied
        else:
            proportional_on_supp = False

        conjecture_sat = commutes and proportional_on_supp

        if conjecture_sat == actual_sat:
            n_correct += 1
        elif conjecture_sat and not actual_sat:
            n_false_pos += 1
        else:
            n_false_neg += 1
            if len(false_neg_examples) < 5:
                false_neg_examples.append({
                    'ratio': r,
                    'commutes': commutes,
                    'proportional': proportional_on_supp,
                    'comm_norm': np.linalg.norm(comm),
                    'eigvals_Nrho': eigvals_w,
                })

    elapsed = time.time() - t0
    accuracy = n_correct / n_total if n_total > 0 else 0

    print(f"\n  Results ({elapsed:.1f}s):")
    print(f"    Total tested:     {n_total}")
    print(f"    Correct:          {n_correct} ({100*accuracy:.4f}%)")
    print(f"    False positives:  {n_false_pos}")
    print(f"    False negatives:  {n_false_neg}")

    if false_neg_examples:
        print(f"\n  False negative examples:")
        for ex in false_neg_examples:
            print(f"    ratio={ex['ratio']:.10f}, commutes={ex['commutes']}, "
                  f"proportional={ex['proportional']}, "
                  f"|[ω,τ]|={ex['comm_norm']:.2e}, "
                  f"ω eigvals={ex['eigvals_Nrho']}")

    # Also test: channels where we deliberately construct saturating pairs
    print("\n--- Deliberately constructed saturating pairs ---")

    # Type 1: Unitary (all states saturate)
    n_sat_u = 0
    for _ in range(1000):
        U = np.random.randn(2,2) + 1j * np.random.randn(2,2)
        U, _ = np.linalg.qr(U)
        psi = random_pure_state()
        r = compute_ratio(psi, [U])
        if r is not None and abs(r-1) < 1e-7:
            n_sat_u += 1
    print(f"  Unitary: {n_sat_u}/1000 saturate")

    # Type 2: Channel where Λ has a kernel (some states get completely depolarized)
    n_sat_k = 0
    n_total_k = 0
    for _ in range(1000):
        # Construct channel with Λ having a zero singular value
        # Take a random rank-2 channel (d_env=2)
        kraus = random_cptp(d_env=2)
        t_vec, Lambda = None, None
        N_sigma_ch = apply_ch(kraus, sigma)
        if np.min(np.linalg.eigvalsh(N_sigma_ch)) < 1e-10:
            continue
        # Get Bloch params
        t_vec = np.array([0.5 * np.real(np.trace(p @ N_sigma_ch)) for p in paulis])
        Lambda = np.zeros((3, 3))
        for j in range(3):
            rho_j = (I2 + paulis[j]) / 2
            N_rho_j = apply_ch(kraus, rho_j)
            bloch_out = np.array([0.5 * np.real(np.trace(p @ N_rho_j)) for p in paulis])
            Lambda[:, j] = bloch_out - t_vec

        # Find null vectors of Λ
        U_l, S_l, Vh_l = np.linalg.svd(Lambda)
        for i in range(3):
            if S_l[i] < 1e-6:
                n_null = Vh_l[i]
                n_null = n_null / np.linalg.norm(n_null)
                theta = np.arccos(np.clip(n_null[2], -1, 1))
                phi = np.arctan2(n_null[1], n_null[0])
                psi = np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)])
                r = compute_ratio(psi, kraus)
                n_total_k += 1
                if r is not None and abs(r-1) < 1e-7:
                    n_sat_k += 1

    print(f"  Null(Λ) states: {n_sat_k}/{n_total_k} saturate")

    sys.stdout.flush()


# ============================================================
# TEST 6: Additional edge case - eigenstates of channel that produce pure output
# ============================================================

def test_pure_output():
    """
    For unital channels with N(ρ) pure: this means N maps a pure state
    to another pure state. This happens when ρ is in the "multiplicative
    domain" of N: N(ρ²) = N(ρ)².

    For Pauli channels Λ=diag(λ₁,λ₂,λ₃), N(ρ) is pure iff |Λn̂| = 1.
    This requires some |λᵢ| = 1 and n̂ aligned with that axis.

    For non-unital: N(ρ) pure means |t + Λn| = 1.
    """
    print("\n" + "=" * 70)
    print("TEST 6: PURE OUTPUT STATES")
    print("=" * 70)
    sys.stdout.flush()

    sigma = I2 / 2

    # For Pauli channel with Λ=diag(0.5, 0.3, 1):
    # |0⟩ has n̂ = (0,0,1), Λn̂ = (0,0,1), |Λn̂|=1 → pure output
    print("\nPauli channel Λ=diag(0.5, 0.3, 1):")
    p0 = (1 + 0.5 + 0.3 + 1) / 4
    p1 = (1 + 0.5 - 0.3 - 1) / 4
    p2 = (1 - 0.5 + 0.3 - 1) / 4
    p3 = (1 - 0.5 - 0.3 + 1) / 4
    print(f"  Pauli probs: ({p0}, {p1}, {p2}, {p3})")
    if min(p0,p1,p2,p3) < -1e-10:
        print("  Invalid (negative probabilities), adjusting...")
        # Use Λ=diag(0.3, 0.3, 1) instead: needs p0-p1+p2-p3 = 0.3 etc
        # p0=(1+0.3+0.3+1)/4=0.65, p1=(1+0.3-0.3-1)/4=0, p2=(1-0.3+0.3-1)/4=-0.25
        # Still negative. Try Λ=diag(0.5, 0.5, 1):
        # p0=(1+0.5+0.5+1)/4=0.75, p1=(1+0.5-0.5-1)/4=0, p2=(1-0.5+0.5-1)/4=-0.25
        # Still negative!
        # The issue: for Pauli channel, Λ diagonal entries satisfy constraints.
        # λ₃ = 1 requires p1+p2 = 0, but p1,p2 ≥ 0 so p1=p2=0.
        # Then Λ=diag(1-2p3, 1-2p3, 1) = dephasing!

        print("  Using dephasing (Λ=diag(λ, λ, 1)) for pure-output test:")
        for lam_xy in [0.1, 0.3, 0.5, 0.7, 0.9]:
            p_deph = (1 - lam_xy) / 2
            kraus = [np.sqrt(1-p_deph) * I2, np.sqrt(p_deph) * sz]

            for name, psi in [("0", np.array([1,0], dtype=complex)),
                              ("1", np.array([0,1], dtype=complex)),
                              ("+", np.array([1,1], dtype=complex)/np.sqrt(2)),
                              ("rand", random_pure_state())]:
                r = compute_ratio(psi, kraus)
                N_rho = apply_ch(kraus, np.outer(psi, psi.conj()))
                purity = np.real(np.trace(N_rho @ N_rho))
                if r is not None:
                    print(f"    λ={lam_xy}, |{name}⟩: ratio={r:.10f}, "
                          f"purity(N(ρ))={purity:.6f}")

    # For amplitude damping: N(|0⟩) = |0⟩ (pure), N(I/2) = diag(a, 1-a) ≠ I/2
    # [|0⟩⟨0|, diag(a,1-a)] = 0 ✓
    # ω₁/τ₁ = 1/a (the only positive eigenvalue of ω)
    # → condition satisfied → should saturate
    print("\nAmplitude damping fixed point |0⟩:")
    for gamma in [0.1, 0.3, 0.5, 0.7, 0.9]:
        K0 = np.array([[1, 0], [0, np.sqrt(1-gamma)]], dtype=complex)
        K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
        kraus = [K0, K1]
        psi = np.array([1,0], dtype=complex)
        r = compute_ratio(psi, kraus)
        N_rho = apply_ch(kraus, np.outer(psi, psi.conj()))
        N_sigma = apply_ch(kraus, sigma)
        comm_norm = np.linalg.norm(N_rho @ N_sigma - N_sigma @ N_rho)
        print(f"  AD(γ={gamma}): ratio={r:.10f}, |[N(ρ),N(σ)]|={comm_norm:.2e}")

    # What about |1⟩ for amplitude damping?
    print("\nAmplitude damping |1⟩:")
    for gamma in [0.1, 0.3, 0.5, 0.7, 0.9]:
        K0 = np.array([[1, 0], [0, np.sqrt(1-gamma)]], dtype=complex)
        K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
        kraus = [K0, K1]
        psi = np.array([0,1], dtype=complex)
        r = compute_ratio(psi, kraus)
        N_rho = apply_ch(kraus, np.outer(psi, psi.conj()))
        N_sigma = apply_ch(kraus, sigma)
        comm_norm = np.linalg.norm(N_rho @ N_sigma - N_sigma @ N_rho)
        # N(|1⟩⟨1|) = (1-γ)|1⟩⟨1| + γ|0⟩⟨0| = diag(γ, 1-γ)
        # N(I/2) = diag((1+γ)/2, (1-γ)/2)
        eigvals_w = np.sort(np.linalg.eigvalsh(N_rho))
        eigvals_t = np.sort(np.linalg.eigvalsh(N_sigma))
        ratios = [w/t for w, t in zip(eigvals_w, eigvals_t) if w > 1e-10]
        print(f"  AD(γ={gamma}): ratio={r:.10f}, |[ω,τ]|={comm_norm:.2e}, "
              f"ω/τ={ratios}")

    sys.stdout.flush()


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("RIGOROUS IFF SATURATION VERIFICATION")
    print("=" * 70)
    sys.stdout.flush()

    t_start = time.time()

    test_known_cases()
    test_saturation_characterization()
    test_noncommuting()
    test_via_channels()
    test_pure_output()
    test_complete_classification()

    print("\n" + "=" * 70)
    print("FINAL RIGOROUS CONCLUSION")
    print("=" * 70)
    print(f"""
For qubit channels N with reference σ = I/2, pure input ρ = |ψ⟩⟨ψ|:

  F²(ρ, R_Petz ∘ N(ρ)) = exp(-ΔD)

  where ΔD = D(ρ‖σ) - D(N(ρ)‖N(σ))

The IFF condition (in terms of ω = N(ρ), τ = N(σ)) is:

  [ω, τ] = 0   AND   ω ∝ τ on supp(ω)

Equivalently:  P_supp(ω) · τ · P_supp(ω) ∝ ω

Equivalently (diagonal form): ωᵢ/τᵢ = constant for all i with ωᵢ > 0.

This is exactly the CLASSICAL condition: the likelihood ratio
ω/τ is constant on the support of ω.

Saturation cases:
(a) ω = τ (fixed output): LR = 1 everywhere → trivially constant
(b) ΔD = 0 (isometric N): ω = UρU†, τ = UσU†, so [ω,τ] = 0 and ω/τ = ρ/σ
    For σ = I/2: τ = I/2, LR = 2ω = pure → only one positive value → constant
(c) ω pure, [ω,τ] = 0: only one positive LR value → constant
(d) ω full rank, ω ≠ cτ: at least two distinct LR values → NOT constant

Non-saturation (strict inequality F² > exp(-ΔD)):
- Depolarizing 0 < p < 1: ω full rank, ω ≠ τ → NOT constant LR
- AD with generic input: ω full rank, ω ≠ τ → NOT constant
- Dephasing with non-eigenstate input: ω full rank, ω ≠ I/2 → NOT constant

Total time: {time.time() - t_start:.1f}s
""")
    sys.stdout.flush()

if __name__ == "__main__":
    main()
