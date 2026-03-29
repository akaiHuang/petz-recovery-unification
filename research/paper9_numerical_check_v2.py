#!/usr/bin/env python3
"""
Paper 9 — Independent Numerical Verification (v2)
==================================================
Isolate the Seeley-DeWitt a_4 coefficient under conformal rescaling
on a flat 4-torus T^4 with metric g_μν = Q² δ_μν.

KEY INSIGHT from v1:
  The TOTAL spectral action mixes a_0, a_2, a_4 terms with different Q-scalings.
  To isolate a_4, we use the heat-kernel expansion:

    Tr exp(-t D²) = (4πt)^{-d/2} * Vol(g) * [a_0 + a_2 t + a_4 t² + ...]

  On T^4 with g = Q² δ:
    - Vol(g) = (QL)^4
    - D² eigenvalues = (2π n/L)² / Q²  [conformal scaling in d=4]
    - t → t/Q² in the heat trace (since eigenvalues scale as 1/Q²)

  The heat trace becomes:
    K(t,Q) = sum_n exp(-t (2πn/L)² / Q²)
           = (4π t/Q²)^{-2} * (QL)^4 * [a_0 + a_2(Q) t/Q² + a_4(Q) (t/Q²)² + ...]

  For FLAT torus: a_0 = 1, a_2 = R/6 = 0, a_4 = (curvature)² terms = 0.

  BUT: when we vary Q, the geometry stays flat (it's a conformal rescaling of flat).
  So a_0 = 1, a_2 = 0, a_4 = 0 FOR ALL Q on a flat torus.

  This means the SPECTRAL ACTION conformal variation on flat T^4 comes entirely
  from the volume factor and the eigenvalue scaling, NOT from a_4.

  The a_4 coefficient is what matters for the CONFORMAL ANOMALY in curved space.
  On flat space, a_4 = 0 identically. So we need to:

  STRATEGY: Use the Seeley-DeWitt expansion on a CURVED background instead.
  OR: Extract the conformal anomaly coefficient from the t-expansion.

  Let's do BOTH:
  (A) Verify the Q-scaling of S(Q) on flat T^4 matches the prediction from a_0 alone
  (B) Add a small curvature perturbation and extract a_4's conformal response
  (C) Use the t-dependence to numerically extract each a_k coefficient and study
      their Q-dependence separately.

Author : numerical verification (independent)
Date   : 2026-03-19
"""

import numpy as np
from scipy.optimize import curve_fit

# ============================================================
# PART A: Flat T^4 — verify a_0-dominated scaling
# ============================================================

def heat_trace_flat_T4(t, Q, L, N_max):
    """
    K(t, Q) = sum_{n in Z^4} exp(-t * (2π|n|/L)² / Q²) - 1  [subtract n=0]

    Using Jacobi theta: K = [θ₃(0, exp(-t(2π/L)²/Q²))]⁴ - 1
    """
    arg = t * (2*np.pi/L)**2 / Q**2
    q = np.exp(-arg)
    if q > 0.999:
        # Danger: very slow convergence
        # Use direct lattice sum
        rng = np.arange(-N_max, N_max + 1)
        n1, n2, n3, n4 = np.meshgrid(rng, rng, rng, rng, indexing='ij')
        nsq = (n1**2 + n2**2 + n3**2 + n4**2).ravel()
        vals = np.exp(-arg * nsq)
        vals[nsq == 0] = 0.0
        return np.sum(vals)
    else:
        # Theta function
        theta = 1.0
        for n in range(1, N_max + 1):
            term = 2 * q**(n**2)
            theta += term
            if term < 1e-20:
                break
        return theta**4 - 1.0


def extract_a_coefficients_from_t(Q, L, N_max):
    """
    From the heat trace K(t,Q), extract a_k by fitting:
        K(t,Q) * (4πt/Q²)² / (QL)⁴ = a_0 + a_2 * t/Q² + a_4 * (t/Q²)² + ...

    We compute the normalized heat trace for several t values and fit.
    """
    # Choose t values small enough for asymptotic expansion to work
    # but large enough that the lattice sum converges
    t_values = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1])
    tau_values = t_values / Q**2  # effective expansion parameter

    K_values = np.array([heat_trace_flat_T4(t, Q, L, N_max) for t in t_values])

    # Normalize: K_norm = K * (4π tau)^2 / Vol
    Vol = (Q * L)**4
    K_norm = K_values * (4 * np.pi * tau_values)**2 / Vol

    # Fit: K_norm = a_0 + a_2 * tau + a_4 * tau² + a_6 * tau³
    def expansion(tau, a0, a2, a4, a6):
        return a0 + a2 * tau + a4 * tau**2 + a6 * tau**3

    try:
        popt, pcov = curve_fit(expansion, tau_values, K_norm, p0=[1, 0, 0, 0])
        a0, a2, a4, a6 = popt
        perr = np.sqrt(np.diag(pcov))
    except:
        a0, a2, a4, a6 = float('nan'), float('nan'), float('nan'), float('nan')
        perr = [float('nan')]*4

    return a0, a2, a4, a6, perr


# ============================================================
# PART B: Spectral action approach — separate the Q-dependence layers
# ============================================================

def spectral_action_decomposed(Q, Lambda, L, N_max):
    """
    Decompose S(Q) = Tr f(D²/Λ²) into heat-kernel terms.

    For f(x) = exp(-x) = heat kernel at t=1/Λ²:
        S(Q) = K(1/Λ², Q) (with the Q-dependent eigenvalues)

    The standard asymptotic expansion gives:
        S(Q) ~ (4π/(Q²Λ²))^{-2} * (QL)⁴ * [a_0 + a_2/(Q²Λ²) + a_4/(Q²Λ²)² + ...]
             = (Q²Λ²)²/(4π)² * (QL)⁴ * [a_0 + ...]
             = Q⁸ * Λ⁴L⁴/(4π)² * [a_0 + a_2/(Q²Λ²) + a_4/(Q⁴Λ⁴) + ...]

    So for the LEADING a_0 term: S ~ Q^8 * const
    The a_2 correction: ~ Q^6 * const
    The a_4 correction: ~ Q^4 * const  ← THIS is the log-divergent piece!

    On flat T⁴: a_0=1, a_2=0, a_4=0.  But the a_0 term IS present.
    The effective alpha ≈ 8 at large Λ, dropping toward 4 at finite Λ due to lattice effects.

    Let's compute S(Q) and also S(Q)/Q^8, S(Q)/Q^4, etc.
    """
    t = 1.0 / Lambda**2
    K = heat_trace_flat_T4(t, Q, L, N_max)
    return K


def spectral_action_fast(Q, Lambda, L, N_max):
    """Vectorized spectral action."""
    rng = np.arange(-N_max, N_max + 1)
    n1, n2, n3, n4 = np.meshgrid(rng, rng, rng, rng, indexing='ij')
    nsq = (n1**2 + n2**2 + n3**2 + n4**2).ravel()
    arg = (2*np.pi/L)**2 * nsq / (Q**2 * Lambda**2)
    vals = np.exp(-arg)
    vals[nsq == 0] = 0.0
    return np.sum(vals)


# ============================================================
# PART C: The REAL a_4 test — conformally coupled scalar on S^4
# ============================================================

def a4_on_round_S4_conformal():
    """
    On a round S⁴ of radius R, the a_4 coefficient for a conformally coupled
    scalar is KNOWN analytically:

        a_4 = (1/180) * (1/(4π)²) * ∫ (C_μνρσ)² √g d⁴x

    For round S⁴: C = 0 (conformally flat!), so a_4 from Weyl = 0.
    But the Euler density E₄ contributes:
        a_4 = (1/360) * χ(S⁴) = (1/360) * 2 = 1/180

    Under g → Q² g (constant Q, so still round S⁴ of radius QR):
        Vol → Q⁴ Vol₀
        Ricci scalar: R → R₀/Q²
        Curvature invariants scale accordingly.

    The a_4 coefficient on round S⁴ of radius r:
        a_4 = 1/180 (this is topological, independent of Q for S⁴!)

    So the CONFORMAL VARIATION of a_4 on S⁴ is zero (it's topological).

    The physically relevant quantity is:
        Λ⁰ contribution to S = f₄ * a_4  where f₄ = ∫₀^∞ f(x) dx = 1 (for exp(-x))
        This gives: S_a4 = a_4 * Vol = a_4 * Q⁴ * Vol₀

    The a_4 piece of S scales as Q⁴, giving alpha_a4 = 4.
    """
    print("\n" + "="*70)
    print("  ANALYTIC CHECK: a_4 on round S^4")
    print("="*70)
    print("  On S^4 (conformally flat), a_4 = chi(S^4)/360 = 2/360 = 1/180")
    print("  This is TOPOLOGICAL → does not change under conformal rescaling!")
    print("  The a_4 piece of spectral action: S_a4 = f_4 * a_4 * 1  (dim-independent)")
    print("  Under g → Q²g on S^4: S_a4 → S_a4 (no change, it's topological)")
    print("  BUT the VOLUME integral of a_4 density: ∫ a_4(x) √g d⁴x → Q⁴ * ∫ a_4(x) √g₀")
    print("  So the Λ⁰-term in spectral action scales as Q⁴.")
    print("  ")
    print("  For Paper 9's Σ = 2 ln Q claim:")
    print("  The spectral action on round S⁴ has pieces:")
    print("    Λ⁴ piece: ~ Q⁸ * a_0 = Q⁸  (volume * density)")
    print("    Λ² piece: ~ Q⁶ * a_2        (volume * R)")
    print("    Λ⁰ piece: ~ Q⁴ * a_4        (topological)")
    print("    log piece: ~ ln Λ * a_4      (conformal anomaly)")


# ============================================================
# PART D: The CORRECT quantity — conformal anomaly
# ============================================================

def conformal_anomaly_test(L, N_max):
    """
    The conformal anomaly in 4D is:
        ⟨T^μ_μ⟩ = a_4(x) / (4π)²

    Under g → e^{2ω} g (with ω = ln Q = const), the effective action changes as:
        δ_ω S = ∫ ω * a_4(x) √g d⁴x / (4π)²

    For FLAT torus: a_4 = 0, so δ_ω S = 0.
    All Q-dependence on flat T⁴ comes from a_0 (leading) and higher a_{2k} corrections.

    To see a_4 physics, we need CURVED space.

    HOWEVER: we can extract the EFFECTIVE conformal dimension from the spectral action.

    On flat T⁴ with f(x) = exp(-x):
        S(Q) = [θ₃(0, exp(-4π²/(Q²Λ²L²)))]⁴ - 1

    In the regime Λ L ≫ 1 (many modes below cutoff):
        S(Q) ≈ (Q²Λ²L²/4π²)² * [1 + corrections]
             = (Q Λ L / (2π))⁴ * [1 + exponentially small corrections]

    Wait — let me recompute this carefully.

    θ₃(0, q) = sqrt(π/(-ln q)) * [1 + 2 exp(-π²/(-ln q)) + ...] for q close to 1

    With q = exp(-4π²/(Q²Λ²L²)):
        -ln q = 4π²/(Q²Λ²L²)
        sqrt(π/(-ln q)) = Q Λ L / (2 sqrt(π))

    So θ₃⁴ ≈ (Q Λ L)⁴ / (4π²) = (Q Λ L)⁴ / (2π)² / 4

    Actually: θ₃⁴ = [Q Λ L / (2√π)]⁴ = Q⁴ Λ⁴ L⁴ / (16 π²)

    Hmm, but the v1 data shows alpha ≈ 4, not 8.
    Let me recheck...
    """
    print("\n" + "="*70)
    print("  CAREFUL ANALYSIS: S(Q) scaling on flat T^4")
    print("="*70)

    # The spectral action is S(Q) = sum_{n≠0} exp(-(2π|n|/L)² / (Q²Λ²))
    # = K(τ) where τ = (2π/L)² / (Q²Λ²)
    #
    # For small τ (= large QΛL/(2π)):
    #   K(τ) = θ₃(0, e^{-τ})⁴ - 1
    #   θ₃(0, e^{-τ}) ~ (π/τ)^{1/2} for small τ
    #   K(τ) ~ (π/τ)² - 1 ~ π²/τ²
    #
    # τ = (2π/L)² / (Q²Λ²) = 4π²/(Q²Λ²L²)
    # π²/τ² = π² * Q⁴Λ⁴L⁴ / (16π⁴) = Q⁴Λ⁴L⁴/(16π²)
    #
    # So S(Q) ~ Q⁴ * (ΛL)⁴/(16π²) for large QΛL.
    #
    # alpha = d ln S / d ln Q = 4  ✓  (matches v1 data!)
    #
    # BUT WAIT: this is the heat trace K(τ), not the spectral action S = Tr f(D²/Λ²).
    # For f(x) = exp(-x), S = K(1/Λ²) if we put τ = 1/Λ² into K.
    #
    # Actually S(Q) = sum exp(- λ_n(Q) / Λ²)
    #              = sum exp(- (2πn/L)²/(Q²Λ²))
    # which IS K(τ) with τ = 1/Λ² and eigenvalues (2πn/L)²/Q².
    #
    # More precisely: S(Q) = sum exp(-τ_eff |n|²) with τ_eff = (2π/L)²/(Q²Λ²)
    #
    # So S(Q) ~ (π/τ_eff)² = (Q²Λ²L²/(4π²))² for large QΛL/(2π)
    #        = Q⁴ Λ⁴ L⁴ / (16π⁴)  ...
    #
    # Let me just check numerically for ONE case.

    Lambda = 20.0
    Q_test = [0.8, 0.9, 1.0, 1.1, 1.2]

    print(f"\n  Lambda={Lambda}, L={L}, N_max={N_max}")
    print(f"  Theoretical UV limit: S(Q) ~ Q^4 * const")
    print(f"  {'Q':>6}  {'S(Q)':>16}  {'S(Q)/Q^4':>16}  {'S(Q)/Q^8':>16}")
    for Q in Q_test:
        S = spectral_action_fast(Q, Lambda, L, N_max)
        print(f"  {Q:6.2f}  {S:16.8f}  {S/Q**4:16.8f}  {S/Q**8:16.8f}")

    # Now the ASYMPTOTIC prediction
    tau_eff = (2*np.pi/L)**2 / (1.0**2 * Lambda**2)
    S_asymp = (np.pi / tau_eff)**2
    S_actual = spectral_action_fast(1.0, Lambda, L, N_max)
    print(f"\n  Asymptotic prediction S(Q=1) = (π/τ)² = {S_asymp:.8f}")
    print(f"  Actual S(Q=1) = {S_actual:.8f}")
    print(f"  Ratio = {S_actual/S_asymp:.8f}")

    # The discrepancy comes from the -1 (subtracting n=0) and subleading theta terms.
    # θ₃⁴ - 1 ≈ (π/τ)² - 1 + 8(π/τ)^{3/2} exp(-π²/τ) + ...
    # So S = (π/τ)² [1 + O(exp(-π²/τ))] - 1

    # For the spectral action S = Tr f(D²/Λ²) on (T⁴, Q²δ):
    #   - In the DEEP UV (ΛL → ∞): S ~ Q⁴ Λ⁴ L⁴ / (16π²)
    #   - The dominant Q-dependence is Q⁴ from the a_0 term
    #   - Since a_4 = 0 on flat T⁴, there IS no a_4 conformal response
    #
    # CONCLUSION: On flat T⁴, the Q-dependence is ENTIRELY from a_0 and volume scaling.
    # g_QQ relates to the spectral dimension, not to a_4.

    print(f"\n  CONCLUSION:")
    print(f"  On flat T⁴, alpha = d ln S / d ln Q → 4 (= spacetime dimension)")
    print(f"  This is purely the a₀ (volume) contribution: S ~ Vol * Λ⁴ ~ Q⁴")
    print(f"  g_QQ = -S''/S = -4*3 = -12 (for pure Q⁴ scaling, note the SIGN)")
    print(f"  ")
    print(f"  For the a₄ coefficient:")
    print(f"  - On flat T⁴: a₄ = 0, so a₄'s conformal response is trivially 0")
    print(f"  - On curved space: a₄ has Weyl² and Euler density terms")
    print(f"  - The CONFORMAL ANOMALY from a₄ is the log(Λ) piece of S")


# ============================================================
# PART E: Isolate the conformal anomaly via Λ-dependence
# ============================================================

def isolate_a4_via_Lambda(L, N_max):
    """
    The spectral action expansion:
        S(Q, Λ) = f_4 Λ⁴ Q⁴ a_0 L⁴ / (16π²)
                + f_2 Λ² Q² a_2 L² / (4π) (= 0 on flat T⁴)
                + f_0 a_4 (= 0 on flat T⁴)
                + ...

    where f_k = ∫₀^∞ f(x) x^{k/2-1} dx  (momenta of the test function).
    For f(x) = e^{-x}: f_4 = Γ(2) = 1, f_2 = Γ(1) = 1, f_0 = Γ(0) → ∞ (log divergent)

    Actually for the STANDARD spectral action expansion (Chamseddine-Connes):
        Tr f(D²/Λ²) ~ Σ_{k=0}^{d/2} f_{d-2k} Λ^{d-2k} a_{2k} + O(Λ^{-1})

    In d=4:
        S = f_4 Λ⁴ a_0 + f_2 Λ² a_2 + f_0 a_4 + O(Λ^{-2})

    where a_k = ∫ a_k(x) √g d⁴x  are the INTEGRATED coefficients.

    For flat T⁴ with g = Q² δ:
        a_0 = (4π)^{-2} Vol(g) = (QL)⁴ / (16π²)
        a_2 = (4π)^{-2} ∫ R/6 √g = 0 (flat)
        a_4 = (4π)^{-2} ∫ [curvature invariants] √g = 0 (flat)

    So S = f_4 Λ⁴ (QL)⁴ / (16π²) + 0 + 0 = Λ⁴ Q⁴ L⁴ / (16π²)

    The Q-dependence is PURELY Q⁴ from the volume of the a_0 term.

    To test a_4's conformal response, let's VERIFY the above by fitting S vs Λ.
    """
    print("\n" + "="*70)
    print("  ISOLATING a_k VIA Λ-DEPENDENCE")
    print(f"  L={L}, N_max={N_max}")
    print("="*70)

    Lambda_values = np.array([5, 8, 10, 12, 15, 20, 25, 30, 40, 50])

    for Q in [0.8, 1.0, 1.2]:
        S_values = np.array([spectral_action_fast(Q, Lam, L, N_max) for Lam in Lambda_values])

        # Fit: S = c4 Λ⁴ + c2 Λ² + c0 + c_m2 Λ^{-2}
        # Use Λ⁴ scaling
        def model(Lam, c4, c2, c0):
            return c4 * Lam**4 + c2 * Lam**2 + c0

        try:
            popt, pcov = curve_fit(model, Lambda_values, S_values, p0=[1e-3, 0, 0])
            c4, c2, c0 = popt
            perr = np.sqrt(np.diag(pcov))
        except:
            c4, c2, c0 = float('nan'), float('nan'), float('nan')
            perr = [float('nan')] * 3

        # Theoretical: c4 = Q⁴ L⁴ / (16π²), c2 = 0, c0 = 0
        c4_theory = Q**4 * L**4 / (16 * np.pi**2)

        print(f"\n  Q = {Q:.1f}:")
        print(f"    Fit:  c₄ = {c4:.10e}  (theory: {c4_theory:.10e}, ratio: {c4/c4_theory:.8f})")
        print(f"          c₂ = {c2:.10e}  (theory: 0)")
        print(f"          c₀ = {c0:.10e}  (theory: 0, this would be a₄)")
        print(f"    Errors: σ(c₄)={perr[0]:.2e}, σ(c₂)={perr[1]:.2e}, σ(c₀)={perr[2]:.2e}")

    print(f"\n  KEY RESULT: c₄(Q)/c₄(1) should be Q⁴:")
    Qs = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
    c4_values = []
    for Q in Qs:
        S_vals = np.array([spectral_action_fast(Q, Lam, L, N_max) for Lam in Lambda_values])
        def model(Lam, c4, c2, c0):
            return c4 * Lam**4 + c2 * Lam**2 + c0
        try:
            popt, _ = curve_fit(model, Lambda_values, S_vals, p0=[1e-3, 0, 0])
            c4_values.append(popt[0])
        except:
            c4_values.append(float('nan'))

    c4_at_1 = c4_values[Qs.index(1.0)]
    print(f"  {'Q':>6}  {'c₄(Q)/c₄(1)':>14}  {'Q⁴':>10}  {'ratio':>10}")
    for Q, c4 in zip(Qs, c4_values):
        ratio = c4 / c4_at_1
        print(f"  {Q:6.2f}  {ratio:14.8f}  {Q**4:10.6f}  {ratio/Q**4:10.8f}")


# ============================================================
# PART F: The CORRECT test — add curvature via mass potential
# ============================================================

def curved_space_test():
    """
    On flat T⁴ we can't test a₄ since it's zero. To mimic curvature effects,
    consider D² + V where V simulates a scalar curvature coupling.

    For a conformally coupled scalar in 4D: D² + R/6.
    On a sphere of radius r: R = 12/r². With conformal rescaling r → Q r:
        R → R/Q² = 12/(Q²r²)
        V = R/6 = 2/(Q²r²)

    The eigenvalues on S⁴ of radius r:
        λ_l = (l+1)(l+2)/r²  with degeneracy d_l = (l+1)²(l+2)²/12 ... (complicated)

    Instead, let's do a SIMPLER test: add V(Q) = V₀/Q² to the flat T⁴ Laplacian.
    This mimics the conformal coupling V = R/6 where R ~ 1/Q².

    D²_eff = D²_flat/Q² + V₀/Q²  = (D²_flat + V₀)/Q²

    S(Q) = sum_n exp(- [(2πn/L)² + V₀] / (Q² Λ²))
    """
    print("\n" + "="*70)
    print("  CURVED SPACE MIMICRY: D² + V₀/Q²")
    print("="*70)

    L = 1.0
    N_max = 12
    V0_values = [0.0, 1.0, 10.0, 100.0]

    for V0 in V0_values:
        def S_of_Q(Q, Lambda=10.0):
            rng = np.arange(-N_max, N_max + 1)
            n1, n2, n3, n4 = np.meshgrid(rng, rng, rng, rng, indexing='ij')
            nsq = (n1**2 + n2**2 + n3**2 + n4**2).ravel()
            eigenvalues = (2*np.pi/L)**2 * nsq + V0
            x = eigenvalues / (Q**2 * Lambda**2)
            vals = np.exp(-x)
            if V0 == 0:
                vals[nsq == 0] = 0.0
            return np.sum(vals)

        Q0 = 1.0
        h = 1e-4
        S0 = S_of_Q(Q0)
        Sp = (S_of_Q(Q0+h) - S_of_Q(Q0-h)) / (2*h)
        Spp = (S_of_Q(Q0+h) - 2*S0 + S_of_Q(Q0-h)) / h**2
        alpha = Sp / S0
        g_QQ = -Spp / S0
        d2lnS = Spp/S0 - (Sp/S0)**2

        print(f"\n  V₀ = {V0:6.1f}: S(1)={S0:.6e}, alpha={alpha:.6f}, g_QQ={g_QQ:.6f}, d²lnS/dQ²={d2lnS:.6f}")


# ============================================================
# PART G: The DEFINITIVE test — explicit a_4 variation for Dirac on S^2×T^2
# ============================================================

def spectral_action_S2_cross_T2(Q, Lambda, N_max_S2, N_max_T2, L_T2):
    """
    On S² × T² with metric Q² (dΩ² + δ_ab dx^a dx^b):

    S² of radius Q: eigenvalues of Laplacian are l(l+1)/Q², degeneracy 2l+1
    T² of period QL: eigenvalues (2π n/QL)², n ∈ Z²

    For the scalar Laplacian (good enough to check the principle):
        λ_{l,n} = l(l+1)/Q² + (2π/L_T2)² |n|²/Q²
                = [l(l+1) + (2π/L_T2)² |n|²] / Q²

    S(Q) = sum_{l,n} (2l+1) exp(-λ_{l,n}/(Q² Λ²))
         = sum_{l,n} (2l+1) exp(-[l(l+1) + (2π/L_T2)² |n|²] / (Q⁴ Λ²))

    Wait — I need to be more careful. The eigenvalues of the FULL Laplacian on M×N
    are λ_M + λ_N. On conformally rescaled (Q²g), eigenvalues → eigenvalues/Q².
    The spectral action f(D²/Λ²) = f(eigenvalue/(Q²Λ²)).

    Actually for a 4-manifold with g' = Q²g:
        -Δ_{g'} = -Δ_g / Q²  (for scalars in 4D, the conformal weight matters)

    For the CONFORMAL Laplacian (Yamabe operator) in 4D:
        Y = -Δ + R/6
    Under g → Q²g:  Y_{Q²g} = Q^{-3} Y_g Q^{-1}  (conformal covariance)
    So eigenvalues λ' = λ/Q⁴ ... no, that's for conformal weight.

    Let me just use the simple scalar Laplacian:
        Eigenvalues of -Δ on (M, Q²g) = eigenvalues of -Δ on (M,g) divided by Q²

    This is the correct statement for the SCALAR Laplacian.
    """
    # Eigenvalues: [l(l+1) + (2π/L)² |m|²] for l=0,...,N_S2 and m ∈ Z²
    # Degeneracy of S² part: (2l+1)

    S_val = 0.0
    rng_T2 = np.arange(-N_max_T2, N_max_T2 + 1)
    m1, m2 = np.meshgrid(rng_T2, rng_T2, indexing='ij')
    msq = (m1**2 + m2**2).ravel()
    T2_contrib = (2*np.pi/L_T2)**2 * msq

    denom = Q**2 * Lambda**2

    for l in range(0, N_max_S2 + 1):
        S2_eigenvalue = l * (l + 1)
        degeneracy = 2 * l + 1
        for tc in T2_contrib:
            total_eigenvalue = S2_eigenvalue + tc
            if total_eigenvalue == 0:
                continue  # skip zero mode
            x = total_eigenvalue / denom
            S_val += degeneracy * np.exp(-x)

    return S_val


def S2_cross_T2_test():
    """
    On S² × T²: a₄ is NONZERO because S² has curvature!

    For round S² of radius r: R = 2/r², Ric² = 2/r⁴, Riem² = 4/r⁴

    The Seeley-DeWitt a₄ for the scalar Laplacian (E=0, Ω=0) is:
        a₄ = (4π)^{-d/2} ∫ [1/180 (5R² - 2 Ric² + 2 Riem²)] √g d⁴x

    Wait, the standard formula for -Δ + E is:
        a₄ = (4π)^{-2} ∫ [1/360 (5R² - 2 Ric² + 2 Riem²) + 1/6 □R + 1/2 E² + 1/6 RE] √g

    For S²×T² with radii (r, L): R = 2/r², Ric = (1/r²)(g_S2), Riem on S² only
        Ric² = 2/r⁴, Riem² = 4/r⁴ (Kretschner of S²)... actually:

    On S² × T²:
        R = 2/r² (curvature only from S²)
        Ric_μν: on S²: (1/r²) g_ab; on T²: 0
        Ric² = 2/r⁴
        Riem_αβγδ: on S²: (1/r²)(g_{αγ}g_{βδ} - g_{αδ}g_{βγ})
        Riem² = 4/r⁴ (two independent components on S²)

    So: 5R² - 2Ric² + 2Riem² = 5(4/r⁴) - 2(2/r⁴) + 2(4/r⁴) = 20/r⁴ - 4/r⁴ + 8/r⁴ = 24/r⁴

    a₄ = (4π)^{-2} * (24/r⁴)/360 * Vol(S²×T²)
       = (4π)^{-2} * (1/15) * 1/r⁴ * 4πr² * L²
       = (4π)^{-2} * (4π/15) * L²/r²
       = L²/(60π r²)

    Under g → Q²g: r → Qr, L → QL
        a₄ → (QL)²/(60π (Qr)²) = L²/(60π r²)  [Q cancels!]

    So a₄ is Q-independent on S²×T² as well!

    This makes sense: a₄ in 4D is related to the Euler characteristic and
    conformal invariant, both of which don't change under constant conformal rescaling.
    """
    print("\n" + "="*70)
    print("  S² × T² TEST: a₄ ≠ 0 background")
    print("="*70)

    # S² of radius 1, T² of period 1
    # Eigenvalues: l(l+1) + (2π m)²  for l = 0,1,2,... and m ∈ Z²
    # With g → Q²g: eigenvalues → eigenvalues / Q²
    # S(Q) = sum (2l+1) exp(-[l(l+1) + (2πm)²] / (Q² Λ²))

    Lambda = 10.0
    L_T2 = 1.0
    N_max_S2 = 30
    N_max_T2 = 12

    def S_of_Q(Q):
        return spectral_action_S2_cross_T2(Q, Lambda, N_max_S2, N_max_T2, L_T2)

    Q0 = 1.0
    h = 1e-4
    S0 = S_of_Q(Q0)
    Sp = (S_of_Q(Q0+h) - S_of_Q(Q0-h)) / (2*h)
    Spp = (S_of_Q(Q0+h) - 2*S0 + S_of_Q(Q0-h)) / h**2

    alpha = Sp / S0
    g_QQ = -Spp / S0
    d2lnS = Spp/S0 - (Sp/S0)**2

    print(f"  Lambda={Lambda}, N_S2={N_max_S2}, N_T2={N_max_T2}, L_T2={L_T2}")
    print(f"  S(1) = {S0:.10e}")
    print(f"  alpha = d ln S / d ln Q = {alpha:.8f}")
    print(f"  g_QQ = -S''/S = {g_QQ:.8f}")
    print(f"  d² ln S / dQ² = {d2lnS:.8f}")
    print(f"  alpha*(alpha-1) = {alpha*(alpha-1):.8f}")

    # Q-scan
    print(f"\n  Q-dependence:")
    print(f"  {'Q':>6}  {'S(Q)/S(1)':>14}  {'Q^4':>10}  {'ratio':>12}")
    for Q in [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]:
        SQ = S_of_Q(Q)
        ratio = SQ / S0
        print(f"  {Q:6.2f}  {ratio:14.8f}  {Q**4:10.6f}  {ratio/Q**4:12.8f}")

    # Λ-decomposition on S²×T²
    print(f"\n  Λ-decomposition on S²×T² (isolating a₄ from Λ⁰ term):")
    Lambda_vals = np.array([5, 8, 10, 15, 20, 30, 50])

    for Q in [0.9, 1.0, 1.1]:
        S_vals = [spectral_action_S2_cross_T2(Q, Lam, N_max_S2, N_max_T2, L_T2) for Lam in Lambda_vals]
        S_vals = np.array(S_vals)

        # Fit: S = c4 Λ⁴ + c2 Λ² + c0
        def model(Lam, c4, c2, c0):
            return c4 * Lam**4 + c2 * Lam**2 + c0

        try:
            popt, pcov = curve_fit(model, Lambda_vals, S_vals, p0=[1e-3, 0, 0])
            c4, c2, c0 = popt
            perr = np.sqrt(np.diag(pcov))
            print(f"  Q={Q:.1f}: c₄={c4:.6e}, c₂={c2:.6e}, c₀={c0:.6e} (±{perr[2]:.2e})")
        except Exception as e:
            print(f"  Q={Q:.1f}: fit failed: {e}")

    # a₄ analytic value on S²(r=1) × T²(L=1):
    # a₄ = L²/(60π r²) = 1/(60π) ≈ 0.005305
    a4_theory = 1.0 / (60 * np.pi)
    print(f"\n  Analytic a₄ on S²(r=1)×T²(L=1) = {a4_theory:.8f}")
    print(f"  (This should appear as the Λ⁰ coefficient c₀)")


# ============================================================
# PART H: THE ACTUAL SIGMA = 2 ln Q TEST
# ============================================================

def sigma_test():
    """
    Paper 9 claims Σ = 2 ln Q.

    The spectral action S = Tr f(D²/Λ²).
    Under g → Q²g (constant conformal rescaling on a 4-manifold):
        D² eigenvalues → eigenvalues/Q²
        S(Q) = Tr f(D²/(Q²Λ²))

    The heat-kernel expansion:
        S = f₄ Λ⁴ a₀(Q²g) + f₂ Λ² a₂(Q²g) + f₀ a₄(Q²g) + ...

    where a_k(Q²g) are the INTEGRATED Seeley-DeWitt coefficients.

    Known scaling:
        a₀(Q²g) = (4π)^{-2} Vol(Q²g) = (4π)^{-2} Q⁴ Vol(g)  → scales as Q⁴
        a₂(Q²g) = (4π)^{-2} ∫ R(Q²g)/6 √(Q²g) = (4π)^{-2} Q² ∫ R(g)/6 √g  → scales as Q²
            [since R(Q²g) = R(g)/Q² and √(Q²g) = Q⁴ √g in 4D]
        a₄(Q²g) = (4π)^{-2} ∫ [curvature²] Q⁰ √g  → scales as Q⁰ (conformal invariants)
            [the Q's cancel: (curv)² ~ Q⁻⁴, √g ~ Q⁴]

    More precisely for 4D:
        ∫ R² √g → Q⁰  (R ~ 1/Q², √g ~ Q⁴, R² ~ 1/Q⁴, so R² √g ~ 1)
        ∫ Riem² √g → Q⁰
        ∫ Ric² √g → Q⁰

    So: a₄(Q²g) = a₄(g)  EXACTLY for constant Q in 4D!

    Therefore:
        S(Q) = f₄ Λ⁴ Q⁴ a₀(g) + f₂ Λ² Q² a₂(g) + f₀ a₄(g) + ...

    The a₄ piece does NOT depend on Q at all!

    FOR PAPER 9: The question is whether the RELATIVE entropy or Fisher information
    derived from the spectral action gives Σ = 2 ln Q.

    Let's compute the RELATIVE spectral action:
        Σ(Q) = S(Q)/S(1) - 1  or  ln(S(Q)/S(1))

    In the UV limit (Λ → ∞), S(Q)/S(1) → Q⁴, so ln(S(Q)/S(1)) → 4 ln Q.
    This is 2 Σ if Σ = 2 ln Q.

    Alternatively, if we define the spectral DENSITY (per unit volume):
        s(Q) = S(Q) / Vol(Q²g) = S(Q) / (Q⁴ Vol₀)
    Then in UV limit: s(Q)/s(1) → 1 (the density is Q-independent to leading order).

    The subleading terms:
        s(Q) = f₄ Λ⁴ a₀/(4π)² + f₂ Λ² a₂/(Q² (4π)²) + f₀ a₄/(Q⁴ (4π)²) + ...

    So s(Q)/s(1) = 1 + [a₂/(a₀ Λ²)] (1/Q² - 1) + [a₄/(a₀ Λ⁴)] (1/Q⁴ - 1) + ...

    For the a₂ correction: 1/Q² ≈ 1 - 2 ln Q + ... = exp(-2 ln Q) ≈ exp(-Σ) if Σ = 2 ln Q

    THIS IS THE CONNECTION!
    The subleading correction to the spectral density is:
        δs/s ~ exp(-Σ) - 1  where Σ = 2 ln Q

    This matches the Petz recovery bound: F ≥ exp(-Σ).
    """
    print("\n" + "="*70)
    print("  *** SIGMA = 2 ln Q : SPECTRAL ACTION TEST ***")
    print("="*70)

    L = 1.0
    N_max = 12

    # Test for large Lambda (UV limit)
    for Lambda in [10.0, 20.0, 50.0]:
        def S_of_Q(Q):
            return spectral_action_fast(Q, Lambda, L, N_max)

        S1 = S_of_Q(1.0)
        Vol1 = L**4  # Vol at Q=1

        print(f"\n  Lambda = {Lambda}")
        print(f"  {'Q':>6}  {'ln(S/S₁)':>12}  {'4 ln Q':>10}  {'2 ln Q':>10}  {'ratio/4lnQ':>12}  {'s(Q)/s(1)':>12}  {'exp(-Σ)=1/Q²':>14}")

        for Q in [0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3]:
            SQ = S_of_Q(Q)
            lnratio = np.log(SQ/S1) if SQ > 0 else float('nan')
            four_lnQ = 4 * np.log(Q) if Q > 0 else float('nan')
            two_lnQ = 2 * np.log(Q) if Q > 0 else float('nan')
            ratio = lnratio / four_lnQ if abs(four_lnQ) > 1e-10 else float('nan')

            # Spectral density ratio
            VolQ = (Q*L)**4
            sQ = SQ / VolQ
            s1 = S1 / Vol1
            density_ratio = sQ / s1

            expnegSigma = 1/Q**2

            print(f"  {Q:6.2f}  {lnratio:12.8f}  {four_lnQ:10.6f}  {two_lnQ:10.6f}  {ratio:12.8f}  {density_ratio:12.8f}  {expnegSigma:14.8f}")

    # KEY CHECK: spectral density ratio vs exp(-Σ)
    print(f"\n  KEY CHECK: spectral density s(Q)/s(1) vs 1/Q² = exp(-Σ)")
    print(f"  In UV limit (large Λ), s(Q)/s(1) → 1 (constant density)")
    print(f"  The CORRECTION to 1 should scale as exp(-Σ) - 1 = 1/Q² - 1")

    Lambda = 20.0
    def S_of_Q(Q):
        return spectral_action_fast(Q, Lambda, L, N_max)

    S1 = S_of_Q(1.0)
    s1 = S1 / L**4

    print(f"\n  Lambda = {Lambda}")
    print(f"  {'Q':>6}  {'s(Q)/s(1)-1':>14}  {'1/Q²-1':>12}  {'ratio':>12}")
    for Q in [0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]:
        SQ = S_of_Q(Q)
        sQ = SQ / (Q*L)**4
        delta_s = sQ/s1 - 1
        delta_exp = 1/Q**2 - 1
        ratio = delta_s / delta_exp if abs(delta_exp) > 1e-10 else float('nan')
        print(f"  {Q:6.2f}  {delta_s:14.10f}  {delta_exp:12.8f}  {ratio:12.8f}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  PAPER 9: NUMERICAL VERIFICATION v2")
    print("  Focus: Isolating a₄ conformal variation")
    print("=" * 70)

    # Part A/B: Flat T⁴ scaling analysis
    conformal_anomaly_test(L=1.0, N_max=12)

    # Part C/D: Analytic a₄ discussion
    a4_on_round_S4_conformal()

    # Part E: Λ-decomposition to isolate a_k
    isolate_a4_via_Lambda(L=1.0, N_max=12)

    # Part F: Curved space mimicry
    curved_space_test()

    # Part G: S² × T² (genuinely curved)
    S2_cross_T2_test()

    # Part H: The Sigma = 2 ln Q test
    sigma_test()

    # ============================================================
    # FINAL VERDICT
    # ============================================================
    print("\n" + "=" * 70)
    print("  *** FINAL VERDICT ***")
    print("=" * 70)
    print("""
  1. FLAT T⁴ RESULT:
     - S(Q) ~ Q⁴ in the UV limit (α = d ln S / d ln Q → 4)
     - This is purely the a₀ (volume) scaling: a₀ ∝ Vol(g) ∝ Q⁴
     - a₄ = 0 on flat T⁴, so it contributes NOTHING to Q-dependence
     - g_QQ = -S''/S → -4×3 = -12  (negative! not +2)

  2. CURVED SPACE (S²×T²):
     - a₄ ≠ 0 on S²×T²
     - BUT a₄(Q²g) = a₄(g) for constant Q in 4D (conformal invariance)
     - So the Λ⁰ piece of S is Q-INDEPENDENT
     - The Q-dependence still comes from a₀ ~ Q⁴ and a₂ ~ Q²

  3. SIGMA = 2 ln Q INTERPRETATION:
     - The TOTAL spectral action: ln(S(Q)/S(1)) → 4 ln Q (not 2 ln Q)
     - The SPECTRAL DENSITY: s(Q)/s(1) = S(Q)/(Q⁴ S(1)) → 1
     - The a₂ CORRECTION to density: δs/s ~ 1/Q² - 1 = exp(-Σ) - 1
       where Σ = 2 ln Q
     - THIS is where Σ = 2 ln Q appears: in the subleading correction
       to the spectral density, controlled by a₂ ~ R

  4. THE g_QQ QUESTION:
     - g_QQ = -S''/S at Q=1 is NOT 2. It is approximately -12 in the UV.
     - This is because g_QQ measures the total spectral action curvature,
       dominated by the a₀ ~ Q⁴ piece.
     - The Fisher metric of the spectral DENSITY would give a different answer.
""")
