#!/usr/bin/env python3
"""
Paper 9 — Independent Numerical Verification
=============================================
Conformal variation of the Seeley-DeWitt a_4 coefficient
on a flat 4-torus T^4 with metric g_uv = Q^2 delta_uv.

We compute the spectral action  S(Q) = sum_n f( (2pi|n|/L)^2 / (Q^2 Lambda^2) )
for f(x) = exp(-x)  [heat-kernel regulator] and extract:

    g_QQ  =  -S''(Q) / S(Q)   evaluated at Q = 1.

If Sigma = 2 ln Q is the correct spectral-action response to a conformal
rescaling, the Fisher metric g_QQ should equal a definite number that we
determine here WITHOUT assuming the answer.

Author : numerical verification script (independent of analytic calculation)
Date   : 2026-03-19
"""

import numpy as np
from itertools import product
import sys

# ============================================================
# 1.  Core computation: spectral action on T^4
# ============================================================

def spectral_action_T4(Q, Lambda, L, N_max, f_func, m2=0.0):
    """
    Compute S(Q) = sum_{n in Z^4, 0<|n|<=N_max} f( [(2pi|n|/L)^2 + m^2] / (Q^2 Lambda^2) )

    Parameters
    ----------
    Q       : conformal factor (metric = Q^2 delta)
    Lambda  : energy cutoff
    L       : torus period
    N_max   : lattice cutoff on each component of n
    f_func  : cutoff function f(x)
    m2      : mass^2 term added to D^2

    Returns
    -------
    S : the spectral action value (float)
    """
    prefactor = (2 * np.pi / L) ** 2
    denom = Q**2 * Lambda**2
    S = 0.0
    # Sum over n in Z^4 with each component in [-N_max, N_max]
    # Exclude n=0 only if m2=0 (massless); if m2>0 n=0 contributes f(m2/denom)
    rng = range(-N_max, N_max + 1)
    for n1 in rng:
        for n2 in rng:
            for n3 in rng:
                for n4 in rng:
                    nsq = n1**2 + n2**2 + n3**2 + n4**2
                    if nsq == 0 and m2 == 0.0:
                        continue  # skip zero mode for massless
                    eigenvalue = prefactor * nsq + m2
                    x = eigenvalue / denom
                    S += f_func(x)
    return S


def spectral_action_T4_fast(Q, Lambda, L, N_max, f_func, m2=0.0):
    """
    Vectorized version for speed.
    """
    rng = np.arange(-N_max, N_max + 1)
    # Build all 4-tuples
    n1, n2, n3, n4 = np.meshgrid(rng, rng, rng, rng, indexing='ij')
    nsq = n1**2 + n2**2 + n3**2 + n4**2
    nsq_flat = nsq.ravel()

    prefactor = (2 * np.pi / L) ** 2
    denom = Q**2 * Lambda**2

    eigenvalues = prefactor * nsq_flat + m2
    x = eigenvalues / denom

    vals = f_func(x)

    if m2 == 0.0:
        # Zero out the n=0 contribution
        mask = (nsq_flat == 0)
        vals[mask] = 0.0

    return np.sum(vals)


# ============================================================
# 2.  Numerical derivatives
# ============================================================

def numerical_first_derivative(func, x0, h=1e-5):
    """Central difference for f'(x0)."""
    return (func(x0 + h) - func(x0 - h)) / (2 * h)


def numerical_second_derivative(func, x0, h=1e-5):
    """Central difference for f''(x0)."""
    return (func(x0 + h) - 2 * func(x0) + func(x0 - h)) / h**2


# ============================================================
# 3.  Main verification
# ============================================================

def run_verification(Lambda, L, N_max, m2=0.0, label=""):
    """Run full verification for given parameters."""

    f_func = lambda x: np.exp(-x)

    def S_of_Q(Q):
        return spectral_action_T4_fast(Q, Lambda, L, N_max, f_func, m2=m2)

    Q0 = 1.0
    S0 = S_of_Q(Q0)
    S_prime = numerical_first_derivative(S_of_Q, Q0, h=1e-4)
    S_double_prime = numerical_second_derivative(S_of_Q, Q0, h=1e-4)

    # g_QQ defined as -S''(Q)/S(Q) at Q=1
    g_QQ_fisher = -S_double_prime / S0

    # Also compute d^2(ln S)/dQ^2 at Q=1
    def lnS_of_Q(Q):
        return np.log(S_of_Q(Q))

    d2_lnS = numerical_second_derivative(lnS_of_Q, Q0, h=1e-4)

    # Analytic prediction from heat-kernel asymptotics:
    # S(Q) ~ Q^d * Lambda^d * (...)  where d=4
    # For f(x)=exp(-x) on T^4: leading term goes as Q^4 * Lambda^4 * L^4 / (2pi)^4 * (...)
    # S'/S at Q=1 should be related to the dimension d
    S_prime_over_S = S_prime / S0

    print(f"\n{'='*60}")
    print(f"  VERIFICATION: {label}")
    print(f"{'='*60}")
    print(f"  Parameters: Lambda={Lambda}, L={L}, N_max={N_max}, m^2={m2}")
    print(f"  S(Q=1)           = {S0:.10e}")
    print(f"  S'(Q=1)          = {S_prime:.10e}")
    print(f"  S''(Q=1)         = {S_double_prime:.10e}")
    print(f"  S'(1)/S(1)       = {S_prime_over_S:.10f}")
    print(f"  -S''(1)/S(1)     = {g_QQ_fisher:.10f}  <-- g_QQ (Fisher)")
    print(f"  d^2(ln S)/dQ^2   = {d2_lnS:.10f}")
    print(f"{'='*60}")

    return {
        'S0': S0, 'S_prime': S_prime, 'S_double_prime': S_double_prime,
        'g_QQ_fisher': g_QQ_fisher, 'd2_lnS': d2_lnS,
        'S_prime_over_S': S_prime_over_S,
        'Lambda': Lambda, 'L': L, 'N_max': N_max, 'm2': m2
    }


# ============================================================
# 4.  Q-dependence scan
# ============================================================

def Q_scan(Lambda, L, N_max, m2=0.0):
    """Compute S(Q)/S(1) for a range of Q values."""
    f_func = lambda x: np.exp(-x)

    Q_values = np.linspace(0.5, 2.0, 31)
    S1 = spectral_action_T4_fast(1.0, Lambda, L, N_max, f_func, m2=m2)

    print(f"\n{'='*60}")
    print(f"  Q-DEPENDENCE SCAN  (Lambda={Lambda}, L={L}, N_max={N_max})")
    print(f"{'='*60}")
    print(f"  {'Q':>6s}  {'S(Q)/S(1)':>14s}  {'Q^4':>10s}  {'Q^8':>12s}  {'ratio/Q^4':>12s}  {'ln(S/S1)/lnQ':>14s}")
    print(f"  {'-'*6}  {'-'*14}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*14}")

    ratios = []
    for Q in Q_values:
        SQ = spectral_action_T4_fast(Q, Lambda, L, N_max, f_func, m2=m2)
        ratio = SQ / S1
        q4 = Q**4
        q8 = Q**8
        ratio_over_q4 = ratio / q4 if q4 > 0 else float('nan')
        if Q != 1.0 and Q > 0:
            ln_ratio_over_lnQ = np.log(ratio) / np.log(Q)
        else:
            ln_ratio_over_lnQ = float('nan')
        ratios.append((Q, ratio, q4, q8, ratio_over_q4, ln_ratio_over_lnQ))
        print(f"  {Q:6.3f}  {ratio:14.8f}  {q4:10.6f}  {q8:12.6f}  {ratio_over_q4:12.8f}  {ln_ratio_over_lnQ:14.6f}")

    return ratios


# ============================================================
# 5.  Analytic cross-check: exact heat-kernel on T^4
# ============================================================

def analytic_heat_kernel_T4(Q, Lambda, L, N_max):
    """
    For f(x) = exp(-x), the spectral action is:
        S(Q) = sum_{n != 0} exp( -(2pi)^2 |n|^2 / (L^2 Q^2 Lambda^2) )

    Let t = (2pi)^2 / (L^2 Q^2 Lambda^2).  Then S = sum exp(-t |n|^2) - 1.

    Using Jacobi theta: sum_{n in Z} exp(-t n^2) = theta_3(0, exp(-t))
    For 4D:  sum_{n in Z^4} exp(-t |n|^2) = [theta_3(0, exp(-t))]^4

    So S(Q) = [theta_3(0, q)]^4 - 1  where q = exp(-t).

    For small t (large Q*Lambda): theta_3 ~ sqrt(pi/t) * [1 + 2 exp(-pi^2/t) + ...]
    So S ~ (pi/t)^2 - 1 ~ (L Q Lambda / 2)^4 / (something)...

    Let's compute this directly.
    """
    t = (2 * np.pi)**2 / (L**2 * Q**2 * Lambda**2)
    q_param = np.exp(-t)

    # Jacobi theta_3(0, q) = 1 + 2*sum_{n=1}^{inf} q^{n^2}
    theta3 = 1.0
    for n in range(1, N_max + 1):
        theta3 += 2 * q_param**(n**2)

    S_analytic = theta3**4 - 1.0  # subtract the n=0 term
    return S_analytic


# ============================================================
# 6.  Power-law extraction
# ============================================================

def extract_power_law(Lambda, L, N_max):
    """
    Fit S(Q) ~ A * Q^alpha near Q=1 to extract the exponent alpha.
    """
    f_func = lambda x: np.exp(-x)
    Q_vals = np.array([0.9, 0.95, 1.0, 1.05, 1.1])
    S_vals = np.array([spectral_action_T4_fast(Q, Lambda, L, N_max, f_func) for Q in Q_vals])

    # ln S = ln A + alpha * ln Q
    lnQ = np.log(Q_vals)
    lnS = np.log(S_vals)

    # Linear fit
    coeffs = np.polyfit(lnQ, lnS, 2)  # quadratic in ln Q

    print(f"\n{'='*60}")
    print(f"  POWER-LAW EXTRACTION  (Lambda={Lambda}, L={L}, N_max={N_max})")
    print(f"{'='*60}")
    print(f"  Fit: ln S = {coeffs[0]:.6f} (ln Q)^2 + {coeffs[1]:.6f} (ln Q) + {coeffs[2]:.6f}")
    print(f"  Linear power:  alpha = {coeffs[1]:.6f}")
    print(f"  Curvature:     beta  = {coeffs[0]:.6f}")
    print(f"  Note: S(Q) ~ Q^alpha * exp(beta (ln Q)^2)")
    print(f"{'='*60}")

    return coeffs


# ============================================================
# 7.  Heat-kernel coefficient extraction
# ============================================================

def extract_a_coefficients(L, N_max_lattice):
    """
    The heat-kernel expansion:  S(Q) = sum_n f(lambda_n / (Q^2 Lambda^2))
    For f(x) = exp(-x), this is Tr exp(-D^2/(Q^2 Lambda^2)).

    Setting tau = 1/(Q^2 Lambda^2), the heat trace on flat T^4 is:
        K(tau) = sum exp(-tau lambda_n) ~ (4pi tau)^{-2} * Vol * [a_0 + a_2 tau + a_4 tau^2 + ...]

    On flat T^4: a_0 = 1, a_2 = R/6 = 0, a_4 = 0 (flat).
    But S depends on Q via tau = 1/(Q^2 Lambda^2), so:
        S(Q) = K(1/(Q^2 Lambda^2))

    The Q-dependence comes purely from the tau(Q) substitution.
    dS/dQ = dK/dtau * dtau/dQ = dK/dtau * (-2/(Q^3 Lambda^2))

    Let's just compute d ln S / d ln Q numerically at Q=1 for various Lambda.
    """
    f_func = lambda x: np.exp(-x)

    Lambda_values = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    print(f"\n{'='*60}")
    print(f"  HEAT-KERNEL COEFFICIENT EXTRACTION")
    print(f"  L={L}, N_max={N_max_lattice}")
    print(f"{'='*60}")
    print(f"  {'Lambda':>8s}  {'S(1)':>14s}  {'d ln S/d ln Q':>14s}  {'d2 ln S/(d ln Q)^2':>20s}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*14}  {'-'*20}")

    for Lam in Lambda_values:
        def S_of_Q(Q):
            return spectral_action_T4_fast(Q, Lam, L, N_max_lattice, f_func)

        def lnS_of_lnQ(lnQ):
            Q = np.exp(lnQ)
            return np.log(S_of_Q(Q))

        S1 = S_of_Q(1.0)

        # d ln S / d ln Q at Q=1 (i.e., ln Q = 0)
        h = 1e-4
        dlnS_dlnQ = (lnS_of_lnQ(h) - lnS_of_lnQ(-h)) / (2*h)
        d2lnS_dlnQ2 = (lnS_of_lnQ(h) - 2*lnS_of_lnQ(0) + lnS_of_lnQ(-h)) / h**2

        print(f"  {Lam:8.1f}  {S1:14.6e}  {dlnS_dlnQ:14.8f}  {d2lnS_dlnQ2:20.8f}")


# ============================================================
# 8.  Conformal scaling theory check
# ============================================================

def conformal_theory_check(L, N_max):
    """
    On a conformally flat 4-torus with g = Q^2 delta:

    The scalar Laplacian eigenvalues scale as lambda_n / Q^2.
    The Dirac operator eigenvalues scale as lambda_n / Q.
    So D^2 eigenvalues scale as lambda_n^2 / Q^2.

    For heat kernel with parameter t = 1/(Q^2 Lambda^2):
        K(t) = sum exp(-t * (2pi n/L)^2)

    Note that t depends on Q as t = 1/(Q^2 Lambda^2).
    BUT the eigenvalues of D^2 on the rescaled torus are (2pi n/L)^2 / Q^2.
    So the argument of the exponential is:
        (2pi n/L)^2 / Q^2  *  1/Lambda^2  =  (2pi n/L)^2 / (Q^2 Lambda^2)

    This means S(Q) = K( 1/(Q^2 Lambda^2) )  where K(t) = sum exp(-t (2pi n/L)^2).

    By Poisson summation / Jacobi theta:
        K(t) = [theta_3(0, exp(-t (2pi/L)^2))]^4 - 1

    For small t (large Q Lambda):
        theta_3(0, exp(-t (2pi/L)^2)) ~ L/(2pi) * sqrt(pi/t) * [1 + O(exp(-L^2/(4t(2pi)^2)))]
        So K(t) ~ (L^2/(4 t))^2 * [1 + ...]  - 1  ~  (L^2 Q^2 Lambda^2 / 4)^2

    In this regime:  S(Q) ~ (L Lambda Q / 2)^4 * (1/4)  proportional to Q^4.
    So d ln S / d ln Q = 4 in the UV limit.

    Let's verify this and compute the exact scaling.
    """
    f_func = lambda x: np.exp(-x)

    print(f"\n{'='*60}")
    print(f"  CONFORMAL SCALING THEORY CHECK")
    print(f"  L={L}, N_max={N_max}")
    print(f"{'='*60}")

    # Test for various Lambda (from IR to UV)
    Lambda_values = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]

    print(f"\n  {'Lambda':>8s}  {'d ln S/d ln Q':>15s}  {'g_QQ=-S\"\"/S':>15s}  {'d2lnS/dQ2':>15s}")
    print(f"  {'-'*8}  {'-'*15}  {'-'*15}  {'-'*15}")

    for Lam in Lambda_values:
        def S_of_Q(Q):
            return spectral_action_T4_fast(Q, Lam, L, N_max, f_func)

        S0 = S_of_Q(1.0)

        h = 1e-4
        S_prime = (S_of_Q(1.0+h) - S_of_Q(1.0-h)) / (2*h)
        S_dprime = (S_of_Q(1.0+h) - 2*S0 + S_of_Q(1.0-h)) / h**2

        g_QQ = -S_dprime / S0
        d2lnS = S_dprime/S0 - (S_prime/S0)**2

        dlnS_dlnQ = S_prime / S0  # = Q * S'/S at Q=1 = S'/S

        print(f"  {Lam:8.1f}  {dlnS_dlnQ:15.8f}  {g_QQ:15.8f}  {d2lnS:15.8f}")

    # Interpretation
    print(f"\n  INTERPRETATION:")
    print(f"  - d ln S / d ln Q = d (dimension of spectral action).")
    print(f"    In UV limit: should approach 4 (= spacetime dimension).")
    print(f"  - g_QQ = -S''/S is the Fisher metric component.")
    print(f"  - d2 ln S / dQ^2 is the connected 2-point function of the scaling.")


# ============================================================
# 9.  The key test: does g_QQ relate to 2?
# ============================================================

def key_test():
    """
    THE CRITICAL TEST.

    If Sigma = 2 ln Q is the spectral action's conformal response,
    then we expect:
        S(Q) = S(1) * exp(c * ln Q)  for some c,
    or more generally S(Q) ~ Q^alpha with some definite alpha.

    The Fisher metric g_QQ = -S''/S at Q=1.
    For S ~ Q^alpha:  S' = alpha S/Q,  S'' = alpha(alpha-1) S/Q^2
    So g_QQ = -alpha(alpha-1)  at Q=1.

    We need to check what alpha is and whether g_QQ has a universal value.
    """
    print(f"\n{'='*70}")
    print(f"  *** KEY TEST: UNIVERSAL g_QQ VALUE ***")
    print(f"{'='*70}")

    # Test with multiple parameter choices
    test_cases = [
        # (Lambda, L, N_max, m2, label)
        (5.0,  1.0, 12, 0.0, "Baseline: Lambda=5, L=1, N=12"),
        (10.0, 1.0, 12, 0.0, "Higher cutoff: Lambda=10, L=1, N=12"),
        (20.0, 1.0, 12, 0.0, "High cutoff: Lambda=20, L=1, N=12"),
        (5.0,  2.0, 12, 0.0, "Larger torus: Lambda=5, L=2, N=12"),
        (5.0,  0.5, 12, 0.0, "Smaller torus: Lambda=5, L=0.5, N=12"),
        (10.0, 2.0, 12, 0.0, "Lambda=10, L=2, N=12"),
        (5.0,  1.0, 12, 1.0, "With mass m^2=1"),
        (5.0,  1.0, 12, 5.0, "With mass m^2=5"),
        (10.0, 1.0, 15, 0.0, "Larger lattice: Lambda=10, L=1, N=15"),
    ]

    results = []
    for Lambda, L, N_max, m2, label in test_cases:
        r = run_verification(Lambda, L, N_max, m2=m2, label=label)
        results.append(r)

    # Summary table
    print(f"\n{'='*70}")
    print(f"  *** SUMMARY TABLE ***")
    print(f"{'='*70}")
    print(f"  {'Label':>45s}  {'g_QQ':>10s}  {'alpha':>8s}  {'d2lnS':>10s}")
    print(f"  {'-'*45}  {'-'*10}  {'-'*8}  {'-'*10}")

    for (Lambda, L, N_max, m2, label), r in zip(test_cases, results):
        alpha = r['S_prime_over_S']  # = S'/S at Q=1 = d ln S / d ln Q at Q=1
        print(f"  {label:>45s}  {r['g_QQ_fisher']:10.6f}  {alpha:8.4f}  {r['d2_lnS']:10.6f}")

    return results


# ============================================================
# 10. Theta function cross-check
# ============================================================

def theta_crosscheck():
    """
    Cross-check using Jacobi theta function (analytic).
    S(Q) = [theta_3(0, q)]^4 - 1  where q = exp(-(2pi/L)^2 / (Q^2 Lambda^2)).

    d/dQ of S: chain rule through q(Q).
    """
    print(f"\n{'='*60}")
    print(f"  JACOBI THETA CROSS-CHECK")
    print(f"{'='*60}")

    L = 1.0
    N_terms = 200  # terms in theta series

    for Lambda in [5.0, 10.0, 20.0]:
        def S_theta(Q):
            t = (2*np.pi/L)**2 / (Q**2 * Lambda**2)
            q = np.exp(-t)
            # theta_3(0,q) = 1 + 2 sum q^{n^2}
            th = 1.0
            for n in range(1, N_terms+1):
                th += 2 * q**(n**2)
            return th**4 - 1.0

        S1 = S_theta(1.0)
        h = 1e-5
        Sp = (S_theta(1+h) - S_theta(1-h)) / (2*h)
        Spp = (S_theta(1+h) - 2*S1 + S_theta(1-h)) / h**2

        g_QQ = -Spp / S1
        alpha = Sp / S1
        d2lnS = Spp/S1 - (Sp/S1)**2

        print(f"  Lambda={Lambda:5.1f}: S(1)={S1:.6e}, alpha={alpha:.8f}, g_QQ={g_QQ:.8f}, d2lnS/dQ2={d2lnS:.8f}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  PAPER 9: NUMERICAL VERIFICATION OF CONFORMAL a_4 VARIATION")
    print("  Spectral action on flat T^4 with conformal factor Q")
    print("=" * 70)

    # --- Test 1: Key universality test ---
    results = key_test()

    # --- Test 2: Q-dependence scan ---
    Q_scan(Lambda=10.0, L=1.0, N_max=12)

    # --- Test 3: Conformal scaling theory ---
    conformal_theory_check(L=1.0, N_max=12)

    # --- Test 4: Power-law extraction ---
    extract_power_law(Lambda=10.0, L=1.0, N_max=12)

    # --- Test 5: Heat-kernel coefficient extraction ---
    extract_a_coefficients(L=1.0, N_max_lattice=12)

    # --- Test 6: Theta function cross-check ---
    theta_crosscheck()

    # --- Final verdict ---
    print(f"\n{'='*70}")
    print(f"  *** FINAL ANALYSIS ***")
    print(f"{'='*70}")

    # Extract the key numbers from baseline
    r0 = results[0]
    alpha = r0['S_prime_over_S']
    g_QQ = r0['g_QQ_fisher']
    d2lnS = r0['d2_lnS']

    print(f"\n  At Q=1 (baseline Lambda=5, L=1, N=12):")
    print(f"    S'(1)/S(1) = alpha = {alpha:.8f}")
    print(f"    g_QQ = -S''(1)/S(1) = {g_QQ:.8f}")
    print(f"    d^2(ln S)/dQ^2      = {d2lnS:.8f}")
    print(f"")
    print(f"  For S(Q) ~ Q^alpha:")
    print(f"    -S''/S = alpha*(alpha-1) = {alpha*(alpha-1):.8f}")
    print(f"    Actual g_QQ             = {g_QQ:.8f}")
    print(f"    Match: {'YES' if abs(g_QQ - alpha*(alpha-1)) < 0.01 else 'NO'}")
    print(f"")
    print(f"  CHECK g_QQ = 2?   Value = {g_QQ:.8f},  difference from 2 = {g_QQ - 2:.8f}")
    print(f"  CHECK alpha = 4?  Value = {alpha:.8f},  difference from 4 = {alpha - 4:.8f}")
    print(f"")
    print(f"  If alpha ~ 4 (= dim), then g_QQ = alpha*(alpha-1) = 4*3 = 12.")
    print(f"  If alpha ~ 2, then g_QQ = 2*1 = 2.")
    print(f"")
    print(f"  IMPORTANT: The spectral action S(Q) counts eigenvalues below Q*Lambda.")
    print(f"  On a flat d-torus with conformal metric g=Q^2 delta:")
    print(f"    D^2 eigenvalues scale as 1/Q^2")
    print(f"    So f(D^2/(Q^2 Lambda^2)) has argument ~ |n|^2 / (Q^4 Lambda^2 L^2/(2pi)^2)")
    print(f"    Wait -- let's be very careful about this!")
    print(f"")
    print(f"  ACTUAL SETUP:")
    print(f"    Eigenvalues of D^2 on (T^4, Q^2 delta): lambda_n = (2pi|n|/L)^2 / Q^2")
    print(f"    Spectral action: S = sum f(lambda_n / Lambda^2)")
    print(f"                       = sum f( (2pi|n|/L)^2 / (Q^2 Lambda^2) )")
    print(f"    This is just K(t) with t = 1/(Q^2 Lambda^2).")
    print(f"    For large Lambda: K(t) ~ Vol * (4 pi t)^{{-d/2}}")
    print(f"                           ~ (L Q)^4 * (Q^2 Lambda^2 / (4 pi))^2")
    print(f"                           ~ Q^8 * (Lambda L)^4 / (4 pi)^2")
    print(f"")
    print(f"  So in the DEEP UV limit, S(Q) ~ Q^8 (!), giving alpha=8, g_QQ=56.")
    print(f"  But this counts Vol(g) ~ Q^4 times the density ~ (Q^2 Lambda^2)^2 = Q^4.")
    print(f"  In practice with finite lattice, the power depends on Lambda*L.")
    print(f"")
    print(f"  The DENSITY (per unit coordinate volume) S/Vol(g) = S/Q^4 should scale as Q^4.")
    print(f"  The INTEGRAND density scales as (Q^2 Lambda^2)^2 ~ Q^4.")
    print(f"")
    print(f"  The PHYSICALLY meaningful conformal variation for Sigma is:")
    print(f"    delta S / delta (ln Q) = alpha * S")
    print(f"    The a_4 coefficient controls the Q-INDEPENDENT part (log divergence).")
    print(f"")

    print(f"\n  {'='*60}")
    print(f"  *** BOTTOM LINE ***")
    print(f"  {'='*60}")
    print(f"  Numerical alpha (= d ln S / d ln Q at Q=1) = {alpha:.6f}")
    print(f"  g_QQ (= -S''/S at Q=1) = {g_QQ:.6f}")
    print(f"  This is a NUMERICAL FACT independent of any theory.")
    print(f"  {'='*60}")
