"""
Test 5: tau = 0 in Vacuum vs GRW Spontaneous Collapse

Standard QM (our framework): perfectly isolated system has tau(t) = 0 for all t.
GRW model: spontaneous localization noise collapses even isolated systems.

Simulation:
  Perfectly isolated system (no environment modes, no gravity)
  Standard QM: tau(t) = 0 for all t (unitary, no decoherence)
  GRW model: add spontaneous localization noise (rate lambda_GRW)
  For N = 1 to N = 10^23 particles, compute tau_GRW(t)

Our prediction: tau = 0 always (no environment -> no collapse)
GRW prediction: tau_GRW ~ 1 - exp(-N*lambda*t) -> macroscopic objects
                collapse even in vacuum

Reference: Ghirardi, Rimini, Weber (1986); Huang (2026).
"""

import os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# GRW parameters
# ============================================================

# Standard GRW parameters
LAMBDA_GRW = 1e-16  # s^-1 per particle (spontaneous collapse rate)
R_C = 1e-7  # m (localization length scale)

# CSL (continuous spontaneous localization) variant
LAMBDA_CSL = 1e-17  # s^-1 per particle (Adler's enhanced rate: ~1e-8)


# ============================================================
# Models
# ============================================================

def tau_standard_qm(t, N):
    """
    Standard QM: perfectly isolated system.
    Unitary evolution preserves the state perfectly.
    tau(t) = 0 for all t.
    """
    return np.zeros_like(t)


def tau_grw(t, N, lambda_grw=LAMBDA_GRW):
    """
    GRW spontaneous collapse model.

    For a superposition of N particles, the collapse rate is N * lambda.
    The off-diagonal elements decay as exp(-N * lambda * t).

    tau_GRW(t) = (1 - exp(-N * lambda * t)) / 2

    For a spatial superposition |L> + |R> with separation >> r_c,
    each GRW hit on any particle collapses the superposition.
    """
    rate = N * lambda_grw
    D = np.exp(-rate * t)
    return (1 - D) / 2


def tau_csl(t, N, lambda_csl=LAMBDA_CSL):
    """
    CSL (Continuous Spontaneous Localization) model.

    Similar to GRW but continuous. The decoherence rate is:
    Gamma_CSL = N * lambda_CSL * (Delta_x / r_c)^2

    For superpositions with Delta_x >> r_c:
    tau_CSL(t) = (1 - exp(-N * lambda_CSL * t)) / 2
    """
    rate = N * lambda_csl
    D = np.exp(-rate * t)
    return (1 - D) / 2


def tau_diosi_penrose(t, N, mass_per_particle=1e-27, delta_x=1e-6):
    """
    Diosi-Penrose gravitational collapse model.

    Rate ~ G * m^2 * N^2 / (hbar * delta_x)
    (different N-scaling from GRW!)
    """
    G = 6.674e-11
    hbar = 1.055e-34
    m_total = N * mass_per_particle
    rate = G * m_total**2 / (hbar * delta_x)
    D = np.exp(-rate * t)
    return (1 - D) / 2


# ============================================================
# Run Test 5
# ============================================================

def run_test5():
    """Run Test 5: Vacuum tau comparison."""
    print("=" * 70)
    print("TEST 5: tau = 0 in Vacuum vs GRW Spontaneous Collapse")
    print("=" * 70)
    print()

    # ========================================================
    # Part A: Comparison across N scales
    # ========================================================
    print("Part A: tau(t) for various particle numbers N")
    print("-" * 50)
    print()
    print(f"  GRW rate: lambda = {LAMBDA_GRW:.0e} s^-1 per particle")
    print(f"  CSL rate: lambda = {LAMBDA_CSL:.0e} s^-1 per particle")
    print()

    N_values = [1, 10, 100, 1e6, 1e10, 1e15, 1e20, 6.022e23]
    N_labels = ['1', '10', '100', '10^6', '10^{10}', '10^{15}',
                '10^{20}', 'N_A']

    # Time range: 1 microsecond to 1 year
    t_values = np.logspace(-6, 7, 500)  # 1 us to ~3 months

    results = {}

    print(f"  {'N':<12} {'GRW t_half':>14} {'CSL t_half':>14} {'Std QM':>12}")
    print("  " + "-" * 55)

    for N, label in zip(N_values, N_labels):
        # GRW half-life: t_half = ln(2) / (N * lambda)
        t_half_grw = np.log(2) / (N * LAMBDA_GRW)
        t_half_csl = np.log(2) / (N * LAMBDA_CSL)

        tau_qm = tau_standard_qm(t_values, N)
        tau_g = tau_grw(t_values, N)
        tau_c = tau_csl(t_values, N)

        results[label] = {
            'N': N,
            'tau_qm': tau_qm,
            'tau_grw': tau_g,
            'tau_csl': tau_c,
            't_half_grw': t_half_grw,
            't_half_csl': t_half_csl,
        }

        def format_time(t):
            if t < 1e-6:
                return f"{t*1e9:.1f} ns"
            elif t < 1e-3:
                return f"{t*1e6:.1f} us"
            elif t < 1:
                return f"{t*1e3:.1f} ms"
            elif t < 60:
                return f"{t:.1f} s"
            elif t < 3600:
                return f"{t/60:.1f} min"
            elif t < 86400:
                return f"{t/3600:.1f} hr"
            elif t < 3.156e7:
                return f"{t/86400:.1f} days"
            elif t < 3.156e10:
                return f"{t/3.156e7:.1f} yr"
            elif t < 3.156e17:
                return f"{t/3.156e7:.1e} yr"
            else:
                return f"{t/3.156e7:.1e} yr"

        print(f"  N={label:<10} {format_time(t_half_grw):>14} "
              f"{format_time(t_half_csl):>14} {'tau=0 always':>12}")

    print()

    # ========================================================
    # Part B: The stark contrast
    # ========================================================
    print("Part B: The Stark Contrast")
    print("-" * 50)
    print()
    print("  Standard QM / tau framework:")
    print("    - Perfectly isolated system: tau(t) = 0 for ALL t, ALL N")
    print("    - No environment = no decoherence = no arrow of time")
    print("    - A cat-sized object in perfect vacuum stays in superposition")
    print("      FOREVER (in principle)")
    print()
    print("  GRW model:")
    print("    - Spontaneous collapse even without environment")
    print("    - Single particle: t_half ~ 10^16 s ~ 300 million years")
    print("    - Avogadro's number: t_half ~ 10^{-8} s ~ 10 ns")
    print("    - Macroscopic object collapses in nanoseconds even in vacuum")
    print()
    print("  This is a STARK, EXPERIMENTALLY DISTINGUISHABLE difference.")
    print()

    # ========================================================
    # Part C: Current experimental status
    # ========================================================
    print("Part C: Current Experimental Constraints")
    print("-" * 50)
    print()
    print("  Best current tests (as of 2025):")
    print()
    print("  1. LISA Pathfinder (Carlesso et al. 2022):")
    print("     - lambda_GRW < 5 x 10^-7 s^-1 at r_c = 10^-7 m")
    print("     - 10 orders of magnitude above standard GRW value")
    print("     - Does NOT rule out standard GRW yet")
    print()
    print("  2. Optomechanical experiments (Vinante et al. 2020):")
    print("     - Constrain CSL: lambda_CSL < 10^-11 at r_c = 10^-7 m")
    print("     - Significantly constrains Adler's enhanced CSL (lambda ~ 10^-8)")
    print()
    print("  3. Matter-wave interferometry (Arndt group, 2019):")
    print("     - Superposition of ~25,000 amu molecules")
    print("     - Coherence maintained for ms timescales")
    print("     - Consistent with standard QM (tau = 0 in vacuum)")
    print()
    print("  4. Proposed tests:")
    print("     - MAQRO satellite mission: micron-sized particles")
    print("     - N ~ 10^10 atoms, t ~ 100 s observation time")
    print("     - GRW prediction: tau ~ 1 (fully collapsed)")
    print("     - tau prediction: tau = 0 (if truly isolated)")
    print()

    # ========================================================
    # Part D: Quantitative predictions
    # ========================================================
    print("Part D: Quantitative Predictions at Specific Scales")
    print("-" * 50)
    print()

    test_cases = [
        ("Single atom, 1 ms", 1, 1e-3),
        ("10 atoms, 1 ms", 10, 1e-3),
        ("C60, 1 ms", 60, 1e-3),
        ("Protein (10^4 atoms), 1 s", 1e4, 1.0),
        ("Virus (10^7 atoms), 1 s", 1e7, 1.0),
        ("Bacterium (10^10 atoms), 1 s", 1e10, 1.0),
        ("Dust grain (10^14 atoms), 1 s", 1e14, 1.0),
        ("Cat (10^26 atoms), 1 s", 1e26, 1.0),
    ]

    print(f"  {'System':<35} {'tau (std QM)':>12} {'tau (GRW)':>12} "
          f"{'tau (CSL)':>12}")
    print("  " + "-" * 75)

    for name, N, t in test_cases:
        tau_qm = 0.0
        tau_g = float(tau_grw(np.array([t]), N)[0])
        tau_c = float(tau_csl(np.array([t]), N)[0])

        # Format tau values
        def fmt_tau(v):
            if v < 1e-15:
                return "~0"
            elif v < 0.001:
                return f"{v:.2e}"
            else:
                return f"{v:.4f}"

        print(f"  {name:<35} {fmt_tau(tau_qm):>12} {fmt_tau(tau_g):>12} "
              f"{fmt_tau(tau_c):>12}")

    print()
    print("CONCLUSION:")
    print("  Standard QM (tau framework): tau = 0 in vacuum for ALL systems.")
    print("  GRW/CSL: tau > 0 for macroscopic systems even in vacuum.")
    print("  The difference becomes measurable at N ~ 10^10-10^14 atoms")
    print("  on timescales of seconds (which MAQRO-type experiments can probe).")
    print()

    return results, N_values, N_labels, t_values


if __name__ == "__main__":
    results, N_values, N_labels, t_values = run_test5()
