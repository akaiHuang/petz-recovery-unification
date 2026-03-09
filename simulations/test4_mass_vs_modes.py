"""
Test 4: tau vs Penrose -- Mass Same, Internal DOF Different (DISTINGUISHING)

Our framework: tau depends on N (internal modes). Penrose: collapse depends
on E_G only (gravitational self-energy of the superposition).

Simulation:
  Object A: mass m, N_A = small number of internal modes (simple crystal)
  Object B: mass m, N_B = large number of internal modes (complex molecule)
  Same spatial superposition size, same gravitational potential
  E_G is the SAME for both (same mass distribution)

Our prediction: tau_B >> tau_A (more modes = faster decoherence)
Penrose prediction: tau_A = tau_B (same E_G = same collapse time)

This is a GENUINELY distinguishing prediction.

Reference: Huang (2026) vs Penrose (1996).
"""

import os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# Physical constants (SI)
# ============================================================

HBAR = 1.054571817e-34  # J*s
KB = 1.380649e-23  # J/K
C = 2.99792458e8  # m/s
G_GRAV = 6.67430e-11  # m^3/(kg*s^2)
AMU = 1.66053906660e-27  # kg


# ============================================================
# Penrose model
# ============================================================

def penrose_collapse_time(mass, delta_x):
    """
    Penrose collapse time: t_P = hbar / E_G
    E_G = G * m^2 / delta_x  (gravitational self-energy for point-like mass)
    """
    E_G = G_GRAV * mass**2 / delta_x
    t_P = HBAR / E_G
    return t_P, E_G


def penrose_tau(t, t_P):
    """
    Penrose collapse: tau_P(t) = 1 - exp(-t/t_P)
    (exponential collapse model)
    """
    return 1 - np.exp(-t / t_P)


# ============================================================
# Pikovski-tau model
# ============================================================

def pikovski_decoherence_rate(mass, delta_x, T, N_modes):
    """
    Pikovski gravitational dephasing rate per mode.

    Delta_phi = G * m * delta_x / (R^2 * c^2)  -- but for a uniform field,
    the relevant quantity is:

    Gamma = (omega_0 * Delta_phi / c^2)^2 * n_th

    For simplicity, we use the Pikovski scaling:
    Gamma_total = N * Gamma_per_mode
    Gamma_per_mode = (mass * delta_x * G / (c^2 * hbar))^2 * (k_B * T / hbar)

    Actually, the correct scaling from Pikovski et al. 2015:
    The dephasing rate for each internal mode k at frequency omega_k is:
    Gamma_k = (omega_k * g * delta_x / c^2)^2 * coth(hbar*omega_k/(2*k_B*T))

    where g is the gravitational acceleration. For thermal modes:
    Gamma_k ~ (omega_k * g * delta_x / c^2)^2 * (2*k_B*T / (hbar*omega_k))
            = (g * delta_x / c^2)^2 * omega_k * 2*k_B*T / hbar

    For simplicity, we use a dimensionless model:
    Total decoherence = exp(-N * Gamma * t^2)
    where Gamma encodes all the physics (same for both objects since same mass).
    """
    # Gravitational phase difference per mode
    # For gravitational time dilation: Delta_phi/c^2 = g*delta_x/c^2
    g_acc = G_GRAV * mass / delta_x  # self-gravity (rough estimate)

    # For external field (e.g., Earth's gravity):
    g_earth = 9.81  # m/s^2

    # Effective gravitational potential difference
    Delta_phi = g_earth * delta_x  # V = g*h

    # Typical internal frequency
    omega_typical = KB * T / HBAR  # thermal frequency

    # Dephasing rate per mode (Pikovski scaling)
    Gamma_per_mode = (omega_typical * Delta_phi / C**2)**2

    # Total rate
    Gamma_total = N_modes * Gamma_per_mode

    return Gamma_per_mode, Gamma_total


def pikovski_tau(t, N_modes, Gamma_per_mode):
    """
    tau(t) = (1 - |D(t)|) / 2
    where D(t) = exp(-N * Gamma * t^2) in the Gaussian approximation.
    """
    D = np.exp(-N_modes * Gamma_per_mode * t**2)
    return (1 - D) / 2


def pikovski_decoherence_time(N_modes, Gamma_per_mode):
    """t_c = 1 / sqrt(N * Gamma)"""
    return 1.0 / np.sqrt(N_modes * Gamma_per_mode)


# ============================================================
# Dimensionless model for clear comparison
# ============================================================

def run_test4():
    """Run Test 4: tau vs Penrose, same mass, different N."""
    print("=" * 70)
    print("TEST 4: tau vs Penrose -- Same Mass, Different Internal DOF")
    print("=" * 70)
    print()

    # ========================================================
    # Part A: Dimensionless model
    # ========================================================
    print("Part A: Dimensionless Model")
    print("-" * 50)
    print()
    print("  Both objects have the SAME mass m and superposition size delta_x.")
    print("  Therefore E_G = G*m^2/delta_x is IDENTICAL.")
    print()
    print("  Object A: N_A = 10 internal modes (simple crystal)")
    print("  Object B: N_B = 1000 internal modes (complex molecule)")
    print()

    # Dimensionless parameters
    # Set Gamma = 1 (unit of rate per mode)
    Gamma = 1.0

    N_values = {
        'Object A (N=10)': 10,
        'Object B (N=100)': 100,
        'Object C (N=1000)': 1000,
        'Object D (N=10000)': 10000,
    }

    # Penrose: same for all (since same mass)
    # Set t_P = 1 (Penrose timescale)
    t_P = 1.0

    t_max = 3.0
    t_values = np.linspace(0, t_max, 500)

    results = {}

    print(f"  {'Object':<25} {'N':>8} {'t_c (tau)':>12} {'t_P (Penrose)':>14} "
          f"{'t_c/t_P':>10}")
    print("  " + "-" * 70)

    for name, N in N_values.items():
        t_c = pikovski_decoherence_time(N, Gamma)
        tau_t = pikovski_tau(t_values, N, Gamma)
        tau_P = penrose_tau(t_values, t_P)

        results[name] = {
            'N': N,
            'tau': tau_t,
            'tau_penrose': tau_P,
            't_c': t_c,
            't_P': t_P,
        }

        print(f"  {name:<25} {N:>8} {t_c:>12.6f} {t_P:>14.6f} "
              f"{t_c/t_P:>10.6f}")

    print()
    print("  KEY DIFFERENCE:")
    print(f"    Penrose: tau_A(t) = tau_B(t) = tau_C(t) = tau_D(t) for ALL t")
    print(f"    Ours: t_c scales as 1/sqrt(N*Gamma)")
    print(f"    Object A (N=10):    t_c = {results['Object A (N=10)']['t_c']:.6f}")
    print(f"    Object D (N=10000): t_c = {results['Object D (N=10000)']['t_c']:.6f}")
    print(f"    Ratio: {results['Object A (N=10)']['t_c'] / results['Object D (N=10000)']['t_c']:.1f}x")
    print()

    # ========================================================
    # Part B: Physical estimate
    # ========================================================
    print("Part B: Physical Estimate with Real Parameters")
    print("-" * 50)
    print()

    # Example: 1000 amu molecule in 1 micron superposition
    mass = 1000 * AMU  # ~10^-24 kg
    delta_x = 1e-6  # 1 micron
    T = 300  # room temperature

    print(f"  Mass: {mass:.2e} kg ({mass/AMU:.0f} amu)")
    print(f"  Superposition size: {delta_x:.1e} m")
    print(f"  Temperature: {T} K")
    print()

    # Penrose prediction
    t_P_phys, E_G = penrose_collapse_time(mass, delta_x)
    print(f"  Penrose:")
    print(f"    E_G = {E_G:.2e} J")
    print(f"    t_P = hbar/E_G = {t_P_phys:.2e} s")
    print()

    # Pikovski prediction for different N
    physical_objects = {
        'Simple diatomic (N=6)': 6,
        'Small molecule (N=30)': 30,
        'C60 fullerene (N=174)': 174,  # 3*(60-1) - 5 = 174 vibrational modes
        'Protein (N=3000)': 3000,
        'Nanocrystal (N=30000)': 30000,
    }

    print(f"  Pikovski prediction (delta_x = {delta_x:.0e} m, T = {T} K):")
    print(f"  {'Object':<30} {'N':>8} {'Gamma_mode':>14} {'t_c':>14}")
    print("  " + "-" * 70)

    for name, N in physical_objects.items():
        Gamma_mode, Gamma_total = pikovski_decoherence_rate(mass, delta_x, T, N)
        t_c = pikovski_decoherence_time(N, Gamma_mode)
        print(f"  {name:<30} {N:>8} {Gamma_mode:>14.2e} {t_c:>14.2e} s")

    print()
    print("  Note: Pikovski rates are EXTREMELY small for these parameters")
    print("  (decoherence time >> age of universe). The key point is the")
    print("  SCALING: t_c ~ 1/sqrt(N), not the absolute timescale.")
    print()

    # ========================================================
    # Part C: Experimental distinguishability
    # ========================================================
    print("Part C: Experimental Signature")
    print("-" * 50)
    print()
    print("  To distinguish tau framework from Penrose:")
    print("  1. Prepare TWO objects with SAME mass but DIFFERENT internal DOF")
    print("     e.g., C60 (174 modes) vs linear C60-chain (fewer modes)")
    print("  2. Put both in the SAME spatial superposition")
    print("  3. Measure decoherence rate for both")
    print()
    print("  tau prediction: decoherence_rate ~ N * Gamma")
    print("  Penrose prediction: decoherence_rate ~ E_G (same for both)")
    print()
    print("  The ratio of decoherence rates equals the ratio of N values:")
    for name_a, N_a in list(physical_objects.items())[:3]:
        for name_b, N_b in list(physical_objects.items())[3:]:
            ratio = N_b / N_a
            print(f"    {name_b} / {name_a}: N_B/N_A = {ratio:.1f}")
    print()
    print("  A ratio significantly different from 1 would FALSIFY Penrose")
    print("  and CONFIRM the tau framework (or vice versa).")
    print()

    return results, N_values, t_values


if __name__ == "__main__":
    results, N_values, t_values = run_test4()
