#!/usr/bin/env python3
"""
R_τ Calculator for Superconductor Material Screening
=====================================================

Framework: τ_eff = τ_phonon - I_Cooper
Room-temperature SC criterion: R_τ = I_Cooper / τ_phonon(300K) ≥ 1

Based on Paper 1: τ = 1 - F, Petz recovery framework.

Usage:
    python rtau_calculator.py                    # Run built-in database
    python rtau_calculator.py --custom           # Interactive mode
    python rtau_calculator.py --scan             # Design space scan

Author: Sheng-Kai Huang (akai@fawstudio.com)
Date: 2026-03-28
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional
import sys


# ============================================================
# Constants
# ============================================================

k_B_eV = 8.617333262e-5   # eV/K
k_B_meV = k_B_eV * 1e3    # meV/K


# ============================================================
# Data structures
# ============================================================

@dataclass
class Material:
    """Superconductor material parameters."""
    name: str
    lam: float          # electron-phonon coupling constant λ (dimensionless)
    theta_D: float      # Debye temperature Θ_D (K)
    T_c: float          # critical temperature (K)
    Delta: float        # superconducting gap Δ(0) (meV)
    N0: float           # DOS at Fermi level N(0) (states/eV/spin/f.u.)
    pressure: float = 0.0   # pressure (GPa), 0 = ambient
    notes: str = ""


# ============================================================
# Known material database
# ============================================================

MATERIALS = [
    Material("Nb",        lam=0.82, theta_D=276,  T_c=9.3,  Delta=1.55, N0=0.94,
             notes="Conventional BCS, elemental"),
    Material("MgB2",      lam=0.87, theta_D=750,  T_c=39,   Delta=7.1,  N0=0.35,
             notes="Two-gap SC, σ-band gap"),
    Material("YBCO",      lam=2.0,  theta_D=410,  T_c=93,   Delta=25.0, N0=1.50,
             notes="d-wave cuprate, effective λ"),
    Material("Bi2212",    lam=1.5,  theta_D=300,  T_c=85,   Delta=30.0, N0=1.20,
             notes="d-wave cuprate"),
    Material("FeSe/STO",  lam=1.0,  theta_D=300,  T_c=65,   Delta=15.0, N0=0.80,
             notes="Interfacial SC, monolayer"),
    Material("H3S",       lam=2.19, theta_D=1500, T_c=203,  Delta=30.0, N0=0.40,
             pressure=150, notes="Confirmed 2015, Drozdov et al."),
    Material("LaH10",     lam=2.5,  theta_D=1100, T_c=260,  Delta=40.0, N0=0.50,
             pressure=170, notes="Confirmed 2019, Somayazulu et al."),
    Material("CaH6",      lam=2.7,  theta_D=1200, T_c=235,  Delta=35.0, N0=0.45,
             pressure=172, notes="Predicted, clathrate hydride"),
    Material("YH6",       lam=2.5,  theta_D=1300, T_c=220,  Delta=33.0, N0=0.48,
             pressure=160, notes="Predicted, clathrate hydride"),
]


# ============================================================
# Core physics functions
# ============================================================

def tau_phonon(T: float, lam: float, theta_D: float) -> float:
    """
    Compute phonon decoherence τ_phonon(T).

    Parameters
    ----------
    T : float
        Temperature (K)
    lam : float
        Electron-phonon coupling constant
    theta_D : float
        Debye temperature (K)

    Returns
    -------
    float
        Dimensionless phonon decoherence parameter

    Notes
    -----
    High-T (T > Θ_D): τ = 2πλ(T/Θ_D)²
    Low-T  (T < Θ_D): τ = 2πλ(T/Θ_D)⁵
    Interpolation: τ = 2πλ(T/Θ_D)² × [1 + (Θ_D/T - 1)³ × f(T/Θ_D)]
    where f smoothly interpolates between the regimes.
    """
    x = T / theta_D

    if x >= 1.0:
        # High-T limit: classical phonons
        return 2 * np.pi * lam * x**2
    else:
        # Smooth interpolation between T^2 (high-T) and T^5 (low-T)
        # Use the Bloch-Gruneisen-inspired form
        # At x=1: gives 2πλ (correct high-T limit)
        # At x→0: gives 2πλ x^5 (correct low-T limit)
        # Interpolation: τ = 2πλ × x^2 × x^(3(1-x))
        # This smoothly transitions: at x=1, exponent=2; at x=0, exponent→5
        alpha = 2 + 3 * (1 - x)
        return 2 * np.pi * lam * x**alpha


def binary_entropy(p: float) -> float:
    """Binary entropy h(p) = -p ln(p) - (1-p) ln(1-p)."""
    if p <= 0 or p >= 1:
        return 0.0
    return -p * np.log(p) - (1 - p) * np.log(1 - p)


def I_k(xi: float, Delta: float) -> float:
    """
    Per-mode Cooper pair entanglement.

    Parameters
    ----------
    xi : float
        Energy relative to Fermi level (same units as Delta)
    Delta : float
        Superconducting gap (same units as xi)

    Returns
    -------
    float
        Entanglement entropy I_k = h(v_k²)
    """
    E_k = np.sqrt(xi**2 + Delta**2)
    v_k2 = 0.5 * (1.0 - xi / E_k)
    return binary_entropy(v_k2)


def I_Cooper(N0: float, Delta_meV: float, theta_D: float) -> float:
    """
    Total Cooper pair entanglement.

    Parameters
    ----------
    N0 : float
        DOS at Fermi level (states/eV/spin/f.u.)
    Delta_meV : float
        Superconducting gap (meV)
    theta_D : float
        Debye temperature (K) — sets integration cutoff

    Returns
    -------
    float
        Dimensionless total entanglement I_Cooper

    Notes
    -----
    I_Cooper = N(0) × ∫_{-ℏω_D}^{+ℏω_D} I(ξ) dξ
    The integral is computed numerically.
    """
    Delta = Delta_meV  # in meV
    omega_D_meV = k_B_meV * theta_D  # ℏω_D in meV

    # Numerical integration over ξ from -ω_D to +ω_D
    n_points = 10000
    xi_max = omega_D_meV
    xi_arr = np.linspace(-xi_max, xi_max, n_points)
    dxi = xi_arr[1] - xi_arr[0]

    I_arr = np.array([I_k(xi, Delta) for xi in xi_arr])
    # np.trapezoid (numpy >= 2.0) or np.trapz (numpy < 2.0)
    _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    integral = _trapz(I_arr, xi_arr)  # in meV

    # Convert: N(0) in states/eV, integral in meV → multiply by 1e-3
    return N0 * integral * 1e-3


def compute_R_tau(mat: Material, T: float = 300.0) -> dict:
    """
    Compute R_τ for a given material at temperature T.

    Parameters
    ----------
    mat : Material
        Material parameters
    T : float
        Temperature (K), default 300K

    Returns
    -------
    dict with keys: tau_phonon, I_Cooper, R_tau, and intermediates
    """
    tp = tau_phonon(T, mat.lam, mat.theta_D)
    ic = I_Cooper(mat.N0, mat.Delta, mat.theta_D)

    R = ic / tp if tp > 0 else np.inf

    return {
        "material": mat.name,
        "T": T,
        "lambda": mat.lam,
        "theta_D": mat.theta_D,
        "T_c": mat.T_c,
        "Delta_meV": mat.Delta,
        "N0": mat.N0,
        "tau_phonon": tp,
        "I_Cooper": ic,
        "R_tau": R,
        "pressure_GPa": mat.pressure,
    }


# ============================================================
# Design equation solver
# ============================================================

def required_Tc_for_Rtau1(lam: float, theta_D: float, N0: float,
                           T: float = 300.0) -> float:
    """
    Compute the minimum T_c needed for R_τ = 1 at temperature T.

    Assumes BCS relation: Δ = 1.76 k_B T_c

    Parameters
    ----------
    lam : float
        Electron-phonon coupling
    theta_D : float
        Debye temperature (K)
    N0 : float
        DOS (states/eV/spin/f.u.)
    T : float
        Target temperature (K)

    Returns
    -------
    float
        Required T_c (K), or np.inf if impossible
    """
    tp = tau_phonon(T, lam, theta_D)

    # Use binary search to find T_c such that I_Cooper(N0, BCS_Delta(T_c), theta_D) = tp
    # BCS: Δ = 1.76 k_B T_c (in meV)
    # We solve I_Cooper(N0, 1.76 * k_B_meV * T_c, theta_D) = tp

    from scipy.optimize import brentq

    def residual(Tc_trial):
        Delta_trial = 1.76 * k_B_meV * Tc_trial  # meV
        ic = I_Cooper(N0, Delta_trial, theta_D)
        return ic - tp

    # Check if achievable at all (T_c up to 1000K)
    ic_max = I_Cooper(N0, 1.76 * k_B_meV * 1000, theta_D)
    if ic_max < tp:
        T_c_required = np.inf
    else:
        try:
            T_c_required = brentq(residual, 0.1, 5000, xtol=0.1)
        except ValueError:
            T_c_required = np.inf

    return T_c_required


def design_space_scan(T: float = 300.0):
    """
    Scan parameter space and identify regions where R_τ ≥ 1.

    Prints a grid of (Θ_D, λ) combinations showing required T_c.
    """
    print(f"\n{'='*70}")
    print(f"DESIGN SPACE SCAN: Required T_c (K) for R_τ = 1 at T = {T}K")
    print(f"Assuming N(0) = 0.5 states/eV/spin/f.u.")
    print(f"{'='*70}")

    N0 = 0.5  # reasonable assumption

    theta_D_values = [300, 500, 750, 1000, 1200, 1500, 2000]
    lam_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

    # Header
    print(f"\n{'Θ_D \\ λ':>10}", end="")
    for lam in lam_values:
        print(f"  {lam:>8.1f}", end="")
    print()
    print("-" * (10 + 10 * len(lam_values)))

    for theta_D in theta_D_values:
        print(f"{theta_D:>8} K", end="")
        for lam in lam_values:
            T_c_req = required_Tc_for_Rtau1(lam, theta_D, N0, T)
            if T_c_req > 1000:
                print(f"  {'> 1000':>8}", end="")
            elif T_c_req < 0:
                print(f"  {'N/A':>8}", end="")
            else:
                marker = " *" if T_c_req <= 300 else ""
                print(f"  {T_c_req:>6.0f}{marker}", end="")
        print()

    print()
    print("* = T_c ≤ 300K (physically achievable in principle)")
    print("Note: BCS T_c max ≈ Θ_D/5 to Θ_D/3 for strong coupling")


# ============================================================
# Numerical verification of the I_k integral prefactor
# ============================================================

def verify_integral_prefactor():
    """
    Numerically compute ∫ I(ξ) dξ / Δ for typical BCS parameters
    to verify the 1.14 prefactor.
    """
    print("\n" + "="*50)
    print("VERIFICATION: ∫ I(ξ)dξ / Δ for various ℏω_D/Δ ratios")
    print("="*50)

    ratios = [5, 10, 20, 50, 100, 200]

    for ratio in ratios:
        Delta = 1.0  # arbitrary units
        omega_D = ratio * Delta

        n_points = 50000
        xi_arr = np.linspace(-omega_D, omega_D, n_points)
        I_arr = np.array([I_k(xi, Delta) for xi in xi_arr])
        _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        integral = _trapz(I_arr, xi_arr)
        prefactor = integral / Delta

        print(f"  ℏω_D/Δ = {ratio:>5d}  →  ∫I dξ / Δ = {prefactor:.4f}")

    print("\nFor typical BCS: ℏω_D/Δ ≈ 10-100, prefactor C ≈ 2.7-3.1")
    print("Full numerical integration used for all material calculations.")


# ============================================================
# Main execution
# ============================================================

def print_results_table():
    """Compute and display R_τ for all materials in database."""

    T = 300.0

    print("="*90)
    print(f"R_τ MATERIAL SCREENING — Room Temperature ({T}K) Superconductor Criterion")
    print(f"Framework: τ_eff = τ_phonon - I_Cooper;  R_τ = I_Cooper / τ_phonon ≥ 1")
    print("="*90)

    # Header
    print(f"\n{'Material':<12} {'λ':>5} {'Θ_D':>6} {'T_c':>6} {'Δ':>6} {'N(0)':>6} "
          f"{'τ_ph':>10} {'I_Coop':>10} {'R_τ':>10} {'P(GPa)':>7}  Status")
    print("-"*100)

    results = []
    for mat in MATERIALS:
        r = compute_R_tau(mat, T)
        results.append(r)

        # Status indicator
        if r["R_tau"] >= 1.0:
            status = "*** R_τ ≥ 1 ***"
        elif r["R_tau"] >= 0.1:
            status = "approaching"
        elif r["R_tau"] >= 0.01:
            status = "moderate deficit"
        else:
            status = "far"

        print(f"{r['material']:<12} {r['lambda']:>5.2f} {r['theta_D']:>6.0f} "
              f"{r['T_c']:>6.1f} {r['Delta_meV']:>6.1f} {r['N0']:>6.2f} "
              f"{r['tau_phonon']:>10.4f} {r['I_Cooper']:>10.6f} "
              f"{r['R_tau']:>10.4f} {r['pressure_GPa']:>7.0f}  {status}")

    # Summary
    print(f"\n{'='*90}")
    print("SUMMARY")
    print(f"{'='*90}")

    above = [r for r in results if r["R_tau"] >= 1.0]
    below = [r for r in results if r["R_tau"] < 1.0]

    print(f"\nMaterials with R_τ ≥ 1 (room-T SC candidates):")
    for r in sorted(above, key=lambda x: -x["R_tau"]):
        print(f"  {r['material']:<12}  R_τ = {r['R_tau']:.3f}  "
              f"(T_c = {r['T_c']:.0f}K, P = {r['pressure_GPa']:.0f} GPa)")

    print(f"\nMaterials with R_τ < 1:")
    for r in sorted(below, key=lambda x: -x["R_tau"]):
        deficit = 1.0 / r["R_tau"]
        print(f"  {r['material']:<12}  R_τ = {r['R_tau']:.4f}  "
              f"(need {deficit:.0f}× improvement)")

    # Design equation
    print(f"\n{'='*90}")
    print("DESIGN EQUATION")
    print(f"{'='*90}")
    print(f"\nFor R_τ = 1 at {T}K:")
    print(f"  N(0) × Δ ≥ τ_phonon({T}K) / 1.14")
    print(f"  Equivalently: N(0) × T_c × Θ_D^α / λ ≥ constant")
    print(f"  where α ≈ 2 (T > Θ_D) or α ≈ 5 (T ≪ Θ_D)")
    print(f"\nRequired T_c for ambient-pressure candidates (assuming λ, Θ_D, N(0)):")

    ambient_candidates = [
        ("MgB2-engineered",   1.2, 900,  0.50),
        ("B-doped diamond",   0.5, 1860, 0.30),
        ("Li@graphene",       1.0, 1500, 0.80),
        ("CaC6-family",       0.8, 1200, 0.60),
        ("LaBeH8 (low-P)",    2.0, 1000, 0.50),
        ("BC2N compound",     1.0, 1400, 0.40),
        ("Metallic H (meta)", 2.5, 2000, 0.60),
    ]

    print(f"\n  {'Candidate':<22} {'λ':>5} {'Θ_D':>6} {'N(0)':>6} {'T_c needed':>10}  Feasible?")
    print(f"  {'-'*70}")
    for name, lam, theta_D, N0 in ambient_candidates:
        Tc_req = required_Tc_for_Rtau1(lam, theta_D, N0, T)
        Tc_max_est = theta_D / 4  # rough strong-coupling upper bound
        feasible = "YES" if Tc_req <= Tc_max_est else "STRETCH" if Tc_req <= theta_D else "NO"
        print(f"  {name:<22} {lam:>5.1f} {theta_D:>6} {N0:>6.2f} {Tc_req:>8.0f} K  {feasible}")

    return results


def interactive_mode():
    """Interactive mode: user inputs parameters, gets R_τ."""
    print("\n" + "="*50)
    print("INTERACTIVE R_τ CALCULATOR")
    print("="*50)
    print("Enter material parameters (or 'q' to quit):\n")

    while True:
        try:
            name = input("Material name: ").strip()
            if name.lower() == 'q':
                break

            lam = float(input("  λ (electron-phonon coupling): "))
            theta_D = float(input("  Θ_D (Debye temperature, K): "))
            T_c = float(input("  T_c (critical temperature, K): "))
            Delta = float(input("  Δ (gap, meV) [0 = use BCS]: "))
            N0 = float(input("  N(0) (DOS, states/eV/spin/f.u.): "))

            if Delta == 0:
                Delta = 1.76 * k_B_meV * T_c
                print(f"  → Using BCS: Δ = {Delta:.2f} meV")

            mat = Material(name, lam, theta_D, T_c, Delta, N0)
            r = compute_R_tau(mat)

            print(f"\n  RESULT: τ_phonon(300K) = {r['tau_phonon']:.6f}")
            print(f"          I_Cooper       = {r['I_Cooper']:.6f}")
            print(f"          R_τ            = {r['R_tau']:.4f}")
            if r['R_tau'] >= 1:
                print(f"          → ROOM-T SC CANDIDATE!")
            else:
                print(f"          → Need {1/r['R_tau']:.1f}× improvement")
            print()

        except (ValueError, EOFError):
            print("Invalid input. Try again or 'q' to quit.\n")


def main():
    if "--custom" in sys.argv:
        interactive_mode()
    elif "--scan" in sys.argv:
        design_space_scan()
    elif "--verify" in sys.argv:
        verify_integral_prefactor()
    else:
        results = print_results_table()
        verify_integral_prefactor()
        design_space_scan()


if __name__ == "__main__":
    main()
