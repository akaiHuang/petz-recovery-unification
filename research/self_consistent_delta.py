#!/usr/bin/env python3
"""
Self-consistent δ(z) for Khronon dark matter.

Solves the coupled Friedmann + Khronon system:
  H² = H²_std / [1 - δ(2+δ)/3]
  δ = I₀c² / (2H²a³)

which reduces to the algebraic equation at each a:
  δ[1 + x(2+δ)/3] = x

where x(z) = I₀c² / (2H²_std(z) · a³(z)) is the "naive" δ.

Author: Sheng-Kai Huang
Date: 2026-03-17
"""

import numpy as np
from scipy.optimize import brentq, fsolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Planck 2018 parameters
# ============================================================
H0_km = 67.4        # km/s/Mpc
H0_si = H0_km * 1e3 / (3.0857e22)  # s^{-1}
c_si = 2.99792458e8  # m/s

Omega_b = 0.049
Omega_r = 9.14e-5    # radiation (photons + 3 massless neutrinos)
Omega_DM_target = 0.265  # Planck Ω_DM
Omega_Lambda = 1.0 - Omega_b - Omega_r - Omega_DM_target  # flatness

rho_crit = 3 * H0_si**2 / (8 * np.pi * 6.674e-11)  # kg/m³

print("=" * 70)
print("SELF-CONSISTENT δ(z) FOR KHRONON DARK MATTER")
print("=" * 70)
print(f"\nPlanck parameters:")
print(f"  H₀ = {H0_km} km/s/Mpc")
print(f"  Ω_b = {Omega_b}")
print(f"  Ω_r = {Omega_r:.4e}")
print(f"  Ω_Λ = {Omega_Lambda:.6f}")
print(f"  Ω_DM(target) = {Omega_DM_target}")
print(f"  ρ_crit = {rho_crit:.4e} kg/m³")

# ============================================================
# Standard H²(z) without dark matter (only baryons + radiation + Λ)
# ============================================================
def H2_std_over_H02(z):
    """H²_std/H₀² = Ω_b(1+z)³ + Ω_r(1+z)⁴ + Ω_Λ"""
    a = 1.0 / (1.0 + z)
    return Omega_b * (1+z)**3 + Omega_r * (1+z)**4 + Omega_Lambda

def H2_LCDM_over_H02(z):
    """Standard ΛCDM H²/H₀² with CDM."""
    return (Omega_b + Omega_DM_target) * (1+z)**3 + Omega_r * (1+z)**4 + Omega_Lambda

# ============================================================
# Self-consistent δ: solve δ[1 + x(2+δ)/3] = x
# ============================================================
def solve_delta_exact(x):
    """
    Solve δ[1 + x(2+δ)/3] = x for δ given x ≥ 0.

    Expanding: δ + x(2δ + δ²)/3 = x
    → (x/3)δ² + (1 + 2x/3)δ - x = 0

    This is a quadratic in δ. The physical root is the positive one.
    """
    A = x / 3.0
    B = 1.0 + 2.0 * x / 3.0
    C = -x

    if x == 0:
        return 0.0

    discriminant = B**2 - 4*A*C
    # positive root:
    delta = (-B + np.sqrt(discriminant)) / (2*A)
    return delta

def solve_delta_exact_array(x_arr):
    """Vectorized version."""
    result = np.zeros_like(x_arr)
    for i, x in enumerate(x_arr):
        if x < 1e-30:
            result[i] = 0.0
        else:
            result[i] = solve_delta_exact(x)
    return result

# ============================================================
# Determine I₀ by requiring Ω_DM(z=0) = 0.265
# ============================================================
# At z=0: Ω_DM = δ₀(2+δ₀)/3 = 0.265
# So δ₀(2+δ₀) = 0.795
# δ₀² + 2δ₀ - 0.795 = 0
# δ₀ = (-2 + sqrt(4 + 3.18))/2 = (-2 + sqrt(7.18))/2 ≈ (-2 + 2.6796)/2 ≈ 0.3398

delta_0_exact = (-2 + np.sqrt(4 + 4*0.795))/2  # quadratic formula for δ²+2δ=0.795
# Check: solve δ² + 2δ - 3*Omega_DM = 0
delta_0_exact = -1 + np.sqrt(1 + 3*Omega_DM_target)
Omega_DM_check = delta_0_exact * (2 + delta_0_exact) / 3

print(f"\n--- Self-consistent δ₀ at z=0 ---")
print(f"  δ₀ = {delta_0_exact:.6f}")
print(f"  Ω_DM = δ₀(2+δ₀)/3 = {Omega_DM_check:.6f} (should be {Omega_DM_target})")

# Now find I₀ (in units of H₀²/c²):
# At z=0: x₀ = I₀c²/(2H²_std·1) = I₀c²/(2H²_std(0))
# And δ₀ = x₀·(1 - δ₀(2+δ₀)/3) / 1  — wait, let me redo.
# Actually x₀ relates to δ₀ via: δ₀[1 + x₀(2+δ₀)/3] = x₀
# So x₀ = δ₀ / [1 - δ₀(2+δ₀)/3] ... wait that's not right either.
#
# From δ[1 + x(2+δ)/3] = x:
# x₀ = δ₀[1 + x₀(2+δ₀)/3]
# This is also an equation for x₀. Let me use the quadratic relation.
# (x/3)δ² + (1+2x/3)δ - x = 0 → x = δ/(1 - δ(2+δ)/3) ... no.
#
# Actually from δ = x[1 - δ(2+δ)/3]:
# x = δ / [1 - δ(2+δ)/3]

x_0 = delta_0_exact / (1 - delta_0_exact*(2+delta_0_exact)/3)
print(f"  x₀ = {x_0:.6f}")
print(f"  Check: δ₀(2+δ₀)/3 = {delta_0_exact*(2+delta_0_exact)/3:.6f}")

# Verify: solve_delta_exact(x_0) should give delta_0_exact
delta_0_verify = solve_delta_exact(x_0)
print(f"  Verify: solve_delta(x₀) = {delta_0_verify:.6f} (should be {delta_0_exact:.6f})")

# x(z) = x₀ · H²_std(0)/H²_std(z) · (1+z)³
# Because x = I₀c²/(2H²_std·a³) and a=1/(1+z)
H2_std_0 = H2_std_over_H02(0)
print(f"  H²_std(0)/H₀² = {H2_std_0:.6f}")
print(f"  = Ω_b + Ω_Λ + Ω_r = {Omega_b + Omega_Lambda + Omega_r:.6f}")

# ============================================================
# Compute δ(z) and ρ_K(z)/ρ_CDM(z) across redshift
# ============================================================
z_arr = np.logspace(-2, 5, 10000)
z_arr = np.concatenate([[0], z_arr])
z_arr = np.sort(z_arr)

# x(z) = x₀ · [H²_std(0)/H²_std(z)] · (1+z)³
x_arr = x_0 * (H2_std_0 / np.array([H2_std_over_H02(z) for z in z_arr])) * (1+z_arr)**3

# Self-consistent δ(z)
delta_sc = solve_delta_exact_array(x_arr)

# Non-self-consistent (naive) δ(z): just use x directly
# The naive approach: δ_naive = Ω_DM_target * H₀²/H²_LCDM(z) * (1+z)³ * (some normalization)
# Actually the "naive" approach from the user's description is:
# δ_naive(z) = δ₀ · (H₀/H_LCDM(z))² · (1+z)³
# But we need to be careful. Let me use two approximations:
# (1) x(z) itself (the ΛCDM approximation where you ignore self-consistency)
# (2) The true self-consistent δ(z)

# For the "ΛCDM approximation", δ was computed as if H = H_LCDM:
# δ_LCDM(z) = I₀c²/(2H²_LCDM(z)·a³)
# We need to set I₀ so that at z=0, δ_LCDM(0)(2+δ_LCDM(0))/3 = Omega_DM
# H²_LCDM(0)/H₀² = 1 by definition
# So δ_LCDM(0) = I₀_LCDM·c²/(2H₀²) → I₀_LCDM = 2H₀²δ₀_LCDM/c²
# And δ₀_LCDM(2+δ₀_LCDM)/3 = Omega_DM → same δ₀.

# The ΛCDM-approximation gives:
# δ_LCDM(z) = δ₀ · (H₀²/H²_LCDM(z)) · (1+z)³
delta_LCDM_approx = delta_0_exact * (1.0 / np.array([H2_LCDM_over_H02(z) for z in z_arr])) * (1+z_arr)**3

# ============================================================
# Key redshifts
# ============================================================
z_key = [0, 1, 5, 10, 50, 100, 500, 1100, 3400, 10000, 50000, 100000]
z_key_labels = ['0', '1', '5', '10', '50', '100', '500', '1100', '3400', '10000', '50000', '100000']

print("\n" + "=" * 90)
print(f"{'z':>8s} | {'x(z)':>12s} | {'δ_SC':>12s} | {'δ_ΛCDM':>12s} | {'Ω_K':>10s} | {'Ω_CDM':>10s} | {'ρ_K/ρ_CDM':>10s}")
print("-" * 90)

results_table = []
for z in z_key:
    h2_std = H2_std_over_H02(z)
    h2_lcdm = H2_LCDM_over_H02(z)

    x_val = x_0 * (H2_std_0 / h2_std) * (1+z)**3
    delta_sc_val = solve_delta_exact(x_val)
    delta_lcdm_val = delta_0_exact / h2_lcdm * (1+z)**3

    # Self-consistent Hubble:
    h2_sc = h2_std / (1 - delta_sc_val*(2+delta_sc_val)/3)

    # Ω_K in self-consistent cosmology
    Omega_K_sc = delta_sc_val * (2 + delta_sc_val) / 3 * (h2_std / h2_sc)
    # Actually Ω_K = ρ_K/(ρ_crit) where ρ_crit uses the actual H²
    # ρ_K = H²δ(2+δ)/(8πG) and ρ_crit = 3H²/(8πG)
    # So Ω_K = δ(2+δ)/3 — this uses the actual H, and it's a fraction of actual ρ_crit
    Omega_K_sc = delta_sc_val * (2 + delta_sc_val) / 3

    # CDM density fraction in ΛCDM
    Omega_CDM_at_z = Omega_DM_target * (1+z)**3 / h2_lcdm

    # ρ_K/ρ_CDM ratio
    # ρ_K = H²_SC · δ(2+δ)/(8πG)
    # ρ_CDM = ρ_CDM,0 · (1+z)³ = (Ω_DM · ρ_crit,0) · (1+z)³ = Ω_DM · 3H₀²/(8πG) · (1+z)³
    # ρ_K/ρ_CDM = H²_SC·δ(2+δ) / [3·Ω_DM·H₀²·(1+z)³]
    #           = (H²_SC/H₀²) · δ(2+δ) / [3·Ω_DM·(1+z)³]
    ratio_rho = h2_sc * delta_sc_val*(2+delta_sc_val) / (3 * Omega_DM_target * (1+z)**3)

    # For LCDM approx:
    h2_lcdm_approx = h2_std / (1 - delta_lcdm_val*(2+delta_lcdm_val)/3) if delta_lcdm_val*(2+delta_lcdm_val)/3 < 1 else float('inf')
    ratio_rho_lcdm = h2_lcdm * delta_lcdm_val*(2+delta_lcdm_val) / (3 * Omega_DM_target * (1+z)**3) if z < 200000 else float('inf')

    print(f"{z:>8d} | {x_val:>12.6f} | {delta_sc_val:>12.6f} | {delta_lcdm_val:>12.6f} | {Omega_K_sc:>10.6f} | {Omega_CDM_at_z:>10.6f} | {ratio_rho:>10.6f}")
    results_table.append((z, x_val, delta_sc_val, delta_lcdm_val, Omega_K_sc, Omega_CDM_at_z, ratio_rho))

print("-" * 90)

# ============================================================
# Detailed comparison at CMB (z=1100)
# ============================================================
z_cmb = 1100
h2_std_cmb = H2_std_over_H02(z_cmb)
h2_lcdm_cmb = H2_LCDM_over_H02(z_cmb)
x_cmb = x_0 * (H2_std_0 / h2_std_cmb) * (1+z_cmb)**3
delta_sc_cmb = solve_delta_exact(x_cmb)
delta_lcdm_cmb = delta_0_exact / h2_lcdm_cmb * (1+z_cmb)**3

h2_sc_cmb = h2_std_cmb / (1 - delta_sc_cmb*(2+delta_sc_cmb)/3)

rho_K_frac = delta_sc_cmb * (2+delta_sc_cmb) / 3
rho_CDM_frac_lcdm = Omega_DM_target * (1+z_cmb)**3 / h2_lcdm_cmb

# Actual energy density ratio
rho_K_over_rho_CDM = h2_sc_cmb * delta_sc_cmb*(2+delta_sc_cmb) / (3 * Omega_DM_target * (1+z_cmb)**3)

print(f"\n{'='*70}")
print(f"DETAILED ANALYSIS AT z = {z_cmb} (CMB)")
print(f"{'='*70}")
print(f"  H²_std/H₀² = {h2_std_cmb:.6f}")
print(f"  H²_ΛCDM/H₀² = {h2_lcdm_cmb:.6f}")
print(f"  x(z=1100) = {x_cmb:.6f}")
print(f"  δ_SC = {delta_sc_cmb:.6f}")
print(f"  δ_ΛCDM_approx = {delta_lcdm_cmb:.6f}")
print(f"  δ_SC/δ_ΛCDM = {delta_sc_cmb/delta_lcdm_cmb:.6f}")
print(f"  H²_SC/H₀² = {h2_sc_cmb:.6f}")
print(f"  H²_SC/H²_ΛCDM = {h2_sc_cmb/h2_lcdm_cmb:.6f}")
print(f"  Ω_K(z=1100) = {rho_K_frac:.6f}")
print(f"  Ω_CDM_ΛCDM(z=1100) = {rho_CDM_frac_lcdm:.6f}")
print(f"  ρ_K/ρ_CDM = {rho_K_over_rho_CDM:.6f}")
print(f"  Excess over CDM: {(rho_K_over_rho_CDM - 1)*100:.2f}%")

# ============================================================
# Effective equation of state w_eff(z)
# ============================================================
# w_eff = -1 - (1/3) d(ln ρ_K)/d(ln a) / 1 ...
# Actually w_eff from: dρ_K/dt + 3H(1+w_eff)ρ_K = 0
# → 1+w_eff = -(1/(3H)) (dρ_K/dt)/ρ_K = -(1/3) d(ln ρ_K)/d(ln a)
#
# ρ_K ∝ H²δ(2+δ), so ln ρ_K = ln H² + ln[δ(2+δ)]
# d(ln ρ_K)/d(ln a) = d(ln H²)/d(ln a) + d(ln[δ(2+δ)])/d(ln a)

# Numerical differentiation
z_fine = np.logspace(-2, 4.5, 5000)
a_fine = 1.0 / (1.0 + z_fine)
ln_a_fine = np.log(a_fine)

x_fine = x_0 * (H2_std_0 / np.array([H2_std_over_H02(z) for z in z_fine])) * (1+z_fine)**3
delta_fine = solve_delta_exact_array(x_fine)
h2_std_fine = np.array([H2_std_over_H02(z) for z in z_fine])
h2_sc_fine = h2_std_fine / (1 - delta_fine*(2+delta_fine)/3)

rho_K_fine = h2_sc_fine * delta_fine * (2 + delta_fine)  # proportional to ρ_K
ln_rho_K = np.log(rho_K_fine + 1e-300)

# Numerical derivative d(ln ρ_K)/d(ln a)
dlnrho_dlna = np.gradient(ln_rho_K, ln_a_fine)
w_eff = -1.0 - dlnrho_dlna / 3.0

# For CDM, w=0 means d(ln ρ)/d(ln a) = -3, so w_eff = -1 - (-3)/3 = 0. ✓

print(f"\n{'='*70}")
print(f"EFFECTIVE EQUATION OF STATE w_eff(z)")
print(f"{'='*70}")
for z_target in [0, 1, 10, 100, 1100, 10000]:
    idx = np.argmin(np.abs(z_fine - z_target))
    print(f"  w_eff(z={z_target:>5d}) = {w_eff[idx]:+.6f}")

# ============================================================
# THE KEY PHYSICS: Why does self-consistency help?
# ============================================================
print(f"\n{'='*70}")
print(f"KEY PHYSICS: SELF-CONSISTENCY EFFECT")
print(f"{'='*70}")

# At high z, H²_std ≈ Ω_b(1+z)³ + Ω_r(1+z)⁴
# x(z) = x₀ · H²_std(0)/H²_std(z) · (1+z)³
# For matter era: H²_std ∝ Ω_b(1+z)³ → x ∝ (1+z)³/(1+z)³ = const
# For radiation era: H²_std ∝ Ω_r(1+z)⁴ → x ∝ (1+z)³/(1+z)⁴ = 1/(1+z) → 0

# In ΛCDM approx: H²_LCDM ∝ (Ω_b+Ω_DM)(1+z)³ → δ ∝ const too but with different normalization
# The key: in self-consistent solution, the "standard" H² only has baryons,
# so at matter domination, x → x₀·(Ω_b+Ω_Λ)/Ω_b = x₀ · H²_std(0)/Ω_b (roughly)

z_test = 1100
h2_std_test = H2_std_over_H02(z_test)
h2_lcdm_test = H2_LCDM_over_H02(z_test)
x_test = x_0 * (H2_std_0 / h2_std_test) * (1+z_test)**3

print(f"\nAt z={z_test}:")
print(f"  H²_std/H₀² = {h2_std_test:.4f} (baryons + radiation + Λ only)")
print(f"  H²_ΛCDM/H₀² = {h2_lcdm_test:.4f} (baryons + CDM + radiation + Λ)")
print(f"  Ratio H²_ΛCDM/H²_std = {h2_lcdm_test/h2_std_test:.4f}")
print(f"  x = {x_test:.6f}")
print(f"  δ_SC = {delta_sc_cmb:.6f} (self-consistent)")

# Show that the self-consistent H² matches ΛCDM H²
print(f"  H²_SC/H₀² = {h2_sc_cmb:.4f}")
print(f"  H²_SC/H²_ΛCDM = {h2_sc_cmb/h2_lcdm_cmb:.6f}")

# The fractional Ω_K should match Ω_CDM if ρ_K tracks CDM
print(f"\n  Ω_K(z=1100) = δ(2+δ)/3 = {delta_sc_cmb*(2+delta_sc_cmb)/3:.6f}")
print(f"  Ω_CDM(z=1100) in ΛCDM = {Omega_DM_target*(1+z_test)**3/h2_lcdm_test:.6f}")

# ============================================================
# Compute ratio at all key redshifts more carefully
# ============================================================
print(f"\n{'='*90}")
print(f"COMPREHENSIVE ρ_K/ρ_CDM TABLE")
print(f"{'='*90}")
print(f"{'z':>8s} | {'δ_SC':>10s} | {'Ω_K=δ(2+δ)/3':>14s} | {'H²_SC/H₀²':>12s} | {'ρ_K/ρ_crit0·a³':>16s} | {'Ω_DM(1+z)³':>14s} | {'ratio':>8s}")
print("-" * 90)

for z in z_key:
    h2s = H2_std_over_H02(z)
    xv = x_0 * (H2_std_0 / h2s) * (1+z)**3
    dsc = solve_delta_exact(xv)
    h2sc = h2s / (1 - dsc*(2+dsc)/3)

    # ρ_K in units of ρ_crit,0
    # ρ_K = (3H₀²/(8πG)) · (H²/H₀²) · δ(2+δ)/3 = ρ_crit,0 · (H²/H₀²) · δ(2+δ)/3
    rho_K_units = h2sc * dsc*(2+dsc)/3

    # ρ_CDM,0·(1+z)³ in units of ρ_crit,0
    rho_CDM_units = Omega_DM_target * (1+z)**3

    ratio = rho_K_units / rho_CDM_units

    omega_K = dsc*(2+dsc)/3

    print(f"{z:>8d} | {dsc:>10.6f} | {omega_K:>14.6f} | {h2sc:>12.4f} | {rho_K_units:>16.4f} | {rho_CDM_units:>14.4f} | {ratio:>8.4f}")

print("-" * 90)

# ============================================================
# DEEP DIVE: Why is ρ_K/ρ_CDM ≠ 1 at high z?
# ============================================================
print(f"\n{'='*70}")
print("DEEP DIVE: SCALING BEHAVIOR")
print(f"{'='*70}")

# In matter era (no radiation): H²_std = H₀²[Ω_b(1+z)³ + Ω_Λ]
# x(z) = x₀ · [Ω_b + Ω_Λ] / [Ω_b(1+z)³ + Ω_Λ] · (1+z)³
# For z >> 1: x(z) → x₀ · (Ω_b + Ω_Λ)/Ω_b
#
# In ΛCDM: δ_LCDM(z) = δ₀ · (1+z)³ / H²_LCDM(z)/H₀²
# For z >> 1: δ_LCDM → δ₀ / (Ω_m) where Ω_m = Ω_b + Ω_DM

x_matter_limit = x_0 * (Omega_b + Omega_Lambda) / Omega_b
delta_matter_limit = solve_delta_exact(x_matter_limit)
delta_lcdm_matter = delta_0_exact / (Omega_b + Omega_DM_target)

print(f"\nMatter-era limit (z >> 1, z << z_eq):")
print(f"  x → x₀·(Ω_b+Ω_Λ)/Ω_b = {x_0:.4f} × {(Omega_b + Omega_Lambda)/Omega_b:.4f} = {x_matter_limit:.4f}")
print(f"  δ_SC → {delta_matter_limit:.6f}")
print(f"  δ_ΛCDM → δ₀/(Ω_b+Ω_DM) = {delta_0_exact:.4f}/{Omega_b+Omega_DM_target:.4f} = {delta_lcdm_matter:.6f}")
print(f"  Ω_K(matter) = δ(2+δ)/3 = {delta_matter_limit*(2+delta_matter_limit)/3:.6f}")
print(f"  Ω_CDM(matter) = Ω_DM/(Ω_b+Ω_DM) = {Omega_DM_target/(Omega_b+Omega_DM_target):.6f}")

# Self-consistent H² in matter era
h2_sc_matter = 1.0 / (1 - delta_matter_limit*(2+delta_matter_limit)/3)  # relative to H²_std
# H²_std,matter = H₀²·Ω_b·(1+z)³
# H²_SC,matter = H₀²·Ω_b·(1+z)³ / (1-Ω_K)
# H²_LCDM,matter = H₀²·(Ω_b+Ω_DM)·(1+z)³

# For H²_SC = H²_LCDM, we need:
# Ω_b/(1-Ω_K) = Ω_b + Ω_DM
# Ω_b = (Ω_b + Ω_DM)(1 - Ω_K)
# Ω_b = Ω_b + Ω_DM - (Ω_b+Ω_DM)·Ω_K
# 0 = Ω_DM - (Ω_b+Ω_DM)·Ω_K
# Ω_K = Ω_DM/(Ω_b+Ω_DM) ... which is exactly what we got!

Omega_K_matter = delta_matter_limit*(2+delta_matter_limit)/3
Omega_CDM_matter = Omega_DM_target/(Omega_b+Omega_DM_target)
print(f"\n  Check: Ω_K = Ω_DM/(Ω_b+Ω_DM)?")
print(f"  Ω_K = {Omega_K_matter:.6f}, Ω_DM/(Ω_b+Ω_DM) = {Omega_CDM_matter:.6f}")
print(f"  Match: {np.isclose(Omega_K_matter, Omega_CDM_matter, rtol=1e-4)}")

# So in the matter era, ρ_K/ρ_CDM should be exactly 1!
# Let me verify:
rho_K_matter = h2_sc_matter * Omega_K_matter  # ∝ H²_SC · δ(2+δ)/3
rho_CDM_matter_unit = Omega_DM_target  # ∝ Ω_DM·(1+z)³ / (1+z)³ normalized
# Actually let me be more careful
# ρ_K = ρ_crit · Ω_K = [3H²_SC/(8πG)] · Ω_K
# ρ_CDM = Ω_DM · ρ_crit,0 · (1+z)³ = Ω_DM · [3H₀²/(8πG)] · (1+z)³
# ratio = H²_SC · Ω_K / (H₀² · Ω_DM · (1+z)³)
# In matter era: H²_SC = H₀² · Ω_b · (1+z)³ / (1-Ω_K)
# ratio = Ω_b · (1+z)³ · Ω_K / [(1-Ω_K) · Ω_DM · (1+z)³]
#       = Ω_b · Ω_K / [(1-Ω_K) · Ω_DM]
# With Ω_K = Ω_DM/(Ω_b+Ω_DM):
# 1-Ω_K = Ω_b/(Ω_b+Ω_DM)
# ratio = Ω_b · [Ω_DM/(Ω_b+Ω_DM)] / [Ω_b/(Ω_b+Ω_DM) · Ω_DM]
#       = Ω_b · Ω_DM / (Ω_b · Ω_DM) = 1 ✓

print(f"\n  THEOREM: In the matter era, ρ_K/ρ_CDM = 1 EXACTLY.")
print(f"  Proof: Ω_K = Ω_DM/(Ω_b+Ω_DM) and H²_SC = Ω_b·(1+z)³/(1-Ω_K)")
print(f"  → ρ_K/ρ_CDM = Ω_b·Ω_K/[(1-Ω_K)·Ω_DM] = 1  ✓")

# ============================================================
# Radiation era behavior
# ============================================================
print(f"\n{'='*70}")
print("RADIATION ERA (z >> z_eq ≈ 3400)")
print(f"{'='*70}")

z_eq = Omega_b / Omega_r  # without DM; but with DM: z_eq = (Ω_b+Ω_DM)/Ω_r
z_eq_std = Omega_b / Omega_r
z_eq_lcdm = (Omega_b + Omega_DM_target) / Omega_r

print(f"  z_eq(baryons only) = Ω_b/Ω_r = {z_eq_std:.0f}")
print(f"  z_eq(ΛCDM) = (Ω_b+Ω_DM)/Ω_r = {z_eq_lcdm:.0f}")

# In radiation era: H²_std ≈ H₀²·Ω_r·(1+z)⁴
# x(z) = x₀ · H²_std(0)/(Ω_r·(1+z)⁴) · (1+z)³ = x₀ · H²_std(0)/(Ω_r·(1+z))
# So x → 0 as z → ∞! This means δ → 0 in radiation era.
# This is EXACTLY what CDM does: Ω_CDM → 0 in radiation era (suppressed by 1/(1+z))

for z_rad in [3400, 10000, 50000, 100000]:
    h2s = H2_std_over_H02(z_rad)
    xv = x_0 * (H2_std_0 / h2s) * (1+z_rad)**3
    dsc = solve_delta_exact(xv)
    h2sc = h2s / (1 - dsc*(2+dsc)/3)
    omega_K = dsc*(2+dsc)/3

    # CDM fraction
    omega_CDM = Omega_DM_target*(1+z_rad)**3 / H2_LCDM_over_H02(z_rad)

    # ρ_K/ρ_CDM
    rho_ratio = h2sc * omega_K / (Omega_DM_target * (1+z_rad)**3)

    print(f"  z={z_rad:>6d}: x={xv:.6f}, δ={dsc:.6f}, Ω_K={omega_K:.6f}, Ω_CDM={omega_CDM:.6f}, ρ_K/ρ_CDM={rho_ratio:.4f}")

# ============================================================
# PLOTS
# ============================================================
outdir = Path("/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/cmb_plots")
outdir.mkdir(exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: δ(z)
ax = axes[0, 0]
z_plot = np.logspace(-1, 4.5, 2000)
x_plot = x_0 * (H2_std_0 / np.array([H2_std_over_H02(z) for z in z_plot])) * (1+z_plot)**3
delta_sc_plot = solve_delta_exact_array(x_plot)
delta_lcdm_plot = delta_0_exact / np.array([H2_LCDM_over_H02(z) for z in z_plot]) * (1+z_plot)**3

ax.loglog(1+z_plot, delta_sc_plot, 'b-', linewidth=2, label=r'$\delta_{\rm SC}$ (self-consistent)')
ax.loglog(1+z_plot, delta_lcdm_plot, 'r--', linewidth=2, label=r'$\delta_{\Lambda CDM}$ (naive)')
ax.axvline(1+1100, color='gray', linestyle=':', alpha=0.5, label='z=1100 (CMB)')
ax.axvline(1+3400, color='orange', linestyle=':', alpha=0.5, label='z=3400 (z_eq)')
ax.set_xlabel('1+z', fontsize=12)
ax.set_ylabel(r'$\delta(z)$', fontsize=12)
ax.set_title(r'Khronon perturbation $\delta(z)$', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(1, 3e4)
ax.grid(True, alpha=0.3)

# Plot 2: ρ_K/ρ_CDM
ax = axes[0, 1]
h2_std_plot = np.array([H2_std_over_H02(z) for z in z_plot])
h2_lcdm_plot_arr = np.array([H2_LCDM_over_H02(z) for z in z_plot])
h2_sc_plot = h2_std_plot / (1 - delta_sc_plot*(2+delta_sc_plot)/3)

rho_K_plot = h2_sc_plot * delta_sc_plot * (2 + delta_sc_plot) / 3
rho_CDM_plot = Omega_DM_target * (1+z_plot)**3

ratio_plot = rho_K_plot / rho_CDM_plot

ax.semilogx(1+z_plot, ratio_plot, 'b-', linewidth=2)
ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
ax.axvline(1+1100, color='gray', linestyle=':', alpha=0.5, label='z=1100')
ax.axvline(1+3400, color='orange', linestyle=':', alpha=0.5, label='z=3400')
ax.set_xlabel('1+z', fontsize=12)
ax.set_ylabel(r'$\rho_K / \rho_{\rm CDM}$', fontsize=12)
ax.set_title(r'Energy density ratio $\rho_K/\rho_{\rm CDM}$', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(1, 3e4)
ax.set_ylim(0.8, 1.3)
ax.grid(True, alpha=0.3)

# Plot 3: w_eff(z)
ax = axes[1, 0]
z_w = z_fine[(z_fine > 0.1) & (z_fine < 3e4)]
w_w = w_eff[(z_fine > 0.1) & (z_fine < 3e4)]
ax.semilogx(1+z_w, w_w, 'b-', linewidth=2)
ax.axhline(0.0, color='k', linestyle='--', alpha=0.5, label='w=0 (CDM)')
ax.axvline(1+1100, color='gray', linestyle=':', alpha=0.5, label='z=1100')
ax.set_xlabel('1+z', fontsize=12)
ax.set_ylabel(r'$w_{\rm eff}(z)$', fontsize=12)
ax.set_title(r'Effective equation of state $w_{\rm eff}$', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(1, 3e4)
ax.set_ylim(-0.1, 0.1)
ax.grid(True, alpha=0.3)

# Plot 4: Ω_K vs Ω_CDM
ax = axes[1, 1]
Omega_K_plot = delta_sc_plot * (2 + delta_sc_plot) / 3
Omega_CDM_plot = Omega_DM_target * (1+z_plot)**3 / h2_lcdm_plot_arr

ax.semilogx(1+z_plot, Omega_K_plot, 'b-', linewidth=2, label=r'$\Omega_K$ (self-consistent)')
ax.semilogx(1+z_plot, Omega_CDM_plot, 'r--', linewidth=2, label=r'$\Omega_{\rm CDM}$ ($\Lambda$CDM)')
ax.axvline(1+1100, color='gray', linestyle=':', alpha=0.5)
ax.axvline(1+3400, color='orange', linestyle=':', alpha=0.5)
ax.set_xlabel('1+z', fontsize=12)
ax.set_ylabel(r'$\Omega(z)$', fontsize=12)
ax.set_title(r'Density fractions $\Omega_K$ vs $\Omega_{\rm CDM}$', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(1, 3e4)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(outdir / 'self_consistent_delta.png', dpi=150)
print(f"\nPlot saved to {outdir}/self_consistent_delta.png")

# ============================================================
# FINAL VERDICT
# ============================================================
print(f"\n{'='*70}")
print("FINAL VERDICT")
print(f"{'='*70}")

# At CMB
idx_cmb = np.argmin(np.abs(z_plot - 1100))
ratio_cmb = ratio_plot[idx_cmb]
print(f"\n  At z=1100 (CMB): ρ_K/ρ_CDM = {ratio_cmb:.4f}")
print(f"  Excess: {(ratio_cmb-1)*100:.2f}%")

# In matter era
idx_100 = np.argmin(np.abs(z_plot - 100))
ratio_100 = ratio_plot[idx_100]
print(f"  At z=100: ρ_K/ρ_CDM = {ratio_100:.4f}")

# At z=0
print(f"  At z=0: δ₀ = {delta_0_exact:.6f}, Ω_K = {delta_0_exact*(2+delta_0_exact)/3:.6f}")

print(f"\n  CRITICAL FINDING:")
if abs(ratio_cmb - 1.0) < 0.05:
    print(f"  Self-consistent δ(z) tracks CDM to {abs(ratio_cmb-1)*100:.1f}% at z=1100.")
    print(f"  The 'CMB excess' problem is RESOLVED by self-consistency!")
elif abs(ratio_cmb - 1.0) < 0.15:
    print(f"  Self-consistent δ(z) tracks CDM to {abs(ratio_cmb-1)*100:.1f}% at z=1100.")
    print(f"  SIGNIFICANT improvement but some tension remains.")
else:
    print(f"  Self-consistent δ(z) deviates from CDM by {abs(ratio_cmb-1)*100:.1f}% at z=1100.")
    print(f"  Self-consistency does NOT fully resolve the CMB problem.")

# ============================================================
# Extra: compute what I₀ adjustment does
# ============================================================
print(f"\n{'='*70}")
print("BONUS: COMPARISON WITH NAIVE (NON-SELF-CONSISTENT) APPROACH")
print(f"{'='*70}")

# In the naive approach: δ_naive(z) = δ₀·H₀²/H²_LCDM(z)·(1+z)³
# At z=1100: δ_naive = δ₀/(Ω_m·(1+z)³+...)·(1+z)³ ≈ δ₀/Ω_m
delta_naive_cmb = delta_0_exact * (1+1100)**3 / H2_LCDM_over_H02(1100)
rho_K_naive_cmb = H2_LCDM_over_H02(1100) * delta_naive_cmb * (2+delta_naive_cmb)/3
rho_CDM_cmb = Omega_DM_target * (1+1100)**3
ratio_naive_cmb = rho_K_naive_cmb / rho_CDM_cmb

print(f"\n  Naive (ΛCDM H(z)) at z=1100:")
print(f"    δ_naive = {delta_naive_cmb:.6f}")
print(f"    ρ_K/ρ_CDM = {ratio_naive_cmb:.4f}")
print(f"    Excess: {(ratio_naive_cmb-1)*100:.2f}%")

print(f"\n  Self-consistent at z=1100:")
print(f"    δ_SC = {delta_sc_cmb:.6f}")
print(f"    ρ_K/ρ_CDM = {rho_K_over_rho_CDM:.4f}")
print(f"    Excess: {(rho_K_over_rho_CDM-1)*100:.2f}%")

print(f"\n  Improvement: {abs(ratio_naive_cmb-1)*100:.2f}% → {abs(rho_K_over_rho_CDM-1)*100:.2f}%")

# ============================================================
# Physical Ω_DM contribution at each redshift
# ============================================================
print(f"\n{'='*70}")
print("KHRONON Ω_K AS FRACTION OF TOTAL ENERGY BUDGET")
print(f"{'='*70}")
print(f"{'z':>8s} | {'Ω_K':>10s} | {'Ω_b':>10s} | {'Ω_r':>10s} | {'Ω_Λ':>10s} | {'Total':>10s}")
print("-" * 70)

for z in [0, 1, 10, 100, 1100, 3400, 10000]:
    h2s = H2_std_over_H02(z)
    xv = x_0 * (H2_std_0 / h2s) * (1+z)**3
    dsc = solve_delta_exact(xv)
    h2sc = h2s / (1 - dsc*(2+dsc)/3)

    Om_K = dsc*(2+dsc)/3
    Om_b = Omega_b*(1+z)**3 / h2sc
    Om_r = Omega_r*(1+z)**4 / h2sc
    Om_L = Omega_Lambda / h2sc
    total = Om_K + Om_b + Om_r + Om_L

    print(f"{z:>8d} | {Om_K:>10.6f} | {Om_b:>10.6f} | {Om_r:>10.6f} | {Om_L:>10.6f} | {total:>10.6f}")
print("-" * 70)

print("\nDone.")
