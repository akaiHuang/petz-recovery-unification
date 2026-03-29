#!/usr/bin/env python3
"""
Self-consistent δ(z) — Version 2: Root cause analysis.

The v1 script showed ρ_K/ρ_CDM ≈ 1.17 at z=1100 even with self-consistency.
Here we trace exactly WHY.

Key insight: the "matter era theorem" (ρ_K/ρ_CDM=1) failed because
the self-consistent solution has a DIFFERENT matter-radiation equality
than ΛCDM. We need to understand this thoroughly.

Author: Sheng-Kai Huang
Date: 2026-03-17
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# Parameters
# ============================================================
Omega_b = 0.049
Omega_r = 9.14e-5
Omega_DM = 0.265
Omega_Lambda = 1.0 - Omega_b - Omega_r - Omega_DM

# δ₀ from Ω_DM = δ₀(2+δ₀)/3
delta_0 = -1 + np.sqrt(1 + 3*Omega_DM)  # = 0.339776

# ============================================================
# Core equations
# ============================================================
def H2_std(z):
    """H²/H₀² without DM (baryons + radiation + Λ only)."""
    return Omega_b*(1+z)**3 + Omega_r*(1+z)**4 + Omega_Lambda

def H2_LCDM(z):
    """Standard ΛCDM H²/H₀²."""
    return (Omega_b + Omega_DM)*(1+z)**3 + Omega_r*(1+z)**4 + Omega_Lambda

def solve_delta(x):
    """Solve (x/3)δ² + (1+2x/3)δ - x = 0, positive root."""
    if x < 1e-30:
        return 0.0
    A = x / 3.0
    B = 1.0 + 2*x/3.0
    return (-B + np.sqrt(B**2 + 4*A*x)) / (2*A)

# I₀ parameter: x₀ = δ₀/(1 - δ₀(2+δ₀)/3) = δ₀/(1-Ω_DM)
x_0 = delta_0 / (1 - Omega_DM)
H2_std_0 = H2_std(0)

def x_of_z(z):
    """The naive δ if self-consistency were ignored."""
    return x_0 * H2_std_0 / H2_std(z) * (1+z)**3

def delta_SC(z):
    """Self-consistent δ(z)."""
    return solve_delta(x_of_z(z))

def H2_SC(z):
    """Self-consistent H²/H₀²."""
    d = delta_SC(z)
    return H2_std(z) / (1 - d*(2+d)/3)

# ============================================================
# ROOT CAUSE ANALYSIS
# ============================================================
print("=" * 70)
print("ROOT CAUSE ANALYSIS: WHY ρ_K/ρ_CDM ≠ 1")
print("=" * 70)

print("\n--- The problem is the QUADRATIC term ---")
print("ρ_K ∝ H²δ(2+δ) = H²(2δ + δ²)")
print("CDM would give ρ ∝ H²·Ω_CDM/Ω_m (linear in Ω_DM)")
print("The δ² term makes ρ_K grow FASTER than linear")

print("\n--- Let's check the pure matter era ---")
print(f"x_matter_limit = x₀·H²_std(0)/Ω_b = {x_0}×{H2_std_0}/{Omega_b} = {x_0*H2_std_0/Omega_b:.6f}")
x_mat = x_0 * H2_std_0 / Omega_b
d_mat = solve_delta(x_mat)
Omega_K_mat = d_mat*(2+d_mat)/3
print(f"δ_matter = {d_mat:.6f}")
print(f"Ω_K = δ(2+δ)/3 = {Omega_K_mat:.6f}")
print(f"1-Ω_K = {1-Omega_K_mat:.6f}")
print(f"H²_SC/H²_std = 1/(1-Ω_K) = {1/(1-Omega_K_mat):.6f}")

# H²_SC in matter era: H₀²·Ω_b·(1+z)³/(1-Ω_K)
# H²_LCDM in matter era: H₀²·(Ω_b+Ω_DM)·(1+z)³
# These are equal iff: Ω_b/(1-Ω_K) = Ω_b+Ω_DM
# → 1-Ω_K = Ω_b/(Ω_b+Ω_DM) → Ω_K = Ω_DM/(Ω_b+Ω_DM) = 0.8440

Omega_K_needed = Omega_DM / (Omega_b + Omega_DM)
print(f"\nFor H²_SC = H²_ΛCDM in matter era, need Ω_K = {Omega_K_needed:.6f}")
print(f"Actual Ω_K = {Omega_K_mat:.6f}")
print(f"MISMATCH: {(Omega_K_mat - Omega_K_needed)/Omega_K_needed * 100:.2f}%")

# Why the mismatch? Because x₀ was calibrated at z=0 where Λ matters,
# not in the pure matter era!

# ============================================================
# ALTERNATIVE: What if we DEFINE the conservation law differently?
# ============================================================
print("\n" + "=" * 70)
print("THE FUNDAMENTAL ISSUE")
print("=" * 70)

print("""
The Khronon conservation law (μ²δ = const/a³) combined with μ = H/c gives:
  H²δ = const·(1+z)³

But ρ_K = c²H²δ(2+δ)/(8πG), which has a QUADRATIC δ dependence.

For CDM: ρ_CDM = const·(1+z)³ (exactly linear scaling)

The quadratic δ² term means:
- At z=0: δ₀ ≈ 0.34, so δ²/2δ ≈ 0.17 → 17% correction
- In matter era: δ ≈ 0.90, so δ²/2δ ≈ 0.45 → 45% correction!

The ratio ρ_K/ρ_CDM ≈ (2δ+δ²)/(2δ) ≈ 1 + δ/2 at leading order.
This is ~1.17 in the matter era and ~1.45 asymptotically.
BUT the δ² also changes H², which compensates partially.

The EXACT statement is: ρ_K/ρ_CDM = 1 requires Ω_K = Ω_CDM in the SAME
Hubble flow. Self-consistency changes H, which changes both sides.
""")

# ============================================================
# Let's compute ρ_K/ρ_CDM analytically in the matter era
# ============================================================
print("=" * 70)
print("EXACT MATTER-ERA ANALYSIS")
print("=" * 70)

# In matter era: H²_std = H₀²Ω_b(1+z)³
# Conservation: H²δ = C(1+z)³ where C = H₀²x₀·H2_std_0/2... let me redo.
#
# Actually the conservation law is: (H²/c²)·δ = I₀/(2a³)
# In our units (H₀=1, c=1): H²·δ = (I₀/2)·(1+z)³
#
# At z=0: δ₀ = I₀/2 → I₀ = 2δ₀ (if H₀=1 at z=0... but H²(z=0)/H₀²=1 only in ΛCDM!)
# Self-consistently: H²_SC(0) = H²_std(0)/(1-Ω_DM) = 0.735/0.735 = 1. OK good.
# So I₀ = 2·H²_SC(0)·δ₀ = 2·1·δ₀ = 2δ₀ = 0.6796 (in units of H₀²)

I0 = 2 * delta_0  # because H²_SC(0)/H₀² = 1

print(f"I₀ = 2δ₀ = {I0:.6f}")

# In matter era: H²_std = Ω_b(1+z)³
# Self-consistent: H² = Ω_b(1+z)³/(1-δ(2+δ)/3)
# Conservation: H²·δ = I₀(1+z)³/2 = δ₀(1+z)³
# So: [Ω_b(1+z)³/(1-δ(2+δ)/3)]·δ = δ₀(1+z)³
# Ω_b·δ/(1-δ(2+δ)/3) = δ₀
# δ = δ₀(1-δ(2+δ)/3)/Ω_b
# This is the SAME equation as before with x = δ₀/Ω_b ≈ 6.93

print(f"x_matter = δ₀/Ω_b = {delta_0/Omega_b:.6f}")
print(f"Matches x_0·H2_std_0/Ω_b = {x_0*H2_std_0/Omega_b:.6f}")
# These should match: x_0·H2_std_0 = δ₀/(1-Ω_DM)·(1-Ω_DM) = δ₀ ✓ (since H2_std_0 = 1-Ω_DM)
print(f"x_0·H2_std_0 = {x_0*H2_std_0:.6f} = δ₀ = {delta_0:.6f}")

# Now ρ_K/ρ_CDM in matter era:
# ρ_K = (H₀²/H₀²)·H²·δ(2+δ)/3 = H²/H₀²·δ(2+δ)/3 [in ρ_crit units]
#     = [Ω_b(1+z)³/(1-δ(2+δ)/3)]·δ(2+δ)/3
#
# ρ_CDM = Ω_DM·(1+z)³
#
# ratio = Ω_b·δ(2+δ)/3 / [(1-δ(2+δ)/3)·Ω_DM]
#       = Ω_b·Ω_K / [(1-Ω_K)·Ω_DM]
#
# For this to be 1: Ω_K(1-Ω_K)⁻¹ = Ω_DM/Ω_b
# → Ω_K = Ω_DM/(Ω_b+Ω_DM) ... but this is a condition, not automatic!
#
# With x = δ₀/Ω_b:
# Solve (x/3)δ² + (1+2x/3)δ - x = 0

d_m = solve_delta(delta_0/Omega_b)
Omega_K_m = d_m*(2+d_m)/3
ratio_exact = Omega_b * Omega_K_m / ((1-Omega_K_m)*Omega_DM)

print(f"\nMatter era exact:")
print(f"  δ = {d_m:.6f}")
print(f"  Ω_K = {Omega_K_m:.6f}")
print(f"  ρ_K/ρ_CDM = Ω_b·Ω_K/[(1-Ω_K)·Ω_DM] = {ratio_exact:.6f}")

# Let's see if we can prove this equals 1 analytically.
# x = δ₀/Ω_b, and δ satisfies δ = x(1-δ(2+δ)/3) = x(1-Ω_K)
# So δ = (δ₀/Ω_b)(1-Ω_K)
# And Ω_K = δ(2+δ)/3
# Substitute: Ω_K = [(δ₀/Ω_b)(1-Ω_K)] · [2+(δ₀/Ω_b)(1-Ω_K)] / 3
#
# For ρ_K/ρ_CDM = 1: Ω_K(1-Ω_K)⁻¹ = Ω_DM/Ω_b
# → Ω_K = Ω_DM/(Ω_b+Ω_DM)
# → 1-Ω_K = Ω_b/(Ω_b+Ω_DM)
# → δ = (δ₀/Ω_b)·Ω_b/(Ω_b+Ω_DM) = δ₀/(Ω_b+Ω_DM)
# → Ω_K = [δ₀/(Ω_b+Ω_DM)]·[2+δ₀/(Ω_b+Ω_DM)]/3

delta_check = delta_0/(Omega_b+Omega_DM)
Omega_K_check = delta_check*(2+delta_check)/3
Omega_K_target = Omega_DM/(Omega_b+Omega_DM)

print(f"\n  If ρ_K/ρ_CDM = 1, then:")
print(f"    δ = δ₀/(Ω_b+Ω_DM) = {delta_check:.6f}")
print(f"    Ω_K = δ(2+δ)/3 = {Omega_K_check:.6f}")
print(f"    Need Ω_K = Ω_DM/(Ω_b+Ω_DM) = {Omega_K_target:.6f}")
print(f"    These match iff δ₀(2+δ₀/(Ω_b+Ω_DM)) = 3Ω_DM")

lhs = delta_0 * (2 + delta_0/(Omega_b+Omega_DM))
rhs = 3*Omega_DM
print(f"    LHS = δ₀(2+δ₀/Ω_m) = {lhs:.6f}")
print(f"    RHS = 3Ω_DM = {rhs:.6f}")
print(f"    BUT δ₀(2+δ₀) = 3Ω_DM by definition!")
print(f"    So LHS = δ₀·2 + δ₀²/Ω_m = {delta_0*2:.6f} + {delta_0**2/0.314:.6f}")
print(f"    And RHS = δ₀·2 + δ₀² = {delta_0*2:.6f} + {delta_0**2:.6f}")
print(f"    LHS - RHS = δ₀²(1/Ω_m - 1) = {delta_0**2*(1/0.314 - 1):.6f}")
print(f"    This is > 0, so Ω_K > Ω_DM/(Ω_b+Ω_DM), meaning ρ_K > ρ_CDM ✓")

# The fractional excess:
# ρ_K/ρ_CDM - 1 ≈ δ₀²(1/Ω_m - 1) / (3Ω_DM) ... rough
excess_estimate = delta_0**2 * (1/(Omega_b+Omega_DM) - 1) / (3*Omega_DM)
print(f"\n  Fractional excess ≈ δ₀²(1/Ω_m-1)/(3Ω_DM) ≈ {excess_estimate:.4f}")
print(f"  Actual excess in matter era: {ratio_exact - 1:.4f}")
print(f"  (Rough estimate captures the order of magnitude)")

# ============================================================
# EXACT analytical formula for ρ_K/ρ_CDM in matter era
# ============================================================
print("\n" + "=" * 70)
print("EXACT FORMULA FOR ρ_K/ρ_CDM IN MATTER ERA")
print("=" * 70)

# Let's just compute it numerically to high precision.
# In matter era, everything is determined by two parameters: δ₀ and Ω_b.
# Given the constraint δ₀(2+δ₀)/3 = Ω_DM:

# x = δ₀/Ω_b, solve quadratic for δ_m:
# (x/3)δ² + (1+2x/3)δ - x = 0
# δ_m = (-1-2x/3 + sqrt((1+2x/3)²+4x²/3))/(2x/3)

x = delta_0/Omega_b
delta_m_analytic = (-(1+2*x/3) + np.sqrt((1+2*x/3)**2 + 4*x**2/3))/(2*x/3)
print(f"  δ_m = {delta_m_analytic:.8f}")

Omega_K_m_analytic = delta_m_analytic*(2+delta_m_analytic)/3
ratio_m_analytic = Omega_b * Omega_K_m_analytic / ((1-Omega_K_m_analytic)*Omega_DM)
print(f"  Ω_K = {Omega_K_m_analytic:.8f}")
print(f"  ρ_K/ρ_CDM = {ratio_m_analytic:.8f}")

# ============================================================
# HOW ABOUT NEAR z_eq? Mixed matter+radiation
# ============================================================
print("\n" + "=" * 70)
print("DETAILED z-DEPENDENCE NEAR MATTER-RADIATION EQUALITY")
print("=" * 70)

# z_eq for self-consistent case is different from ΛCDM!
# Self-consistent: radiation = Ω_K + baryons at z_eq
# In self-consistent: Ω_r(1+z)⁴ = Ω_b(1+z)³ + ρ_K
# But ρ_K depends on H which depends on ρ_K... complicated.

# Let's just find z_eq numerically for both cases.
# ΛCDM: z_eq where Ω_r(1+z) = Ω_b + Ω_DM
z_eq_LCDM = (Omega_b + Omega_DM)/Omega_r - 1
print(f"  z_eq(ΛCDM) = (Ω_b+Ω_DM)/Ω_r - 1 = {z_eq_LCDM:.0f}")

# Self-consistent: z_eq where radiation fraction = matter fraction
# Ω_r_frac(z) = Ω_r(1+z)⁴/H²_SC(z)
# Ω_m_frac(z) = Ω_b(1+z)³/H²_SC(z) + Ω_K(z)
# Equal when: Ω_r(1+z)⁴ = Ω_b(1+z)³ + H²_SC·Ω_K
# = Ω_b(1+z)³ + H²_std·Ω_K/(1-Ω_K)

def matter_minus_rad(z):
    d = delta_SC(z)
    h2 = H2_SC(z)
    omega_K = d*(2+d)/3
    matter = Omega_b*(1+z)**3/h2 + omega_K
    rad = Omega_r*(1+z)**4/h2
    return matter - rad

z_eq_SC = brentq(matter_minus_rad, 100, 10000)
print(f"  z_eq(self-consistent) = {z_eq_SC:.0f}")
print(f"  Ratio z_eq_SC/z_eq_ΛCDM = {z_eq_SC/z_eq_LCDM:.4f}")

# ============================================================
# Detailed ρ_K/ρ_CDM with understanding
# ============================================================
print("\n" + "=" * 70)
print("ρ_K/ρ_CDM WITH PHYSICAL EXPLANATION")
print("=" * 70)

z_list = [0, 0.5, 1, 2, 5, 10, 50, 100, 500, 1100, 2000, 3400, 5000, 10000, 50000]

print(f"\n{'z':>8s} | {'δ_SC':>10s} | {'δ²/(2δ)':>10s} | {'ρ_K/ρ_CDM':>10s} | {'H²_SC/H²_LCDM':>14s} | Note")
print("-" * 85)

for z in z_list:
    d = delta_SC(z)
    h2sc = H2_SC(z)
    h2lcdm = H2_LCDM(z)

    rho_K = h2sc * d*(2+d)/3
    rho_CDM = Omega_DM*(1+z)**3
    ratio = rho_K / rho_CDM

    quad_frac = d / (2)  # δ/(2) is roughly the quadratic correction
    h_ratio = h2sc / h2lcdm

    note = ""
    if z < 1:
        note = "Λ-dominated"
    elif z < 100:
        note = "transition"
    elif z < 3000:
        note = "matter era"
    elif z < 5000:
        note = "near z_eq"
    else:
        note = "radiation era"

    print(f"{z:>8g} | {d:>10.6f} | {quad_frac:>10.4f} | {ratio:>10.4f} | {h_ratio:>14.6f} | {note}")

# ============================================================
# THE KEY DIAGNOSTIC: decompose the ratio
# ============================================================
print("\n" + "=" * 70)
print("DECOMPOSITION: ρ_K/ρ_CDM = (H²_SC/H²_ΛCDM) × [δ(2+δ)/3] / [Ω_CDM/(Ω_b+Ω_DM+Ω_r(1+z)/...)]")
print("=" * 70)

# ρ_K = H²_SC·δ(2+δ)/(3H₀²) × (3H₀²/(8πG)) = (H²_SC/H₀²)·δ(2+δ)/3 × ρ_crit0
# ρ_CDM = Ω_DM(1+z)³ × ρ_crit0
# ratio = (H²_SC/H₀²)·δ(2+δ)/3 / [Ω_DM(1+z)³]
#       = (H²_SC/H²_ΛCDM) × (H²_ΛCDM/H₀²)·Ω_K / [Ω_DM(1+z)³]
#       = (H²_SC/H²_ΛCDM) × Ω_K/Ω_CDM_frac
# where Ω_CDM_frac = Ω_DM(1+z)³ / (H²_ΛCDM/H₀²)

print(f"\n{'z':>6s} | {'H²_SC/H²_ΛCDM':>14s} | {'Ω_K':>10s} | {'Ω_CDM_frac':>11s} | {'Ω_K/Ω_CDM':>10s} | {'product':>10s}")
print("-" * 75)

for z in [0, 10, 100, 1100, 3400, 10000]:
    d = delta_SC(z)
    h2sc = H2_SC(z)
    h2lcdm = H2_LCDM(z)

    Omega_K = d*(2+d)/3
    Omega_CDM_frac = Omega_DM*(1+z)**3 / h2lcdm
    h_fac = h2sc / h2lcdm
    omega_ratio = Omega_K / Omega_CDM_frac
    product = h_fac * omega_ratio

    print(f"{z:>6d} | {h_fac:>14.6f} | {Omega_K:>10.6f} | {Omega_CDM_frac:>11.6f} | {omega_ratio:>10.6f} | {product:>10.6f}")

# ============================================================
# Can we reduce the excess by modifying the Khronon action?
# ============================================================
print("\n" + "=" * 70)
print("WHAT WOULD MAKE ρ_K/ρ_CDM = 1 EXACTLY?")
print("=" * 70)

print("""
For ρ_K/ρ_CDM = 1 exactly, we need ρ_K = ρ_CDM,0·(1+z)³.
This requires ρ_K ∝ (1+z)³ exactly, i.e., w=0 exactly.

Current: ρ_K = H²δ(2+δ)/(8πG).
If δ is small: ρ_K ≈ H²·2δ/(8πG) = μ²δ·c⁴/(4πG).
Conservation: μ²δ = const/a³ → H²δ = const·(1+z)³.
So ρ_K ≈ const·(1+z)³·c⁴/(4πG) → w ≈ 0, but the δ² term gives w ≠ 0.

The δ² contribution: ρ_K = H²(2δ+δ²)/(8πG) = const·(1+z)³/a³ · (1+δ/2)...
No, this isn't right. Let me be more careful.

H²δ = C(1+z)³ (exact conservation)
ρ_K = H²δ(2+δ)/(8πG) = C(1+z)³(2+δ)/(8πG)

So ρ_K = C(1+z)³·(2+δ(z))/(8πG).
For w=0 exactly, need 2+δ(z) = const. But δ varies!

The w ≠ 0 comes from the (2+δ) factor varying with z.
w_eff = -1/3 · d ln(2+δ)/d ln(1+z) ≈ -(1/3)·(dδ/dlna)/(2+δ)
""")

# Compute d(2+δ)/d(ln a) to show the effect
z_fine = np.logspace(-1, 4, 3000)
delta_fine = np.array([delta_SC(z) for z in z_fine])
lna = -np.log(1+z_fine)
factor = 2 + delta_fine
d_factor_dlna = np.gradient(np.log(factor), lna)
w_from_factor = -d_factor_dlna / 3

print(f"\n  w contribution from (2+δ) variation:")
for zt in [0, 10, 100, 1100, 10000]:
    idx = np.argmin(np.abs(z_fine - zt))
    print(f"    z={zt:>5d}: δ={delta_fine[idx]:.4f}, d ln(2+δ)/d ln a = {d_factor_dlna[idx]:+.6f}, w_contrib = {w_from_factor[idx]:+.6f}")

# ============================================================
# Also compute effective w from full ρ_K
# ============================================================
rho_K_fine = np.array([H2_SC(z) for z in z_fine]) * delta_fine * (2+delta_fine) / 3
ln_rho_K = np.log(rho_K_fine)
dlnrho_dlna = np.gradient(ln_rho_K, lna)
w_eff_fine = -1 - dlnrho_dlna / 3

print(f"\n  Full w_eff(z):")
for zt in [0, 1, 10, 100, 500, 1100, 3400, 10000]:
    idx = np.argmin(np.abs(z_fine - zt))
    print(f"    z={zt:>5d}: w_eff = {w_eff_fine[idx]:+.8f}")

# ============================================================
# KEY QUESTION: What's the MAXIMUM acceptable excess?
# ============================================================
print("\n" + "=" * 70)
print("CMB IMPACT ASSESSMENT")
print("=" * 70)

d_cmb = delta_SC(1100)
ratio_cmb = H2_SC(1100) * d_cmb*(2+d_cmb)/3 / (Omega_DM*(1+1100)**3)
omega_K_cmb = d_cmb*(2+d_cmb)/3
omega_CDM_cmb = Omega_DM*(1+1100)**3 / H2_LCDM(1100)

print(f"\n  At z = 1100:")
print(f"    ρ_K/ρ_CDM = {ratio_cmb:.4f} → {(ratio_cmb-1)*100:.2f}% excess")
print(f"    Ω_K = {omega_K_cmb:.4f} vs Ω_CDM(ΛCDM) = {omega_CDM_cmb:.4f}")
print(f"    H²_SC/H²_ΛCDM = {H2_SC(1100)/H2_LCDM(1100):.4f}")

# What matters for CMB:
# 1. Sound horizon r_s: depends on ∫ c_s/H dz, so H(z) matters
# 2. Silk damping: depends on H
# 3. Peak heights: depend on matter/radiation ratio at decoupling

# The 17% excess in ρ_K/ρ_CDM translates to:
# - Different H at recombination → different angular diameter distance
# - Different z_eq → different peak ratios

# Fractional change in H at z=1100:
frac_H_change = np.sqrt(H2_SC(1100)/H2_LCDM(1100)) - 1
print(f"\n    ΔH/H at z=1100: {frac_H_change*100:.2f}%")
print(f"    This shifts θ_s by ~{frac_H_change*100:.1f}%")
print(f"    Planck precision on θ_s: ~0.03%")
print(f"    → {abs(frac_H_change)*100/0.03:.0f}σ tension")

# z_eq difference:
print(f"\n    z_eq(SC) = {z_eq_SC:.0f} vs z_eq(ΛCDM) = {z_eq_LCDM:.0f}")
print(f"    Difference: {abs(z_eq_SC-z_eq_LCDM)/z_eq_LCDM*100:.1f}%")
print(f"    This affects the 1st-to-3rd peak ratio significantly")

# ============================================================
# WHAT IF: we allow δ₀ to be a free parameter?
# ============================================================
print("\n" + "=" * 70)
print("PARAMETER SCAN: Vary δ₀ to minimize CMB tension")
print("=" * 70)

# If we allow Ω_DM to vary, what value gives ρ_K/ρ_CDM(z=1100) closest to 1?
# This means adjusting δ₀ so that the effective DM density at z=1100 matches.

def compute_ratio_at_z(Omega_DM_try, z_target=1100):
    """Compute ρ_K/ρ_CDM at z_target for given Ω_DM."""
    d0 = -1 + np.sqrt(1 + 3*Omega_DM_try)
    OmL = 1 - Omega_b - Omega_r - Omega_DM_try

    def h2std(z):
        return Omega_b*(1+z)**3 + Omega_r*(1+z)**4 + OmL

    x0 = d0 / (1 - Omega_DM_try)
    h2std0 = h2std(0)

    def x_z(z):
        return x0 * h2std0 / h2std(z) * (1+z)**3

    def delta_z(z):
        return solve_delta(x_z(z))

    def h2sc(z):
        d = delta_z(z)
        return h2std(z) / (1 - d*(2+d)/3)

    d = delta_z(z_target)
    rK = h2sc(z_target) * d*(2+d)/3
    rCDM = Omega_DM_try * (1+z_target)**3
    return rK/rCDM, d, h2sc(z_target)

print(f"\n{'Ω_DM':>8s} | {'δ₀':>10s} | {'ρ_K/ρ_CDM(1100)':>16s} | {'excess':>10s}")
print("-" * 55)

for om_try in [0.10, 0.15, 0.20, 0.25, 0.265, 0.30, 0.35, 0.40]:
    d0_try = -1 + np.sqrt(1 + 3*om_try)
    ratio_try, d_cmb_try, _ = compute_ratio_at_z(om_try, 1100)
    print(f"{om_try:>8.3f} | {d0_try:>10.6f} | {ratio_try:>16.6f} | {(ratio_try-1)*100:>+10.2f}%")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY OF SELF-CONSISTENT CALCULATION")
print("=" * 70)

print(f"""
1. SELF-CONSISTENT SOLUTION:
   - The cubic δ[1+x(2+δ)/3] = x reduces to a quadratic (since δ appears
     at most quadratically in the Friedmann equation).
   - At z=0: δ₀ = {delta_0:.6f}, Ω_DM = {Omega_DM:.3f} (by construction)

2. ρ_K/ρ_CDM RATIO:
   - At z=0: 1.0000 (by construction)
   - At z=1100: {ratio_cmb:.4f} ({(ratio_cmb-1)*100:.2f}% excess)
   - At z=3400: ~1.08 (8% excess)
   - At z→∞: ratio < 1 (Khronon is subdominant in radiation era)

3. ROOT CAUSE:
   - ρ_K = H²δ(2+δ)/(8πG) has a QUADRATIC δ dependence
   - Conservation law H²δ = const·(1+z)³ guarantees 2δ part scales as (1+z)³
   - But the δ² part gives an EXTRA contribution that varies with z
   - At z=0: δ₀ ≈ 0.34, so δ² correction is ~17%
   - In deep matter era: δ ≈ 0.90, so δ² correction is ~45%
   - Net effect: ρ_K grows faster than (1+z)³, giving w_eff < 0 slightly

4. IMPACT ON CMB:
   - H(z=1100) is {frac_H_change*100:+.2f}% different from ΛCDM
   - z_eq shifts from {z_eq_LCDM:.0f} to {z_eq_SC:.0f}
   - Both effects are detectable by Planck at high significance
   - Self-consistency REDUCED the naive excess from 20.5% to 17.1%
     but this is NOT sufficient to resolve the CMB tension

5. CRITICAL INSIGHT:
   The excess comes from the NONLINEAR (δ²) term in ρ_K = H²δ(2+δ)/(8πG).
   This is intrinsic to the Khronon kinetic structure.

   Three possible resolutions:
   (a) The Khronon kinetic function K(Y) could have higher-order terms
       that modify the δ² coefficient at high δ
   (b) A perturbation-level treatment (δ << 1) would eliminate this issue
       but then Ω_DM = 2δ/3 requires δ₀ ≈ 0.40, still not small
   (c) The Boltzmann hierarchy treatment may redistribute the excess
       between background and perturbation effects
""")

# ============================================================
# PLOTS v2
# ============================================================
outdir = Path("/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/cmb_plots")
outdir.mkdir(exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

z_plot = np.logspace(-1, 4.5, 2000)

# Compute arrays
delta_sc_arr = np.array([delta_SC(z) for z in z_plot])
delta_lcdm_arr = delta_0 * (1+z_plot)**3 / np.array([H2_LCDM(z) for z in z_plot])
h2_sc_arr = np.array([H2_SC(z) for z in z_plot])
h2_lcdm_arr = np.array([H2_LCDM(z) for z in z_plot])
rho_K_arr = h2_sc_arr * delta_sc_arr*(2+delta_sc_arr)/3
rho_CDM_arr = Omega_DM*(1+z_plot)**3
ratio_arr = rho_K_arr / rho_CDM_arr

# Plot 1: δ(z)
ax = axes[0, 0]
ax.loglog(1+z_plot, delta_sc_arr, 'b-', lw=2, label=r'$\delta_{\rm SC}$')
ax.loglog(1+z_plot, delta_lcdm_arr, 'r--', lw=2, label=r'$\delta_{\Lambda CDM}$')
ax.axvline(1101, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('1+z'); ax.set_ylabel(r'$\delta(z)$')
ax.set_title(r'$\delta(z)$: self-consistent vs naive'); ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: ρ_K/ρ_CDM
ax = axes[0, 1]
ax.semilogx(1+z_plot, ratio_arr, 'b-', lw=2)
ax.axhline(1.0, color='k', ls='--', alpha=0.5)
ax.fill_between(1+z_plot, 0.95, 1.05, alpha=0.1, color='green', label='5% band')
ax.axvline(1101, color='gray', ls=':', alpha=0.5, label='z=1100')
ax.axvline(3401, color='orange', ls=':', alpha=0.5, label='z=3400')
ax.set_xlabel('1+z'); ax.set_ylabel(r'$\rho_K/\rho_{\rm CDM}$')
ax.set_title(r'Energy density ratio'); ax.legend()
ax.set_ylim(0.7, 1.4); ax.grid(True, alpha=0.3)

# Plot 3: w_eff
ax = axes[0, 2]
z_w = z_fine[(z_fine > 0.5) & (z_fine < 2e4)]
w_w = w_eff_fine[(z_fine > 0.5) & (z_fine < 2e4)]
ax.semilogx(1+z_w, w_w, 'b-', lw=2)
ax.axhline(0, color='k', ls='--', alpha=0.5)
ax.axvline(1101, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('1+z'); ax.set_ylabel(r'$w_{\rm eff}$')
ax.set_title(r'Effective equation of state'); ax.set_ylim(-0.05, 0.15)
ax.grid(True, alpha=0.3)

# Plot 4: H_SC/H_LCDM
ax = axes[1, 0]
ax.semilogx(1+z_plot, np.sqrt(h2_sc_arr/h2_lcdm_arr), 'b-', lw=2)
ax.axhline(1.0, color='k', ls='--', alpha=0.5)
ax.axvline(1101, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('1+z'); ax.set_ylabel(r'$H_{\rm SC}/H_{\Lambda CDM}$')
ax.set_title('Hubble rate ratio'); ax.grid(True, alpha=0.3)

# Plot 5: Ω fractions
ax = axes[1, 1]
Om_K_arr = delta_sc_arr*(2+delta_sc_arr)/3
Om_b_arr = Omega_b*(1+z_plot)**3 / h2_sc_arr
Om_r_arr = Omega_r*(1+z_plot)**4 / h2_sc_arr
Om_L_arr = Omega_Lambda / h2_sc_arr

ax.fill_between(1+z_plot, 0, Om_r_arr, alpha=0.3, color='red', label=r'$\Omega_r$')
ax.fill_between(1+z_plot, Om_r_arr, Om_r_arr+Om_b_arr, alpha=0.3, color='blue', label=r'$\Omega_b$')
ax.fill_between(1+z_plot, Om_r_arr+Om_b_arr, Om_r_arr+Om_b_arr+Om_K_arr, alpha=0.3, color='purple', label=r'$\Omega_K$ (Khronon)')
ax.fill_between(1+z_plot, Om_r_arr+Om_b_arr+Om_K_arr, 1, alpha=0.3, color='green', label=r'$\Omega_\Lambda$')
ax.set_xscale('log')
ax.axvline(1101, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('1+z'); ax.set_ylabel(r'$\Omega$')
ax.set_title('Energy density fractions (SC)'); ax.legend(fontsize=8)
ax.set_ylim(0, 1.05); ax.grid(True, alpha=0.3)

# Plot 6: δ_SC(z) × (2+δ) / 3 vs CDM fraction
ax = axes[1, 2]
Om_CDM_lcdm = Omega_DM*(1+z_plot)**3 / h2_lcdm_arr
ax.semilogx(1+z_plot, Om_K_arr, 'b-', lw=2, label=r'$\Omega_K$ (SC)')
ax.semilogx(1+z_plot, Om_CDM_lcdm, 'r--', lw=2, label=r'$\Omega_{\rm CDM}$ ($\Lambda$CDM)')
ax.axvline(1101, color='gray', ls=':', alpha=0.5)
ax.axvline(3401, color='orange', ls=':', alpha=0.5)
ax.set_xlabel('1+z'); ax.set_ylabel(r'$\Omega(z)$')
ax.set_title(r'$\Omega_K$ vs $\Omega_{\rm CDM}$'); ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(outdir / 'self_consistent_delta_v2.png', dpi=150)
print(f"\nPlots saved to {outdir}/self_consistent_delta_v2.png")

# Also save a focused CMB-era plot
fig2, ax2 = plt.subplots(1, 1, figsize=(8, 5))
z_cmb_range = np.logspace(2, 4, 500)
ratio_cmb_range = []
for z in z_cmb_range:
    d = delta_SC(z)
    h2 = H2_SC(z)
    rK = h2 * d*(2+d)/3
    rC = Omega_DM*(1+z)**3
    ratio_cmb_range.append(rK/rC)
ratio_cmb_range = np.array(ratio_cmb_range)

ax2.semilogx(1+z_cmb_range, ratio_cmb_range, 'b-', lw=2.5)
ax2.axhline(1.0, color='k', ls='--', lw=1)
ax2.axhline(1.05, color='green', ls='--', lw=1, alpha=0.5, label='5% threshold')
ax2.axhline(0.95, color='green', ls='--', lw=1, alpha=0.5)
ax2.axvline(1101, color='red', ls=':', lw=2, label='z=1100 (recombination)')
ax2.axvline(3401, color='orange', ls=':', lw=1.5, label=r'z=3400 ($z_{eq}^{\Lambda CDM}$)')

# Mark the actual SC z_eq
ax2.axvline(1+z_eq_SC, color='purple', ls=':', lw=1.5, label=f'z={z_eq_SC:.0f} ($z_{{eq}}^{{SC}}$)')

ax2.set_xlabel('1+z', fontsize=13)
ax2.set_ylabel(r'$\rho_K / \rho_{\rm CDM}$', fontsize=13)
ax2.set_title(r'Khronon vs CDM energy density ratio (self-consistent)', fontsize=14)
ax2.legend(fontsize=10)
ax2.set_ylim(0.9, 1.35)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(outdir / 'self_consistent_ratio_cmb.png', dpi=150)
print(f"Focused CMB plot saved to {outdir}/self_consistent_ratio_cmb.png")

print("\nAll computations complete.")
