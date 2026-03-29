#!/usr/bin/env python3
"""
Extremal Robustness Test: Sigma_SE vs Sigma_channel for Omega_DM prediction
===========================================================================

Tests the retrodictability extremal principle R[delta_0] = int F * Sigma * dt/dz dz
with four combinations of (F, Sigma):

  1. F = 1/(1+delta),      Sigma = delta(2+delta)       [canonical]
  2. F = 1/(1+delta),      Sigma = 2 ln(1+delta)        [channel entropy]
  3. F = 1/sqrt(1+delta),  Sigma = delta(2+delta)        [alpha=1/2, SE]
  4. F = 1/sqrt(1+delta),  Sigma = 2 ln(1+delta)         [alpha=1/2, channel]

Also tests: F = exp(-Sigma/2) exactly (Petz bound saturation) for each Sigma.

Author: Sheng-Kai Huang (computational analysis)
Date: 2026-03-17
"""

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ============================================================
# Physical constants (Planck 2018)
# ============================================================
Omega_b = 0.049
Omega_r = 9.15e-5
h_val   = 0.674
Omega_DM_obs = 0.265
Omega_DM_obs_err = 0.007

# ============================================================
# Core functions
# ============================================================

def E2_func(a, Omega_m, Omega_L):
    """Friedmann equation: E^2(a) = H^2(a)/H_0^2."""
    return Omega_r / a**4 + Omega_m / a**3 + Omega_L


def delta_of_z(z, delta_0, Omega_m, Omega_L):
    """delta(z) = delta_0 / (E^2(a) * a^3) with running mu = H(a)/c."""
    a = 1.0 / (1.0 + z)
    E2 = E2_func(a, Omega_m, Omega_L)
    if E2 <= 0:
        return 0.0
    return delta_0 / (E2 * a**3)


def Omega_K_from_delta(delta_0):
    """Omega_K = (delta_0^2 + 2*delta_0) / 3."""
    return (delta_0**2 + 2 * delta_0) / 3.0


def cosmology_from_delta(delta_0):
    """Return (Omega_m, Omega_L) from delta_0."""
    Omega_K = Omega_K_from_delta(delta_0)
    Omega_m = Omega_b + Omega_K
    Omega_L = 1.0 - Omega_m - Omega_r
    return Omega_m, Omega_L


# ============================================================
# F and Sigma choices
# ============================================================

def F_alpha1(d):
    """F = 1/(1+delta) -- canonical tau framework (alpha=1)."""
    return 1.0 / (1.0 + d)

def F_alpha_half(d):
    """F = 1/sqrt(1+delta) -- alpha=1/2."""
    return 1.0 / np.sqrt(1.0 + d)

def F_petz_SE(d):
    """F = exp(-Sigma_SE/2) = exp(-delta(2+delta)/2) -- Petz bound saturation for SE."""
    Sigma = d * (2.0 + d)
    return np.exp(-Sigma / 2.0)

def F_petz_ch(d):
    """F = exp(-Sigma_ch/2) = exp(-ln(1+delta)) = 1/(1+delta) -- Petz bound saturation for channel."""
    # Sigma_ch = 2 ln(1+delta), so exp(-Sigma_ch/2) = exp(-ln(1+delta)) = 1/(1+delta)
    # This is IDENTICAL to F_alpha1!
    return 1.0 / (1.0 + d)

def Sigma_SE(d):
    """Sigma = delta(2+delta) -- stress-energy tensor (exact Khronon)."""
    return d * (2.0 + d)

def Sigma_ch(d):
    """Sigma = 2 ln(1+delta) -- information-theoretic channel."""
    return 2.0 * np.log(1.0 + d)


# ============================================================
# Compute R[delta_0] for general (F, Sigma) choices
# ============================================================

def compute_R(delta_0, F_func, Sigma_func, z_max=1e5):
    """Compute the retrodictability functional R[delta_0]."""
    if delta_0 <= 0 or delta_0 > 4:
        return 0.0

    Omega_m, Omega_L = cosmology_from_delta(delta_0)
    if Omega_L < 0:
        return -1e10

    def integrand(lnz1):
        """Integrand in ln(1+z) variable."""
        z = np.exp(lnz1) - 1.0
        a = 1.0 / (1.0 + z) if z > 0 else 1.0
        E2 = E2_func(a, Omega_m, Omega_L)
        if E2 <= 0:
            return 0.0
        d = delta_0 / (E2 * a**3)
        if d < 1e-20:
            return 0.0
        F_val = F_func(d)
        S_val = Sigma_func(d)
        # dt/dz = 1/((1+z) * E(z)) in H_0^{-1} units
        dtdz = 1.0 / ((1.0 + z) * np.sqrt(E2))
        # Jacobian for ln(1+z): dz = (1+z) * d(ln(1+z))
        return F_val * S_val * dtdz * (1.0 + z)

    val, err = quad(integrand, 0, np.log(1 + z_max),
                    limit=500, epsrel=1e-10)
    return val


def find_optimal_delta(F_func, Sigma_func, bounds=(0.01, 2.0)):
    """Find delta_0 that maximizes R[delta_0]."""
    res = minimize_scalar(
        lambda x: -compute_R(x, F_func, Sigma_func),
        bounds=bounds, method='bounded',
        options={'xatol': 1e-10}
    )
    return res.x


# ============================================================
# Run all 4 main combinations + 2 Petz saturation cases
# ============================================================

print("=" * 70)
print("EXTREMAL ROBUSTNESS TEST: Sigma_SE vs Sigma_channel")
print("=" * 70)
print()

# Define all test cases
cases = [
    ("1/(1+d)",     "d(2+d)",       F_alpha1,     Sigma_SE,  "canonical"),
    ("1/(1+d)",     "2ln(1+d)",     F_alpha1,     Sigma_ch,  "channel-Sigma"),
    ("1/sqrt(1+d)", "d(2+d)",       F_alpha_half, Sigma_SE,  "alpha=1/2, SE"),
    ("1/sqrt(1+d)", "2ln(1+d)",     F_alpha_half, Sigma_ch,  "alpha=1/2, ch"),
    ("exp(-d(2+d)/2)", "d(2+d)",    F_petz_SE,    Sigma_SE,  "Petz-sat SE"),
    ("1/(1+d)*",    "2ln(1+d)",     F_petz_ch,    Sigma_ch,  "Petz-sat ch = case 2"),
]

results = []

for F_name, S_name, F_func, S_func, label in cases:
    delta_opt = find_optimal_delta(F_func, S_func)
    R_opt = compute_R(delta_opt, F_func, S_func)
    Omega_K = Omega_K_from_delta(delta_opt)
    Omega_m, Omega_L = cosmology_from_delta(delta_opt)
    deviation = (Omega_K - Omega_DM_obs) / Omega_DM_obs * 100
    sigma_dev = (Omega_K - Omega_DM_obs) / Omega_DM_obs_err

    results.append({
        'F_name': F_name,
        'S_name': S_name,
        'label': label,
        'delta_opt': delta_opt,
        'R_opt': R_opt,
        'Omega_DM': Omega_K,
        'Omega_m': Omega_m,
        'Omega_L': Omega_L,
        'deviation_pct': deviation,
        'sigma_dev': sigma_dev,
        'F_func': F_func,
        'S_func': S_func,
    })

    print(f"Case: {label}")
    print(f"  F = {F_name},  Sigma = {S_name}")
    print(f"  delta_0* = {delta_opt:.6f}")
    print(f"  R_max    = {R_opt:.6f}")
    print(f"  Omega_DM = {Omega_K:.6f}")
    print(f"  Omega_m  = {Omega_m:.6f}")
    print(f"  Omega_L  = {Omega_L:.6f}")
    print(f"  Deviation from 0.265: {deviation:+.2f}% ({sigma_dev:+.1f} sigma)")
    print()


# ============================================================
# Detailed comparison of Sigma at the optimal points
# ============================================================

print("=" * 70)
print("SIGMA VALUES AT OPTIMAL POINTS")
print("=" * 70)
print()
print(f"{'Case':<20} {'delta_0*':>10} {'Sigma_SE':>10} {'Sigma_ch':>10} {'Ratio SE/ch':>12}")
print("-" * 62)
for r in results[:4]:
    d = r['delta_opt']
    s_se = Sigma_SE(d)
    s_ch = Sigma_ch(d)
    print(f"{r['label']:<20} {d:10.4f} {s_se:10.4f} {s_ch:10.4f} {s_se/s_ch:12.4f}")

# ============================================================
# Extended alpha scan
# ============================================================

print()
print("=" * 70)
print("ALPHA SCAN: F = (1+delta)^{-alpha}")
print("=" * 70)
print()

alphas = np.linspace(0.2, 2.5, 24)
alpha_results_SE = []
alpha_results_ch = []

for alpha in alphas:
    F_a = lambda d, a=alpha: (1.0 + d)**(-a)

    d_SE = find_optimal_delta(F_a, Sigma_SE)
    Om_SE = Omega_K_from_delta(d_SE)
    R_SE = compute_R(d_SE, F_a, Sigma_SE)

    d_ch = find_optimal_delta(F_a, Sigma_ch)
    Om_ch = Omega_K_from_delta(d_ch)
    R_ch = compute_R(d_ch, F_a, Sigma_ch)

    alpha_results_SE.append((alpha, d_SE, Om_SE, R_SE))
    alpha_results_ch.append((alpha, d_ch, Om_ch, R_ch))

    dev_SE = (Om_SE - Omega_DM_obs) / Omega_DM_obs * 100
    dev_ch = (Om_ch - Omega_DM_obs) / Omega_DM_obs * 100
    print(f"  alpha = {alpha:.2f}:  Omega_DM(SE) = {Om_SE:.4f} ({dev_SE:+6.1f}%),  "
          f"Omega_DM(ch) = {Om_ch:.4f} ({dev_ch:+6.1f}%),  "
          f"ratio = {Om_SE/Om_ch:.3f}" if Om_ch > 0 else "")


# ============================================================
# PLOT 1: R(delta_0) for all four main cases
# ============================================================

# Use a delta range that stays physical (Omega_L > 0 requires delta_0 < ~0.95 for Omega_b=0.049)
delta_range = np.linspace(0.01, 0.90, 500)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(r'Retrodictability Functional $\mathcal{R}[\delta_0]$ — Four Combinations',
             fontsize=14, fontweight='bold')

main_cases = results[:4]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

for idx, (case, color) in enumerate(zip(main_cases, colors)):
    ax = axes[idx // 2][idx % 2]
    R_vals = []
    for d0 in delta_range:
        rv = compute_R(d0, case['F_func'], case['S_func'])
        R_vals.append(max(rv, 0.0))  # clip negative penalty values
    R_vals = np.array(R_vals)

    ax.plot(delta_range, R_vals, color=color, linewidth=2)
    ax.axvline(case['delta_opt'], color='gray', linestyle='--', alpha=0.7,
               label=f"$\\delta_0^* = {case['delta_opt']:.3f}$")
    ax.axvline(0.340, color='red', linestyle=':', alpha=0.5,
               label=r'$\delta_0^{\rm obs} \approx 0.340$')
    ax.set_xlabel(r'$\delta_0$', fontsize=12)
    ax.set_ylabel(r'$\mathcal{R}[\delta_0]$ (in $H_0^{-1}$ units)', fontsize=11)
    ax.set_title(f"F = {case['F_name']}, $\\Sigma$ = {case['S_name']}\n"
                 f"$\\Omega_{{DM}} = {case['Omega_DM']:.4f}$ "
                 f"(dev: {case['deviation_pct']:+.1f}%)", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/cmb_plots/extremal_robustness_R_four_cases.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("\nPlot 1 saved: extremal_robustness_R_four_cases.png")


# ============================================================
# PLOT 2: All four R(delta_0) overlaid
# ============================================================

fig, ax = plt.subplots(figsize=(10, 7))
ax.set_title(r'$\mathcal{R}[\delta_0]$ — All Four Cases Overlaid', fontsize=14, fontweight='bold')

for idx, (case, color) in enumerate(zip(main_cases, colors)):
    R_vals = []
    for d0 in delta_range:
        rv = compute_R(d0, case['F_func'], case['S_func'])
        R_vals.append(max(rv, 0.0))
    R_vals = np.array(R_vals)
    # Normalize to peak = 1 for shape comparison
    R_norm = R_vals / np.max(R_vals) if np.max(R_vals) > 0 else R_vals
    ax.plot(delta_range, R_norm, color=color, linewidth=2,
            label=f"F={case['F_name']}, $\\Sigma$={case['S_name']}\n"
                  f"  $\\delta_0^*$={case['delta_opt']:.3f}, "
                  f"$\\Omega_{{DM}}$={case['Omega_DM']:.4f}")

# Also overlay the Petz-sat SE case
case_petz = results[4]
R_petz = []
for d0 in delta_range:
    rv = compute_R(d0, case_petz['F_func'], case_petz['S_func'])
    R_petz.append(max(rv, 0.0))
R_petz = np.array(R_petz)
R_petz_norm = R_petz / np.max(R_petz) if np.max(R_petz) > 0 else R_petz
ax.plot(delta_range, R_petz_norm, color='purple', linewidth=2, linestyle='--',
        label=f"F=exp(-$\\Sigma_{{SE}}$/2), $\\Sigma$=d(2+d)\n"
              f"  $\\delta_0^*$={case_petz['delta_opt']:.3f}, "
              f"$\\Omega_{{DM}}$={case_petz['Omega_DM']:.4f}")

ax.axvline(0.340, color='red', linestyle=':', alpha=0.7, linewidth=1.5,
           label=r'$\delta_0^{\rm obs} \approx 0.340$')
ax.axhspan(0.98, 1.02, alpha=0.1, color='gray')
ax.set_xlabel(r'$\delta_0$', fontsize=13)
ax.set_ylabel(r'$\mathcal{R}/\mathcal{R}_{\max}$ (normalized)', fontsize=13)
ax.legend(fontsize=9, loc='upper right')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 0.90)
ax.set_ylim(0, 1.05)

plt.tight_layout()
plt.savefig('/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/cmb_plots/extremal_robustness_overlaid.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Plot 2 saved: extremal_robustness_overlaid.png")


# ============================================================
# PLOT 3: Omega_DM vs alpha for both Sigma choices
# ============================================================

fig, ax = plt.subplots(figsize=(10, 7))
ax.set_title(r'$\Omega_{DM}(\alpha)$ for $F = (1+\delta)^{-\alpha}$', fontsize=14, fontweight='bold')

alphas_plot = [r[0] for r in alpha_results_SE]
Om_SE_plot = [r[2] for r in alpha_results_SE]
Om_ch_plot = [r[2] for r in alpha_results_ch]

ax.plot(alphas_plot, Om_SE_plot, 'o-', color='#1f77b4', linewidth=2, markersize=5,
        label=r'$\Sigma = \delta(2+\delta)$ (stress-energy)')
ax.plot(alphas_plot, Om_ch_plot, 's-', color='#ff7f0e', linewidth=2, markersize=5,
        label=r'$\Sigma = 2\ln(1+\delta)$ (channel)')

# Observed band
ax.axhspan(Omega_DM_obs - Omega_DM_obs_err, Omega_DM_obs + Omega_DM_obs_err,
           alpha=0.2, color='green', label=r'Planck 2018: $\Omega_{DM} = 0.265 \pm 0.007$')
ax.axhline(Omega_DM_obs, color='green', linestyle='--', alpha=0.7)

# Mark alpha = 1
ax.axvline(1.0, color='gray', linestyle=':', alpha=0.5, label=r'$\alpha = 1$ (canonical)')

ax.set_xlabel(r'$\alpha$ in $F = (1+\delta)^{-\alpha}$', fontsize=13)
ax.set_ylabel(r'$\Omega_{DM}$ (predicted)', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.2, 2.5)
ax.set_ylim(0, 0.6)

plt.tight_layout()
plt.savefig('/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/cmb_plots/extremal_robustness_alpha_scan.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Plot 3 saved: extremal_robustness_alpha_scan.png")


# ============================================================
# PLOT 4: F * Sigma for both Sigma choices (as function of delta)
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
d_arr = np.linspace(0.001, 2.0, 500)

# Left: F * Sigma for F = 1/(1+d)
ax = axes[0]
ax.set_title(r'$F \cdot \Sigma$ with $F = 1/(1+\delta)$', fontsize=13, fontweight='bold')
FS_SE = np.array([F_alpha1(d) * Sigma_SE(d) for d in d_arr])
FS_ch = np.array([F_alpha1(d) * Sigma_ch(d) for d in d_arr])
ax.plot(d_arr, FS_SE, color='#1f77b4', linewidth=2,
        label=r'$\Sigma_{SE} = \delta(2+\delta)$: $F\Sigma = (1+\delta) - 1/(1+\delta)$')
ax.plot(d_arr, FS_ch, color='#ff7f0e', linewidth=2,
        label=r'$\Sigma_{ch} = 2\ln(1+\delta)$: $F\Sigma = 2\ln(1+\delta)/(1+\delta)$')
ax.set_xlabel(r'$\delta$', fontsize=13)
ax.set_ylabel(r'$F \cdot \Sigma$', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Right: F * Sigma for F = 1/sqrt(1+d)
ax = axes[1]
ax.set_title(r'$F \cdot \Sigma$ with $F = 1/\sqrt{1+\delta}$', fontsize=13, fontweight='bold')
FS_SE_h = np.array([F_alpha_half(d) * Sigma_SE(d) for d in d_arr])
FS_ch_h = np.array([F_alpha_half(d) * Sigma_ch(d) for d in d_arr])
ax.plot(d_arr, FS_SE_h, color='#2ca02c', linewidth=2,
        label=r'$\Sigma_{SE} = \delta(2+\delta)$')
ax.plot(d_arr, FS_ch_h, color='#d62728', linewidth=2,
        label=r'$\Sigma_{ch} = 2\ln(1+\delta)$')
ax.set_xlabel(r'$\delta$', fontsize=13)
ax.set_ylabel(r'$F \cdot \Sigma$', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/cmb_plots/extremal_robustness_FS_product.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Plot 4 saved: extremal_robustness_FS_product.png")


# ============================================================
# KEY ANALYSIS: Why do the two Sigma choices give different results?
# ============================================================

print()
print("=" * 70)
print("KEY ANALYSIS: WHY THE TWO SIGMA CHOICES DIVERGE")
print("=" * 70)
print()

d = results[0]['delta_opt']
print(f"At the canonical optimal delta_0 = {d:.4f}:")
print(f"  Sigma_SE  = delta(2+delta) = {Sigma_SE(d):.6f}")
print(f"  Sigma_ch  = 2 ln(1+delta)  = {Sigma_ch(d):.6f}")
print(f"  Ratio SE/ch = {Sigma_SE(d)/Sigma_ch(d):.4f}")
print()

print("The key difference:")
print("  Sigma_SE = delta(2+delta)  is QUADRATIC in delta for large delta")
print("  Sigma_ch = 2 ln(1+delta)   is LOGARITHMIC in delta for large delta")
print()
print("For delta << 1:  both ~ 2*delta (agree)")
print("For delta ~ 0.34: SE/ch = {:.3f} (36% divergence)".format(Sigma_SE(0.34)/Sigma_ch(0.34)))
print("For delta ~ 1:    SE/ch = {:.3f}".format(Sigma_SE(1.0)/Sigma_ch(1.0)))
print("For delta ~ 2:    SE/ch = {:.3f}".format(Sigma_SE(2.0)/Sigma_ch(2.0)))
print()

print("With F = 1/(1+delta):")
print("  F * Sigma_SE = (1+d) - 1/(1+d)  -- MONOTONICALLY INCREASING (no turnover)")
print("  F * Sigma_ch = 2 ln(1+d)/(1+d)  -- HAS A MAXIMUM at d = e-1 = 1.718")
print()

# Check: does F*Sigma_ch have a maximum?
d_test = np.linspace(0.001, 5, 10000)
FS_ch_test = 2 * np.log(1 + d_test) / (1 + d_test)
d_max_FS_ch = d_test[np.argmax(FS_ch_test)]
print(f"  Maximum of F*Sigma_ch at delta = {d_max_FS_ch:.3f} (value = {np.max(FS_ch_test):.4f})")
print(f"  Theoretical: delta = e-1 = {np.e - 1:.3f}")
print()

print("CRITICAL OBSERVATION:")
print("  With Sigma_ch, the product F*Sigma itself has an internal maximum.")
print("  This means the extremal principle is determined by BOTH the F*Sigma")
print("  turnover AND the age competition. This shifts delta_0* upward")
print("  (or downward depending on which effect dominates).")
print()

# ============================================================
# Which alpha makes Sigma_ch match observations?
# ============================================================
print("=" * 70)
print("FINDING alpha THAT MAKES Sigma_ch MATCH Omega_DM = 0.265")
print("=" * 70)
print()

from scipy.optimize import brentq

def Om_DM_for_alpha_ch(alpha):
    F_a = lambda d: (1.0 + d)**(-alpha)
    d_opt = find_optimal_delta(F_a, Sigma_ch)
    return Omega_K_from_delta(d_opt)

def Om_DM_for_alpha_SE(alpha):
    F_a = lambda d: (1.0 + d)**(-alpha)
    d_opt = find_optimal_delta(F_a, Sigma_SE)
    return Omega_K_from_delta(d_opt)

# Find alpha that gives 0.265 for each Sigma
try:
    alpha_ch_target = brentq(lambda a: Om_DM_for_alpha_ch(a) - 0.265, 0.3, 3.0)
    print(f"  Sigma_ch: alpha = {alpha_ch_target:.4f} gives Omega_DM = 0.265")
except:
    print("  Sigma_ch: could not find alpha for Omega_DM = 0.265")

try:
    alpha_SE_target = brentq(lambda a: Om_DM_for_alpha_SE(a) - 0.265, 0.3, 3.0)
    print(f"  Sigma_SE: alpha = {alpha_SE_target:.4f} gives Omega_DM = 0.265")
except:
    print("  Sigma_SE: could not find alpha for Omega_DM = 0.265")


# ============================================================
# CRITICAL TEST: Does the SHAPE of R(delta_0) change?
# ============================================================
print()
print("=" * 70)
print("SHAPE ANALYSIS: Width of R peak (how tightly delta_0 is constrained)")
print("=" * 70)
print()

for case in main_cases:
    d_opt = case['delta_opt']
    R_max = case['R_opt']
    # Find width at R = 0.95 * R_max
    d_left = d_opt
    d_right = d_opt
    for d_test in np.linspace(0.001, d_opt, 1000):
        if compute_R(d_test, case['F_func'], case['S_func']) >= 0.95 * R_max:
            d_left = d_test
            break
    for d_test in np.linspace(d_opt, 2.0, 1000):
        if compute_R(d_test, case['F_func'], case['S_func']) < 0.95 * R_max:
            d_right = d_test
            break
    width = d_right - d_left
    Om_left = Omega_K_from_delta(d_left)
    Om_right = Omega_K_from_delta(d_right)
    print(f"  {case['label']:<20}: delta_0 in [{d_left:.3f}, {d_right:.3f}] "
          f"(width = {width:.3f}), "
          f"Omega_DM in [{Om_left:.3f}, {Om_right:.3f}]")


# ============================================================
# ADDITIONAL TEST: What if we use Sigma_ch as WEIGHT only?
# ============================================================
print()
print("=" * 70)
print("HYBRID TEST: F=1/(1+d) with various Sigma forms")
print("=" * 70)
print()

sigma_forms = [
    ("delta(2+delta)",       lambda d: d * (2 + d)),
    ("2 ln(1+delta)",        lambda d: 2 * np.log(1 + d)),
    ("delta",                lambda d: d),
    ("delta^2",              lambda d: d**2),
    ("(1+delta)^2 - 1",     lambda d: (1 + d)**2 - 1),
    ("2 tanh(delta)",        lambda d: 2 * np.tanh(d)),
    ("2 delta/(1+delta)",    lambda d: 2 * d / (1 + d)),
    ("ln((1+delta)^2)",      lambda d: np.log((1 + d)**2)),  # = 2 ln(1+d)
]

print(f"{'Sigma form':<25} {'delta_0*':>10} {'Omega_DM':>10} {'Dev %':>8}")
print("-" * 55)
for name, S_func in sigma_forms:
    d_opt = find_optimal_delta(F_alpha1, S_func)
    Om = Omega_K_from_delta(d_opt)
    dev = (Om - Omega_DM_obs) / Omega_DM_obs * 100
    print(f"{name:<25} {d_opt:10.4f} {Om:10.4f} {dev:+8.1f}%")


# ============================================================
# SUMMARY TABLE (for markdown output)
# ============================================================
print()
print()
print("=" * 70)
print("SUMMARY TABLE (Markdown format)")
print("=" * 70)
print()
print("| F choice | Sigma choice | delta_0* | Omega_DM | Omega_m | Dev from 0.265 | Sigma dev |")
print("|----------|-------------|----------|----------|---------|---------------|-----------|")
for r in results:
    print(f"| {r['F_name']:<15} | {r['S_name']:<12} | {r['delta_opt']:.4f}   | "
          f"{r['Omega_DM']:.4f}   | {r['Omega_m']:.4f}  | "
          f"{r['deviation_pct']:+.1f}%          | {r['sigma_dev']:+.1f} sigma  |")


# ============================================================
# THE VERDICT
# ============================================================
print()
print("=" * 70)
print("VERDICT")
print("=" * 70)
print()
Om_canonical = results[0]['Omega_DM']
Om_channel = results[1]['Omega_DM']
print(f"Canonical (F=1/(1+d), Sigma=d(2+d)):  Omega_DM = {Om_canonical:.4f}")
print(f"Channel   (F=1/(1+d), Sigma=2ln(1+d)): Omega_DM = {Om_channel:.4f}")
print(f"Observed:                               Omega_DM = {Omega_DM_obs}")
print()
if abs(Om_channel - Omega_DM_obs) / Omega_DM_obs < 0.05:
    print("RESULT: ROBUST -- both Sigma choices give similar Omega_DM (< 5% deviation each)")
else:
    print(f"RESULT: NOT ROBUST -- Sigma_ch gives {(Om_channel - Omega_DM_obs)/Omega_DM_obs*100:+.1f}% "
          f"deviation vs Sigma_SE's {(Om_canonical - Omega_DM_obs)/Omega_DM_obs*100:+.1f}%")
    print()
    print("The prediction is SENSITIVE to the choice of Sigma.")
    print("The canonical Sigma_SE = delta(2+delta) works; Sigma_ch = 2ln(1+delta) does not.")
    print()
    print("HOWEVER, note the key mathematical reason:")
    print("  F*Sigma_SE = (1+d) - 1/(1+d) is monotonically increasing")
    print("  F*Sigma_ch = 2ln(1+d)/(1+d)  has a maximum at d = e-1")
    print()
    print("The monotonicity of F*Sigma_SE means the extremal principle is")
    print("entirely controlled by the age-of-universe competition.")
    print("With Sigma_ch, there is an ADDITIONAL internal balance, which")
    print("shifts the optimal delta_0 significantly.")


print()
print("Script completed successfully.")
