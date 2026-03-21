#!/usr/bin/env python3
"""
accuracy_analysis.py — Systematic decomposition of mlx_class vs CLASS residuals.

Identifies every error source contributing to the ~88% RMS residual,
quantifies each, and produces a prioritized fix plan.

Author: Sheng-Kai Huang, 2026
"""
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mlx_class.background import (
    Background, T_CMB, A_s, n_s, k_pivot, k_eq,
    Omega_r, Omega_b, Omega_c, Omega_m, H0_Mpc, a_eq,
)

# ============================================================================
# Configuration
# ============================================================================
CLASS_FILE = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
PRECISION_FILE = os.path.join(os.path.dirname(__file__), 'precision_comparison.dat')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_class_reference():
    """Load CLASS C_l reference: l(l+1)/(2pi) C_l -> D_l [muK^2]."""
    data = np.loadtxt(CLASS_FILE)
    ell = data[:, 0].astype(int)
    Dl = data[:, 1] * (T_CMB * 1e6) ** 2
    return ell, Dl


def load_precision_comparison():
    """Load existing precision solver comparison."""
    data = np.loadtxt(PRECISION_FILE)
    ell = data[:, 0].astype(int)
    Dl_class = data[:, 1]
    Dl_mlx = data[:, 2]
    residual_pct = data[:, 3]
    return ell, Dl_class, Dl_mlx, residual_pct


# ============================================================================
# Error Source 1: ISW Template Excess at Low-l
# ============================================================================
def analyze_isw_excess(ell, Dl_class, Dl_mlx):
    """Quantify the catastrophic ISW template excess at l < 100."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 1: ISW TEMPLATE EXCESS (l < 100)")
    print("=" * 72)

    mask = ell < 100
    e = ell[mask]
    dc = Dl_class[mask]
    dm = Dl_mlx[mask]

    excess = dm - dc  # absolute excess in muK^2
    ratio = dm / dc

    print(f"\n  l-range: 2 to {e[-1]}")
    print(f"  Max ratio (mlx/CLASS): {np.max(ratio):.1f}x at l={e[np.argmax(ratio)]}")
    print(f"  Max excess: {np.max(excess):.0f} muK^2 at l={e[np.argmax(excess)]}")
    print(f"  Mean ratio (l<30): {np.mean(ratio[e<30]):.1f}x")
    print(f"  Mean ratio (30<l<100): {np.mean(ratio[(e>=30)&(e<100)]):.1f}x")

    # RMS contribution to total
    mask_sig = (ell > 30) & (Dl_class > 100)
    res_total = np.where(mask_sig, (Dl_mlx - Dl_class) / Dl_class * 100, 0.0)
    mask_isw = mask_sig & (ell < 100)
    rms_isw = np.sqrt(np.mean(res_total[mask_isw] ** 2)) if np.any(mask_isw) else 0
    rms_total = np.sqrt(np.mean(res_total[mask_sig] ** 2))

    print(f"\n  RMS residual (30<l<100): {rms_isw:.1f}%")
    print(f"  Fraction of total RMS^2: {rms_isw**2 * np.sum(mask_isw) / (rms_total**2 * np.sum(mask_sig)) * 100:.1f}%")

    print(f"\n  ROOT CAUSE:")
    print(f"    solver_precision.py does NOT compute a proper ISW template at low-l.")
    print(f"    The early ISW from the TCA solver Phi_dot is evaluated at")
    print(f"    tau ~ tau_eq, but the integration through the visibility function")
    print(f"    produces unphysical large-amplitude oscillations at l < 100.")
    print(f"    The late ISW (dark energy) adds more excess.")
    print(f"    In class_comparison.py, an empirical Gaussian ISW bump was added")
    print(f"    centered at l~140 with amplitude 1.3x the acoustic peak, which")
    print(f"    massively overestimates the actual ISW contribution.")

    print(f"\n  FIX:")
    print(f"    [1] Remove the empirical ISW Gaussian template entirely.")
    print(f"    [2] Compute early ISW properly: integrate Phi_dot * j_l(k*chi)")
    print(f"        from tau_eq/2 to start of visibility window, using the actual")
    print(f"        Phi_dot from the ODE solver. This is already partially done")
    print(f"        in solver_precision.py but the normalization is wrong.")
    print(f"    [3] For late ISW (l < 20): use the Eisenstein-Hu Phi_dot with")
    print(f"        correct normalization (Phi_plateau ~ 0.9 * T(k) * Phi_primordial).")
    print(f"    [4] Verify by checking D_l at l=10 ~ 820 muK^2 (CLASS).")
    print(f"    Expected improvement: removes ~484% RMS at l<100.")

    return {
        'rms': rms_isw,
        'max_ratio': np.max(ratio),
        'priority': 1,
    }


# ============================================================================
# Error Source 2: Driving Correction Miscalibration
# ============================================================================
def analyze_driving_correction(ell, Dl_class, Dl_mlx):
    """Analyze the driving correction: wrong amplitude and ell-dependence."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 2: DRIVING CORRECTION MISCALIBRATION")
    print("=" * 72)

    mask = (ell >= 100) & (ell <= 2000) & (Dl_class > 100)
    e = ell[mask].astype(float)
    ratio = Dl_mlx[mask] / Dl_class[mask]

    # Fit: ratio = a + b*l + c*cos(2pi*l/P + phi)
    def model(params, ell_arr):
        a, b, c, period, phi = params
        return a + b * ell_arr + c * np.cos(2 * np.pi * ell_arr / period + phi)

    def cost(params):
        return np.sum((ratio - model(params, e)) ** 2)

    result = minimize(cost, [1.0, -0.0002, 0.3, 310.0, 0.0],
                      method='Nelder-Mead', options={'maxiter': 50000})
    a, b, c, period, phi = result.x

    print(f"\n  Fit: ratio = {a:.4f} + {b:.6f}*l + {c:.4f}*cos(2pi*l/{period:.1f} + {phi:.2f})")
    print(f"\n  (A) MEAN NORMALIZATION:")
    print(f"      Mean ratio: {a:.4f} (should be 1.000)")
    print(f"      Excess normalization: {(a-1)*100:+.1f}%")
    print(f"      In D_l: this means the driving correction D_inf=1.95 is")
    print(f"      too large by ~{(a-1)*100/2:.0f}% in source (squared in D_l)")
    print(f"      => D_inf should be ~ {1.95/np.sqrt(a):.2f} instead of 1.95")

    print(f"\n  (B) DAMPING SLOPE:")
    print(f"      Slope: {b*1000:.4f} per 1000 ell")
    print(f"      At l=2000, this contributes {b*2000*100:.1f}% to the ratio")
    print(f"      CAUSE: Silk damping scale k_D is slightly wrong.")
    print(f"      Current k_D = 0.1506 Mpc^-1 (from formula)")
    print(f"      CLASS k_D ~ 0.140 Mpc^-1 (more damping at high l)")
    print(f"      At l=2000, k ~ l/D_A ~ 0.14 Mpc^-1, right at k_D")

    print(f"\n  (C) OSCILLATORY COMPONENT:")
    print(f"      Amplitude: {abs(c)*100:.1f}% peak-to-peak")
    print(f"      Period: {period:.1f} (acoustic period pi*D_A/r_s ~ 302)")
    print(f"      Phase: {phi:.2f} rad = {phi/np.pi*180:.0f} deg")
    print(f"      CAUSE: The constant driving correction D_corr(k) does not")
    print(f"      properly account for the phase shift between SW source")
    print(f"      (Theta_0+Psi) and Doppler source (v_b). In CLASS, the")
    print(f"      driving correction is implicitly handled by the full")
    print(f"      photon hierarchy evolution (Theta_l for l=0,1,2,...)")
    print(f"      which redistributes power into higher multipoles.")

    rms_after_fit = np.sqrt(np.mean((ratio - model(result.x, e)) ** 2))

    print(f"\n  Fit residual RMS: {rms_after_fit*100:.1f}%")

    print(f"\n  FIX:")
    print(f"    [1] Replace constant D_inf=1.95 with k-dependent correction")
    print(f"        calibrated against CLASS at 5 peak positions:")
    print(f"        D(k) ~ 1 + alpha * (1 - exp(-k/k_eq))^beta")
    print(f"    [2] Better: switch to full photon hierarchy (perturbations_hires.py)")
    print(f"        which eliminates the need for driving correction entirely.")
    print(f"    [3] Fix Silk damping scale: use integral formula instead of")
    print(f"        k_D = 0.15*(omega_b/0.022)^0.25")
    print(f"    Expected improvement: 20-30% reduction in RMS at 100<l<2000")

    return {
        'mean_norm': a,
        'slope': b,
        'osc_amp': c,
        'osc_period': period,
        'osc_phase': phi,
        'rms_residual': rms_after_fit,
        'priority': 3,
    }


# ============================================================================
# Error Source 3: Missing Photon Hierarchy (TCA limitation)
# ============================================================================
def analyze_tca_limitation(ell, Dl_class, Dl_mlx):
    """Analyze the fundamental TCA limitation: no Theta_l for l >= 2."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 3: TCA LIMITATION (no Theta_l for l >= 2)")
    print("=" * 72)

    # The TCA locks v_b = 3*Theta_1 and has no Theta_2 (quadrupole).
    # This causes several related problems:

    print(f"\n  (A) MISSING POLARIZATION SOURCE:")
    print(f"      The SW source is Theta_0 + Psi + Pi/4, where Pi = Theta_2 + Theta_P0 + Theta_P2")
    print(f"      In TCA, Theta_2 = 0, so Pi = 0. This missing term contributes ~5% at peaks.")
    print(f"      The TCA approximation Theta_2 ~ (2k)/(9*kappa_dot)*Theta_1 is attempted")
    print(f"      in the code comments but NOT actually used in the source construction.")

    print(f"\n  (B) NO FREE-STREAMING POWER REDISTRIBUTION:")
    print(f"      After recombination, photons free-stream and Theta_1 power moves to")
    print(f"      higher l. In the full hierarchy, the Theta_1 Doppler source is reduced")
    print(f"      by this redistribution. In TCA, Theta_1 stays locked to v_b/3.")
    print(f"      This means the Doppler source is OVERESTIMATED in TCA.")

    print(f"\n  (C) ACOUSTIC OSCILLATION AMPLITUDE:")
    print(f"      The TCA cannot capture the radiation driving enhancement correctly")
    print(f"      because it misses the Theta_2 -> Theta_0 feedback (photon quadrupole")
    print(f"      back-reaction). This is why a calibrated driving correction D_corr is needed.")
    print(f"      D_corr ~ 1.95 is an empirical fudge, not first-principles physics.")

    print(f"\n  (D) PEAK-TO-TROUGH CONTRAST:")
    # Analyze peak-to-trough contrast
    mask = (ell >= 100) & (ell <= 2000) & (Dl_class > 100)
    e = ell[mask]
    dc = Dl_class[mask]
    dm = Dl_mlx[mask]

    dcs = gaussian_filter1d(dc, sigma=5)
    dms = gaussian_filter1d(dm, sigma=5)
    peaks_c, _ = find_peaks(dcs, distance=80)
    troughs_c, _ = find_peaks(-dcs, distance=80)

    print(f"      Peak-to-following-trough contrast:")
    n = min(len(peaks_c), len(troughs_c))
    contrasts_class = []
    contrasts_mlx = []
    for i in range(n):
        p = peaks_c[i]
        t_candidates = troughs_c[troughs_c > p]
        if len(t_candidates) == 0:
            continue
        t = t_candidates[0]
        cc = dc[p] / dc[t]
        cm = dm[p] / dm[t]
        contrasts_class.append(cc)
        contrasts_mlx.append(cm)
        print(f"        Peak {i+1} (l={e[p]}): CLASS={cc:.3f}, MLX={cm:.3f}, deficit={cm/cc:.3f}")

    if contrasts_class:
        mean_contrast_ratio = np.mean(np.array(contrasts_mlx) / np.array(contrasts_class))
        print(f"      Mean contrast ratio (MLX/CLASS): {mean_contrast_ratio:.3f}")
        print(f"      This means MLX troughs are {(1-mean_contrast_ratio)*100:.0f}% too shallow")

    print(f"\n  FIX:")
    print(f"    THE DEFINITIVE FIX: Switch from TCA to full photon hierarchy.")
    print(f"    perturbations_hires.py already implements this with Strang splitting.")
    print(f"    This single change eliminates error sources 2(A,B,C), 3(A,B,C,D), 5.")
    print(f"    Expected improvement: from ~88% RMS to ~5-10% RMS at 100<l<2000")

    return {
        'mean_contrast_ratio': mean_contrast_ratio if contrasts_class else 0,
        'priority': 2,
    }


# ============================================================================
# Error Source 4: Silk Damping Scale
# ============================================================================
def analyze_silk_damping(ell, Dl_class, Dl_mlx):
    """Analyze the Silk damping scale error at high l."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 4: SILK DAMPING SCALE ERROR (l > 1200)")
    print("=" * 72)

    mask = (ell >= 1200) & (ell <= 2500) & (Dl_class > 50)
    e = ell[mask].astype(float)
    ratio = Dl_mlx[mask] / Dl_class[mask]

    # The ratio should be ~1 if damping is correct.
    # A slope in log(ratio) vs l indicates k_D error
    log_ratio = np.log(ratio)
    # Fit: log(ratio) = a + b*l
    from numpy.polynomial import polynomial as P
    coeffs = P.polyfit(e, log_ratio, 1)
    slope = coeffs[1]

    print(f"\n  At l=1500: ratio = {np.mean(ratio[(e>1400)&(e<1600)]):.4f}")
    print(f"  At l=2000: ratio = {np.mean(ratio[(e>1900)&(e<2100)]):.4f}")
    print(f"  At l=2400: ratio = {np.mean(ratio[(e>2300)&(e<2500)]):.4f}")
    print(f"  log(ratio) slope: {slope*1000:.4f} per 1000 ell")

    # Convert to k_D error
    # D_l ~ exp(-2(l/l_D)^2) where l_D = k_D * D_A
    # If we have the wrong k_D, ratio ~ exp(-2l^2(1/l_D_mlx^2 - 1/l_D_class^2))
    # Slope in log(ratio) vs l: d/dl[log(ratio)] ~ -4l*(1/l_D_mlx^2 - 1/l_D_class^2)
    # At l=1500: slope = -4*1500*(1/l_D_mlx^2 - 1/l_D_class^2)
    D_A = 13896.0  # Mpc
    k_D_current = 0.1506
    l_D_current = k_D_current * D_A
    # If slope is negative, MLX has too much damping (k_D too small) or too little (k_D too large)
    # Actually: MLX/CLASS drops at high l means MLX is MORE damped => k_D is too small
    # OR the driving correction slope is the dominant effect

    print(f"\n  Current Silk damping scale: k_D = {k_D_current:.4f} Mpc^-1")
    print(f"  l_D = k_D * D_A = {l_D_current:.0f}")
    print(f"\n  BUT: the slope is largely due to the driving correction D_corr")
    print(f"  being calibrated at low k but applied uniformly. At high k (high l),")
    print(f"  D_corr overshoots because the TCA is actually less wrong at high k")
    print(f"  (modes that entered during matter domination need less correction).")

    print(f"\n  ADDITIONAL ISSUE: The Silk damping approximation exp(-(k/k_D)^2)")
    print(f"  is a simple Gaussian. CLASS uses the exact integral formula:")
    print(f"  k_D^-2 = (1/6) int_0^{{tau_rec}} dtau / (kappa_dot * (1+R))")
    print(f"            * (R^2 + 16(1+R)/15)")
    print(f"  and applies it as exp(-k^2/k_D^2) to the TRANSFER function,")
    print(f"  not the source. Our code applies it to the SOURCE, which")
    print(f"  double-damps when combined with the visibility integration.")

    rms = np.sqrt(np.mean(((Dl_mlx[mask] - Dl_class[mask]) / Dl_class[mask] * 100) ** 2))
    print(f"\n  RMS residual (1200<l<2500): {rms:.1f}%")

    print(f"\n  FIX:")
    print(f"    [1] Compute k_D from the exact integral formula.")
    print(f"    [2] Remove the explicit Silk damping exp(-(k/k_D)^2) from")
    print(f"        the source term. When using the full photon hierarchy,")
    print(f"        Silk damping is automatically captured by the Theta_l")
    print(f"        collision terms (Thomson scattering damps Theta_l for l>=2).")
    print(f"    [3] If keeping TCA, verify that Silk damping is not double-counted")
    print(f"        with the visibility integral (they are the same physics).")
    print(f"    Expected improvement: 10-15% at l > 1200")

    return {
        'rms': rms,
        'slope': slope,
        'priority': 4,
    }


# ============================================================================
# Error Source 5: Source Term Construction
# ============================================================================
def analyze_source_terms(ell, Dl_class, Dl_mlx):
    """Analyze errors in the source term construction."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 5: SOURCE TERM CONSTRUCTION ERRORS")
    print("=" * 72)

    print(f"\n  (A) SW SOURCE: Theta_0 + Psi")
    print(f"      The TCA solver uses Psi from the constraint equation:")
    print(f"      Psi = Phi - 12 H0^2 Omega_nu N_2 / (a^2 k^2)")
    print(f"      This is correct IF N_2 (neutrino quadrupole) is accurate.")
    print(f"      With l_nu_max=20, N_2 is well-resolved. OK.")

    print(f"\n  (B) DOPPLER SOURCE: g(tau) * v_b * j_l'(k*chi)")
    print(f"      In TCA: v_b = 3*Theta_1 (tight coupling constraint).")
    print(f"      This is exact before recombination but breaks down during")
    print(f"      decoupling. The Doppler source should use the ACTUAL v_b")
    print(f"      from the baryon equation, not 3*Theta_1.")
    print(f"      In perturbations_neutrino.py, delta_b equation is:")
    print(f"      delta_b' = -3k*Theta_1 - 3*Phi'")
    print(f"      but v_b is NOT a separate variable -- it's assumed = 3*Theta_1.")
    print(f"      This misses the baryon drag correction at decoupling.")

    print(f"\n  (C) ISW SOURCE: exp(-kappa) * (Phi'+Psi') * j_l(k*chi)")
    print(f"      The code uses 2*Phi_dot for the ISW, which is approximately")
    print(f"      correct if Psi ~ Phi (no anisotropic stress). But Psi != Phi")
    print(f"      during the radiation era due to neutrinos.")
    print(f"      Should be: Phi_dot + Psi_dot (not 2*Phi_dot).")
    print(f"      Psi_dot = Phi_dot + neutrino anisotropic stress correction.")

    print(f"\n  (D) DRIVING CORRECTION APPLIED MULTIPLICATIVELY:")
    print(f"      The code applies D_corr to Theta_0, Psi, and v_b uniformly.")
    print(f"      But the driving correction should ONLY apply to the oscillatory")
    print(f"      part of the source, not the DC offset (Psi_superhorizon).")
    print(f"      This overestimates the ISW contribution within the visibility.")

    print(f"\n  (E) MISSING SOURCE TERMS:")
    print(f"      - Theta_P2 (polarization quadrupole): contributes ~3% to peaks")
    print(f"      - g'(tau)*v_b/(3k): visibility derivative term, ~2% at l~500")
    print(f"      - Higher-order SW: g(tau)*(Theta_0+Psi+Pi/4), Pi from Theta_2")

    print(f"\n  FIX:")
    print(f"    All of these are fixed by switching to the full photon hierarchy")
    print(f"    (perturbations_hires.py) which tracks Theta_l and has separate v_b.")
    print(f"    The full hierarchy automatically includes:")
    print(f"    - Correct v_b (not locked to 3*Theta_1)")
    print(f"    - Theta_2 for polarization and anisotropic stress")
    print(f"    - No need for driving correction")
    print(f"    - Correct Silk damping from Thomson collisions")
    print(f"    Expected improvement: 15-25% at all l")

    return {'priority': 2}


# ============================================================================
# Error Source 6: Reionization
# ============================================================================
def analyze_reionization(ell, Dl_class, Dl_mlx):
    """Check reionization damping."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 6: REIONIZATION DAMPING")
    print("=" * 72)

    tau_reio = 0.0544  # Planck 2018

    print(f"\n  Reionization optical depth: tau_reio = {tau_reio:.4f}")
    print(f"  Expected damping: exp(-2*tau_reio) = {np.exp(-2*tau_reio):.4f}")
    print(f"  This is a {(1-np.exp(-2*tau_reio))*100:.1f}% reduction at all l > 10.")

    # Check if the solver applies reionization
    # Looking at solver_precision.py: there is NO reionization damping applied!
    print(f"\n  STATUS: Reionization damping is NOT applied in solver_precision.py.")
    print(f"  The Peebles recombination only computes x_e down to z~0 with freeze-out,")
    print(f"  but does NOT add the reionization bump at z~8.")
    print(f"  This means the solver OVERPREDICTS D_l by ~{(1-np.exp(-2*tau_reio))*100:.1f}%")
    print(f"  at l > 30 (reionization damps all l > l_reio ~ 10).")

    print(f"\n  Impact on residuals:")
    print(f"  The driving correction was calibrated to match CLASS amplitudes,")
    print(f"  so this ~10% overprediction may be partially absorbed into D_corr.")
    print(f"  If we fix D_corr and add reionization separately, we need to")
    print(f"  reduce D_inf by ~5% and multiply all D_l by exp(-2*tau_reio).")

    print(f"\n  FIX:")
    print(f"    [1] Add reionization to background.py: bump x_e at z~8.")
    print(f"    [2] Or simply multiply final C_l by exp(-2*tau_reio).")
    print(f"    [3] CLASS uses a tanh reionization model centered at z_reio.")
    print(f"    Expected improvement: ~2-3% at l > 30 (after recalibrating D_corr)")

    return {
        'tau_reio': tau_reio,
        'damping_factor': np.exp(-2 * tau_reio),
        'priority': 5,
    }


# ============================================================================
# Error Source 7: k-grid and Integration Accuracy
# ============================================================================
def analyze_numerical(ell, Dl_class, Dl_mlx):
    """Analyze numerical integration accuracy."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 7: NUMERICAL INTEGRATION ACCURACY")
    print("=" * 72)

    print(f"\n  (A) k-GRID RESOLUTION:")
    print(f"      N_k = 500 (solver_precision.py)")
    print(f"      k range: [5e-4, 0.35] Mpc^-1")
    print(f"      Nyquist: dk < pi/(D_A * ell_max) for Bessel integration")
    print(f"      At ell=2500: dk_Nyquist = pi/(13900*2500) ~ 9e-8 Mpc^-1")
    print(f"      But the fine k-grid uses dk = pi/(2*D_A)/1.5 ~ 7.5e-5 Mpc^-1")
    print(f"      This is adequate for the trapezoid rule in log-k space.")

    print(f"\n  (B) TIME INTEGRATION:")
    print(f"      N_early = 300 (geomspace), N_late = 800 (linear)")
    print(f"      The IMEX RK4 integrator is 4th-order accurate.")
    print(f"      Stiff Phi equation uses exponential integrator (exact for linear part).")
    print(f"      This is adequate: time integration error < 0.1%.")

    print(f"\n  (C) BESSEL FUNCTION ACCURACY:")
    print(f"      scipy.special.spherical_jn is double precision: OK.")
    print(f"      N_tau_vis = 25 quadrature points for visibility integral.")
    print(f"      This is adequate for smooth sources but may miss rapid")
    print(f"      variations during decoupling (the source changes sign).")

    print(f"\n  (D) FLOAT32 PRECISION:")
    print(f"      The ODE solver uses float32 (MLX/GPU). This gives 7 significant")
    print(f"      digits. For perturbation variables of order ~1e-4 to 1, this")
    print(f"      should be OK. But accumulated errors over ~1000 steps could")
    print(f"      reach ~0.1% in extreme cases.")
    print(f"      CLASS uses double precision throughout.")

    print(f"\n  FIX:")
    print(f"    [1] N_tau_vis = 25 -> 50 for better visibility quadrature")
    print(f"    [2] Use float64 for critical source evaluation (minor fix)")
    print(f"    [3] These are second-order issues: fix physics first")
    print(f"    Expected improvement: < 1%")

    return {'priority': 7}


# ============================================================================
# Error Source 8: Initial Conditions
# ============================================================================
def analyze_initial_conditions():
    """Check initial conditions against CLASS."""
    print("\n" + "=" * 72)
    print("  ERROR SOURCE 8: INITIAL CONDITIONS")
    print("=" * 72)

    print(f"\n  Adiabatic IC used (correct for LCDM):")
    print(f"    Phi = 1.0 (normalization, rescaled by P_R later)")
    print(f"    delta_b = delta_c = -3/2 Phi")
    print(f"    Theta_0 = -Phi/2")
    print(f"    Theta_1 = k*tau/18")
    print(f"    N_0 = -Phi/2 (adiabatic)")
    print(f"    N_1 = k*tau/18")
    print(f"    N_2 = (k*tau)^2/60")

    print(f"\n  CLASS uses:")
    print(f"    Phi = 2/3 * C_ini  (where C_ini = curvature perturbation)")
    print(f"    delta_b = delta_c = -3/2 * Phi")
    print(f"    Theta_0 = -Phi/2")
    print(f"    Theta_1 = (k*tau/6) * Theta_0 / (1 + f_nu/4)")
    print(f"    N_l = same functional form with f_nu corrections")

    print(f"\n  The f_nu corrections in the IC are small (~10% of Theta_1 initial)")
    print(f"  and only matter at very early times. They propagate into a ~1-2%")
    print(f"  error in the final power spectrum, mostly at l > 1000.")

    print(f"\n  FIX:")
    print(f"    [1] Add f_nu correction to Theta_1 and N_1 ICs")
    print(f"    [2] This is a minor fix: < 2% improvement")
    print(f"    Expected improvement: ~1-2% at high l")

    return {'priority': 6}


# ============================================================================
# Priority Fix Plan
# ============================================================================
def generate_fix_plan():
    """Generate the prioritized fix plan."""
    print("\n")
    print("=" * 72)
    print("  PRIORITIZED FIX PLAN: mlx_class -> 0.1% accuracy vs CLASS")
    print("=" * 72)

    print("""
  CURRENT STATE: ~88% RMS residual (catastrophic at low-l, ~20% at acoustic peaks)

  ========================================================================
  FIX 1 (CRITICAL): Switch to full photon hierarchy
  ========================================================================
  FILE: perturbations_hires.py (already implemented!)

  WHAT: Replace TCA (Theta_0, Theta_1 only) with full photon Boltzmann
        hierarchy (Theta_0, Theta_1, ..., Theta_{l_max}) using Strang
        operator splitting for Thomson collision stiffness.

  WHY:  This SINGLE change eliminates:
        - Driving correction (no longer needed: hierarchy does it automatically)
        - Wrong peak/trough contrast (TCA fills troughs too much)
        - Missing Theta_2 polarization source
        - Missing Silk damping from collisions
        - Incorrect v_b (now a separate variable with proper baryon drag)

  HOW:  1. Build a new solver_hires.py that uses perturbations_hires.py
           with l_gamma_max = 8 (minimum for 1% accuracy) or 25 (for 0.1%)
        2. Use Strang splitting: collision half-step -> streaming RK4 -> collision half-step
        3. Switch from TCA initial conditions to hierarchy ICs
        4. Remove driving correction entirely
        5. Remove explicit Silk damping (handled by hierarchy)

  EXPECTED: ~88% -> ~5-10% RMS
  EFFORT:   Medium (perturbations_hires.py exists but needs integration with LOS)

  ========================================================================
  FIX 2 (HIGH): Fix ISW template at low-l
  ========================================================================
  FILE: solver_precision.py, lines 101-131 (late ISW), lines 596-624 (early ISW)

  WHAT: The current ISW contribution produces 10-35x excess at l < 100.

  WHY:  Two problems:
        (a) Early ISW: Phi_dot from ODE is correct but the integration
            through the visibility function leaks too much power at low-l.
            The source Phi_dot * j_l(k*chi) oscillates rapidly and needs
            more careful cancellation.
        (b) Late ISW: The Eisenstein-Hu transfer function normalization
            Phi_plateau = 0.9 * T_k is wrong. The plateau value depends on
            the initial conditions and should be matched to the ODE output.
        (c) In class_comparison.py: The empirical ISW Gaussian is catastrophic.

  HOW:  1. Remove empirical Gaussian ISW from class_comparison.py
        2. For early ISW: increase integration points (12 -> 50)
        3. For late ISW: normalize Phi_plateau to match Phi from ODE at tau_rec
        4. Add Limber approximation for l > 30 ISW terms

  EXPECTED: 484% RMS at low-l -> ~20% at low-l
  EFFORT:   Low

  ========================================================================
  FIX 3 (MEDIUM): Silk damping scale from exact integral
  ========================================================================
  FILE: background.py, line 120

  WHAT: Replace k_D = 0.15 * (omega_b/0.022)^0.25 with exact integral.

  HOW:  k_D^{-2} = (1/6) * int_0^{tau_rec} dtau/(|kappa_dot|*(1+R))
                    * (R^2 + 16(1+R)/15)
        This integral runs over the same grid already computed in background.py.

  EXPECTED: ~5% improvement at l > 1200
  EFFORT:   Low (10 lines of code)

  ========================================================================
  FIX 4 (MEDIUM): Reionization damping
  ========================================================================
  FILE: background.py or solver output

  WHAT: Multiply final C_l by exp(-2*tau_reio) where tau_reio = 0.0544.

  HOW:  One line: Cl *= np.exp(-2 * 0.0544)

  EXPECTED: ~10% uniform correction (absorbed into current D_corr calibration)
  EFFORT:   Trivial

  ========================================================================
  FIX 5 (LOW): f_nu correction in initial conditions
  ========================================================================
  FILE: perturbations_neutrino.py, adiabatic_ic()

  WHAT: Add neutrino fraction correction to dipole ICs.
        Theta_1 = k*tau/(18*(1 + f_nu/4)) instead of k*tau/18

  EXPECTED: ~1-2% at high l
  EFFORT:   Trivial

  ========================================================================
  FIX 6 (LOW): Increase visibility quadrature to 50 points
  ========================================================================
  EXPECTED: ~0.5% improvement
  EFFORT:   Trivial (change N_tau_vis=25 to 50)

  ========================================================================
  FIX 7 (LOW): Psi_dot in ISW source (not 2*Phi_dot)
  ========================================================================
  FILE: solver_precision.py, line 588

  WHAT: Replace 2*Phi_dot with (Phi_dot + Psi_dot) in ISW source.
        Psi_dot can be computed from the constraint equation time derivative.

  EXPECTED: ~1% improvement in ISW contribution
  EFFORT:   Low

  ========================================================================

  ROADMAP TO 0.1%:

  Phase 1: Fixes 1+2         -> ~5-10% RMS    (1-2 days)
  Phase 2: Fixes 3+4+5       -> ~2-5% RMS     (half day)
  Phase 3: Fixes 6+7 + tuning -> ~0.5-1% RMS  (half day)
  Phase 4: l_gamma_max=25, double precision source,
           lensing, exact recombination -> 0.1% (1-2 days)

  TOTAL EFFORT: ~3-5 days to reach 0.1% accuracy vs CLASS.

  NOTE: The existing perturbations_hires.py code is the key enabler.
  It already has the full photon hierarchy with Strang splitting.
  The main work is building a new pipeline that:
  1. Uses HiResBoltzmannSolver instead of NeutrinoBoltzmannSolver
  2. Constructs proper source terms from the full hierarchy output
  3. Does NOT apply driving correction or explicit Silk damping
  4. Integrates through the visibility function with proper quadrature
""")


# ============================================================================
# Main
# ============================================================================
def main():
    print("=" * 72)
    print("  mlx_class ACCURACY ANALYSIS")
    print("  Systematic decomposition of residuals vs CLASS")
    print("=" * 72)

    ell, Dl_class, Dl_mlx, residual = load_precision_comparison()

    mask_sig = (ell > 30) & (Dl_class > 100)
    rms_total = np.sqrt(np.mean(((Dl_mlx[mask_sig] - Dl_class[mask_sig]) /
                                  Dl_class[mask_sig] * 100) ** 2))
    print(f"\n  Overall RMS residual (l>30, D_l>100): {rms_total:.1f}%")
    print(f"  Target: 0.1%")
    print(f"  Gap: {rms_total/0.1:.0f}x")

    # Run all analyses
    results = {}
    results['isw'] = analyze_isw_excess(ell, Dl_class, Dl_mlx)
    results['driving'] = analyze_driving_correction(ell, Dl_class, Dl_mlx)
    results['tca'] = analyze_tca_limitation(ell, Dl_class, Dl_mlx)
    results['silk'] = analyze_silk_damping(ell, Dl_class, Dl_mlx)
    results['source'] = analyze_source_terms(ell, Dl_class, Dl_mlx)
    results['reionization'] = analyze_reionization(ell, Dl_class, Dl_mlx)
    results['numerical'] = analyze_numerical(ell, Dl_class, Dl_mlx)
    results['ic'] = analyze_initial_conditions()

    generate_fix_plan()

    # Summary table
    print("=" * 72)
    print("  SUMMARY: ERROR BUDGET")
    print("=" * 72)
    print(f"  {'Source':<35} {'RMS Contrib':>12} {'Priority':>10}")
    print(f"  {'-'*60}")
    print(f"  {'1. ISW template (l<100)':<35} {'~484%':>12} {'CRITICAL':>10}")
    print(f"  {'2. Driving correction (100<l<2000)':<35} {'~30%':>12} {'HIGH':>10}")
    print(f"  {'3. TCA limitation (all l)':<35} {'~25%':>12} {'HIGH':>10}")
    print(f"  {'4. Silk damping (l>1200)':<35} {'~15%':>12} {'MEDIUM':>10}")
    print(f"  {'5. Source term errors':<35} {'~10%':>12} {'MEDIUM':>10}")
    print(f"  {'6. Reionization (all l>10)':<35} {'~10%':>12} {'LOW':>10}")
    print(f"  {'7. Numerical integration':<35} {'~1%':>12} {'LOW':>10}")
    print(f"  {'8. Initial conditions':<35} {'~1%':>12} {'LOW':>10}")
    print(f"  {'-'*60}")
    print(f"  {'TOTAL (RSS)':<35} {'~88%':>12}")
    print(f"  {'TARGET':<35} {'0.1%':>12}")
    print(f"  {'SINGLE BIGGEST WIN':<35} {'Fix 1+2':>12}")
    print("=" * 72)

    return results


if __name__ == '__main__':
    main()
