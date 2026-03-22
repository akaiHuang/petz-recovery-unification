# CLASS Line-of-Sight Integration & Thermodynamics: Technical Reference

Extracted from `gdm_class_public/source/transfer.c` (5127 lines) and
`gdm_class_public/source/thermodynamics.c` (4344 lines).

---

## 1. Line-of-Sight (LOS) Integration

### 1.1 The Master Formula

CLASS computes the transfer function Delta_l^X(q) via:

```
Delta_l(k) = integral_0^{tau_0} d(tau) S(k, tau) * Phi_l(k * chi)
```

where:
- `S(k, tau)` = source function (from perturbation module)
- `Phi_l(x)` = generalized radial function (spherical Bessel j_l(x) in flat case)
- `chi = tau_0 - tau` = comoving distance (conformal lookback time)
- The integration variable is `tau0_minus_tau`

The code uses the variable `q` throughout the transfer module (related to `k` by
`q^2 = k^2 + K(1+m)` where K = spatial curvature, m=0,1,2 for scalar/vector/tensor;
in flat space q=k).

### 1.2 Quadrature Method: Trapezoidal Rule

The integral is computed by **trapezoidal convolution** in `transfer_integrate()`
(line 3183), which calls `array_trapezoidal_convolution()`:

```c
// tools/arrays.c line 3194
double res = 0.0;
for (i = 0; i < n; i++) {
    res += sources[i] * radial_function[i] * w_trapz[i];
}
```

where `w_trapz[i]` are precomputed trapezoidal weights (half-intervals at endpoints).
This is a simple weighted dot product of source * Bessel * weights.

**Bessel truncation correction**: When the integral is truncated because the Bessel
function is only sampled above some x_min, CLASS corrects the trapezoidal error
at the truncation boundary (line 3334-3337):

```c
if ((index_tau_max != tau_size-1) && (index_tau_max == index_tau_max_Bessel)) {
    *trsf -= 0.5 * (tau0_minus_tau[index_tau_max+1] - tau0_minus_tau_min_bessel)
           * radial_function[index_tau_max] * sources[index_tau_max];
}
```

### 1.3 Source Function Interpolation

Sources are computed on the perturbation module's k-grid and interpolated to the
transfer module's q-grid using **cubic spline interpolation** in k-space
(`transfer_interpolate_sources()`, line 1999):

```c
interpolated_sources[index_tau] =
    a * S[tau, k_i] + b * S[tau, k_{i+1}]
    + ((a^3 - a) * S''[tau, k_i] + (b^3 - b) * S''[tau, k_{i+1}]) * h^2/6
```

where `a = 1-b`, `b = (q - k_i)/(k_{i+1} - k_i)`, standard cubic spline formula.

### 1.4 Limber Approximation

For high multipoles, CLASS switches to the Limber approximation instead of full
integration. The switch is controlled by `transfer_use_limber()` (line 3094):

| Type | Switch condition |
|------|-----------------|
| CMB lensing potential | `l > l_switch_limber` (default: 10) |
| Number counts (local) | `l >= l_switch_limber_for_nc_local_over_z * z_mean` (default factor: 100) |
| Number counts (LOS) | `l >= l_switch_limber_for_nc_los_over_z * z_mean` (default factor: 30) |
| Any type | `q > q_max_bessel` (forced Limber when Bessels not computed) |

The Limber formula for SCALAR_TEMPERATURE_0 (line 3386-3431):

```
tau0_minus_tau_limber = (l + 0.5) / q       [flat case]

Delta_l(q) = sqrt(pi/(2l)) * (1 - 1/(4l) + 1/(32l^2)) * S(tau_limber) / (l + 0.5)
```

The source at `tau_limber` is obtained by **parabolic interpolation** (3-point
polynomial fit) through `transfer_limber_interpolate()` (line 3521). It interpolates
`S * (tau0 - tau)` rather than `S` alone to handle the lensing divergence at tau=tau0.

### 1.5 Number of tau Points

The tau sampling is determined dynamically in `transfer_sources()` and depends on:
- The source type (CMB, lensing, number counts)
- The Bessel oscillation frequency: `tau_size >= (tau_max - tau_min) / (wavelength/l) * selection_sampling_bessel`
- For number count bins: the width and shape of the selection function
- The late-source cut: `ptw->tau0_minus_tau_cut = tau0 - pth->tau_cut` allows ignoring sources at very recent times

### 1.6 OpenMP Parallelization

The main loop over wavenumber index_q is parallelized with OpenMP dynamic scheduling
(line 314, 342):

```c
#pragma omp parallel shared(...) private(ptw, index_q, ...)
{
    // allocate per-thread workspace
    #pragma omp for schedule(dynamic)
    for (index_q = 0; index_q < ptr->q_size; index_q++) { ... }
}
```

Each thread has its own workspace containing sources, tau0_minus_tau, w_trapz, etc.

---

## 2. Bessel Function Computation

### 2.1 Flat Case: Precomputed Hyperspherical Table

In the flat case (K=0), CLASS precomputes a global table of spherical Bessel functions
j_l(x) using the **hyperspherical function infrastructure** (HIS). The table is
created once during initialization (line 267):

```c
hyperspherical_HIS_create(
    sgnK=0,      // flat space
    K=1.,
    l_size_max,  // number of multipoles
    l[],         // multipole values
    hyper_x_min, // default 1e-5
    x_max,       // q_max * tau0
    hyper_sampling_flat,  // default 8.0 points per wavelength 2pi
    l_max + 1,
    hyper_phi_min_abs,    // default 1e-10
    &BIS,        // output: Bessel Interpolation Structure
    error_message
);
```

Key parameters:
- **x_min**: 1e-5 (below this, j_l(x) ~ 0)
- **Sampling**: ~8 points per oscillation wavelength 2pi (controlled by `hyper_sampling_flat = 8.0`)
- **phi_min_abs**: 1e-10 (threshold for first non-negligible point)

### 2.2 Interpolation Method: Hermite Polynomials

The precomputed Bessel values are interpolated using **Hermite polynomial interpolation**
of configurable order (line 3857-3878):

| Case | Hermite Order | When Used |
|------|---------------|-----------|
| Flat (K=0) | HERMITE4 (4th order) | All flat-space q modes |
| Curved, exact | HERMITE6 (6th order) | Low-q modes needing exact hyperspherical functions |
| Curved, flat approx | HERMITE4 | High-q modes where flat approximation suffices |

The interpolation functions are:
- `hyperspherical_Hermite4_interpolation_vector_Phi` for Phi_l(x)
- `hyperspherical_Hermite4_interpolation_vector_dPhi` for Phi_l'(x)
- etc. for combinations (Phi, dPhi, d2Phi)

### 2.3 Curved-Space Handling

For K != 0, the code uses a two-regime strategy:
1. **Low q (nu < hyper_flat_approximation_nu)**: Compute exact hyperspherical Bessel
   functions with HERMITE6 interpolation
2. **High q**: Use flat Bessel functions with a **rescaling trick** (line 3842-3908):
   - Rescale the argument: `chi_rescaled = chi * sqrt(l(l+1)) / chi_turning_point`
   - Rescale the amplitude with a fitted correction factor involving
     `(1 - K*l(l+1)/q^2)^{-1/12}` and empirical atan-based corrections

### 2.4 Radial Function Types

Different source types require different combinations of Bessel functions (line 3926-4041):

| Type | Radial function |
|------|----------------|
| SCALAR_TEMPERATURE_0 | Phi_l(chi) = j_l(k*chi) |
| SCALAR_TEMPERATURE_1 | sqrt(K/k) * dPhi_l/dchi |
| SCALAR_TEMPERATURE_2 | (3K/k^2 * d^2Phi/dchi^2 + Phi) / (2*sqrt(1-3K/k^2)) |
| SCALAR_POLARISATION_E | sqrt(3/8 * (l+2)(l+1)l(l-1)) / s2 * csc^2(chi) * Phi |
| TENSOR_TEMPERATURE_2 | sqrt(3/8 * (l+2)(l+1)l(l-1)) / (si * ssqrt2) * csc^2(chi) * Phi |
| NC_RSD | (K/k^2) * d^2Phi/dchi^2 |

---

## 3. k-Mode (q-Mode) Sampling

### 3.1 Strategy: Log-to-Linear Transition

The q-sampling is defined in `transfer_get_q_list()` (line 1016). The step size
transitions smoothly from logarithmic at small q to linear at large q:

```c
q_next = q_prev + q_period * q_linstep * q_prev / (q_prev + q_linstep / q_logstep_spline);
```

This formula gives:
- **Small q limit** (`q << q_linstep/q_logstep_spline`): delta_q/q = q_period * q_logstep_spline (logarithmic)
- **Large q limit** (`q >> q_linstep/q_logstep_spline`): delta_q = q_period * q_linstep (linear)

### 3.2 Default Precision Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `q_linstep` | 0.45 | Asymptotic linear step in q space |
| `q_logstep_spline` | 170.0 | Initial log step: delta_q/q ~ q_period/170 |
| `q_logstep_trapzd` | 20.0 | Log step for closed-space integer nu |
| `q_logstep_open` | 6.0 | Power-law modification in open models |

### 3.3 Range

- **q_min**: `ppt->k_min` (from perturbation module, typically ~1e-5/Mpc)
- **q_max**: `max over modes of ppt->k[mode][k_size_cl[mode]-1]` (typically ~0.3-1/Mpc for TT, higher for lensing)

### 3.4 Closed Universe Special Treatment

For K > 0 (closed), q = nu * sqrt(K) where nu is an integer at low q. The code:
1. Steps through integer nu values with `q_logstep_trapzd` spacing
2. Transitions to continuous (non-integer) sampling above `hyper_flat_approximation_nu`
3. Uses a smooth `q_numstep_transition` zone to merge the two regimes

### 3.5 k from q

After the q-list is built, k values are derived (line 1273):
```c
k[index_q] = sqrt(q[index_q]^2 - K*(m+1))
```
where m=0 (scalars), 1 (vectors), 2 (tensors).

---

## 4. Recombination

### 4.1 Two Recombination Engines

CLASS supports two recombination codes:

1. **RECFAST v1.5** (default): Based on Seager, Sasselov & Scott (ApJ 523 L1, 1999)
   with fudge updates from Wong, Moss & Scott (2008).
2. **HyRec**: Ali-Haimoud & Hirata, more accurate multi-level atom calculation.
   Enabled with `recombination = HyRec` in input, requires compilation with `-DHYREC`.

### 4.2 RECFAST Implementation Details

Located in `thermodynamics_recombination_with_recfast()` (line 3457).

**ODE system**: Three coupled variables y[0], y[1], y[2]:
- `y[0] = x_H` (hydrogen ionization fraction)
- `y[1] = x_He` (helium ionization fraction)
- `y[2] = T_mat` (matter temperature)

**Initial conditions** (at z = z_initial = 10^4):
```c
y[0] = 1.0;    // fully ionized hydrogen
y[1] = 1.0;    // fully ionized helium
y[2] = T_cmb * (1 + z);  // matter temperature = radiation temperature
```

**Integration**: Nz = `recfast_Nz0` = 20000 uniform steps from z = 10^4 to z = 0.
Step size: delta_z = z_initial / Nz = 0.5. Uses a generic integrator (ndf15-like
stiff ODE solver from `generic_integrator_workspace`).

### 4.3 Hydrogen Recombination (Peebles Equation with Fudge)

The hydrogen evolution equation (line 4045):

```
dx_H/dz = (x * x_H * n * R_down - R_up * (1-x_H) * exp(-CL/T_mat)) * C / (H*(1+z))
```

where the **Peebles C-factor** (line 4020):
```
C = (1 + K * Lambda * n * (1-x_H)) / (1/fu + K * Lambda * n * (1-x_H)/fu + K * R_up * n * (1-x_H))
```

### 4.4 Fudge Factors

| Parameter | Default | Source |
|-----------|---------|--------|
| `recfast_fudge_H` | 1.14 | RECFAST v1.4, hydrogen |
| `recfast_delta_fudge_H` | -0.015 | RECFAST v1.5.2, correction when Hswitch=TRUE |
| **Effective H fudge** | **1.125** | = 1.14 - 0.015 (default with Hswitch=TRUE) |
| `recfast_fudge_He` | 0.86 | RECFAST v1.4, helium |

Additional Gaussian corrections to the Peebles fudge factor (from v1.5):
- AGauss1 = -0.14, center z=7.28, width 0.18
- AGauss2 = 0.079, center z=6.73, width 0.33

### 4.5 Helium Recombination

RECFAST handles helium with multiple levels of sophistication controlled by
`recfast_Heswitch` (default: 6, most accurate):

- Heswitch >= 3: includes triplet states with solenoid feedback
- Heswitch >= 5: includes H absorption of He photons (feedback from H on He lines)
- Separate treatment for He I recombination (first recombination, z ~ 6000-8000)
  and He II recombination (second, very early)

### 4.6 Matter Temperature Evolution

Two regimes for T_mat (line 4083-4103):

**Tight-coupling regime** (`t_Thomson < H_frac * t_Hubble`):
```
dT/dz = T_cmb + epsilon * correction_terms
```
(smoothed transition as in CAMB, Adam Moss suggestion, RECFAST v1.5)

**Decoupled regime**:
```
dT/dz = CT * T_rad^4 * x / (1+x+fHe) * (T_mat - T_rad) / (H*(1+z)) + 2*T_mat/(1+z)
```
First term: Compton cooling/heating. Second term: adiabatic cooling.

### 4.7 HyRec Parameters

When using HyRec, the integration uses `dlna = 8.49e-5` and
`nz = floor(2 + log((1+z_start)/(1+z_end))/dlna)` steps.
HyRec reads tabulated effective recombination rates from external files:
- `hyrec_Alpha_inf_file`
- `hyrec_R_inf_file`
- `hyrec_two_photon_tables_file`

---

## 5. Reionization

### 5.1 Default Model: tanh (CAMB-like)

The default reionization model (`reio_camb`) uses a tanh-shaped transition
(line 1976-1990):

```
argument = ((1+z_reio)^p - (1+z)^p) / (p * (1+z_reio)^(p-1)) / Delta_z

x_e(z) = (x_after - x_before) * (tanh(argument) + 1) / 2 + x_before
```

where:
- `z_reio` = reionization redshift (input or derived from tau_reio)
- `p` = reionization_exponent (default: 1.5)
- `Delta_z` = reionization_width (default: 0.5)
- `x_after` = 1 + Y_He / (4*(1-Y_He)) (H fully ionized + He singly ionized)
- `x_before` = residual x_e from recombination table

### 5.2 Helium Full Reionization

A second tanh for helium double-ionization (line 2001-2007):
```
argument = (z_He_fullreio - z) / width_He_fullreio
x_e += f_He * (tanh(argument) + 1) / 2
```
with `f_He = Y_He / (4*(1-Y_He))` and default `z_He_fullreio = 3.5`, `width = 0.5`.

### 5.3 Determining z_reio from tau_reio

When the input is `tau_reio` (optical depth) rather than `z_reio`, CLASS finds
`z_reio` by **bisection** (dichotomy, line 2331-2415):

```
while (tau_sup - tau_inf > tau_reio * reionization_optical_depth_tol):
    z_mid = (z_sup + z_inf) / 2
    compute tau for z_mid
    update bounds
```

### 5.4 Alternative Reionization Models

CLASS supports several parametrizations:
- `reio_camb`: standard tanh (described above)
- `reio_half_tanh`: half-tanh (no helium double ionization)
- `reio_bins_tanh`: binned history with tanh jumps between bins
- `reio_many_tanh`: multiple tanh jumps
- `reio_inter`: fully interpolated xe(z) from input table

### 5.5 Adaptive Sampling

The reionization table uses **adaptive step size** in z (line 2878-2957):

```c
relative_variation = |delta(dkappa/dz)/kappa'_z| + |delta(dkappa/dtau)/kappa'_tau|

if (relative_variation < reionization_sampling):
    accept step; dz = min(0.9 * (target/variation), 5) * dz
else:
    reject step; dz = 0.9 * (target/variation) * dz
```

with `dz <= dz_max` (inherited from recombination table spacing).

---

## 6. Opacity and Visibility Function

### 6.1 Thomson Scattering Rate kappa'(tau)

The differential optical depth is stored as `dkappa/dtau` (line 125-126):

```
dkappa/dtau = a^{-2} * n_e(today) * x_e * sigma_T * (Mpc/m)
            = (1+z)^2 * n_e * x_e * sigma_T * (Mpc/m)
```

This is computed at each z in both the recombination and reionization tables and
stored in `thermodynamics_table[..., index_th_dkappa]`.

### 6.2 Derivatives of kappa

After the recombination/reionization tables are merged:

1. **d^2 kappa/dtau^2** (ddkappa): Computed by cubic spline differentiation of dkappa
   (line 680-700):
   - First, spline the dkappa column to get second derivatives (d^3 kappa/dtau^3)
   - Then derive the spline to get ddkappa

2. **d^3 kappa/dtau^3** (dddkappa): Used as spline coefficients for dkappa

### 6.3 Optical Depth kappa(tau)

The optical depth is computed by **spline integration** (line 704-714):

```
-kappa(tau) = integral_{tau_today}^{tau} d(tau') * (dkappa/dtau')
```

stored temporarily in the `g` column (later overwritten with the visibility function).
Uses `array_integrate_spline_table_line_to_line()` which performs cubic spline quadrature
(much more accurate than raw trapezoidal).

### 6.4 Visibility Function g(tau)

Computed exactly (line 745-776):

```
g(tau) = (dkappa/dtau) * exp(-kappa(tau))
```

with derivatives:

```
g'  = (kappa'' + (kappa')^2) * exp(-kappa)

g'' = (kappa''' + 3*kappa'*kappa'' + (kappa')^3) * exp(-kappa)
```

Note: these are the exact analytical derivatives, not numerical differentiation.

### 6.5 exp(-kappa) Separately

```
exp(-kappa) is stored separately in index_th_exp_m_kappa
```

### 6.6 Recombination Redshift z_rec

Defined as the **maximum of g(tau)**, found by scanning the table (line 999-1020).
A **cubic interpolation** refines the peak location:

```
z_rec = z[i+1] + 0.5 * dz * (g[i] - g[i+2]) / (g[i] - 2*g[i+1] + g[i+2])
```

### 6.7 z_star (tau = 1 Surface)

Found where exp(-kappa) crosses 1/e (line 1142-1147):
```
scan until exp(-kappa) > 1/e, then linear interpolation
```

### 6.8 Variation Rate (for Perturbation Sampling)

CLASS computes a "variation rate" used to set the time step in the perturbation module
(line 783-788):

```
rate = sqrt( (kappa')^2 + (kappa''/kappa')^2 + |kappa'''/kappa'| )
```

This is smoothed with `array_smooth()` (radius = `thermo_rate_smoothing_radius`).

---

## 7. Damping Scale

The photon diffusion damping scale is computed in thermodynamics_init() (line 599-676):

```
r_d = 2pi * sqrt( integral_{tau_ini}^{tau} dtau * (1/kappa') * (1/6) * (R^2/(1+R) + 16/15) / (1+R) )
```

where R = (3/4) * rho_b / rho_g.

**Boundary correction**: At tau_ini, the analytic contribution is added:
```
r_d = 2pi * sqrt( 16/(15*6*3) * tau_ini/kappa'_ini + numerical_integral )
```

This uses the fact that near tau_ini (deep in radiation domination), R~0 and
kappa' ~ 1/a^2, tau ~ a, giving int_0^{tau_ini} dtau/kappa' = tau_ini/(3*kappa'_ini).

### 7.1 Baryon Drag Epoch

The baryon drag optical depth is computed as:
```
tau_d = integral_{tau_today}^{tau} dtau * (1/R) * kappa'
```

where R = (3/4)*rho_b/rho_g. The drag epoch z_d is found where tau_d = 1.

---

## 8. Numerical Tricks Summary

### 8.1 Spline Everywhere
CLASS uses cubic spline interpolation pervasively:
- Source functions interpolated in k via cubic spline
- Thermodynamic quantities interpolated in z via cubic spline
- Optical depth computed by spline integration (more accurate than trapezoidal)
- Transfer functions interpolated in q via `array_interpolate_two()`

### 8.2 Handling Oscillatory Bessel Integrals
- **Precomputed Bessel table**: ~8 points per 2pi wavelength (hyper_sampling_flat = 8.0)
- **Hermite interpolation**: 4th or 6th order for sub-grid accuracy
- **Truncation at x_min**: j_l(x) ~ 0 for x < x_min(l); skip those tau values
- **Limber approximation**: At high l, replace oscillatory integral with single-point evaluation
- **Bessel truncation correction**: Exact correction for the trapezoidal error at cutoff

### 8.3 Visibility Function Computation
- Compute -kappa by spline integration (not trapezoidal)
- Then g = kappa' * exp(-kappa) with analytical derivatives
- z_rec found by cubic peak interpolation
- Smoothing of variation rate for perturbation time-stepping

### 8.4 Source * (tau0 - tau) Trick
In Limber interpolation, the code interpolates S*(tau0-tau) instead of S alone.
This product is regular at tau = tau0, while S alone can diverge for lensing sources.

### 8.5 Non-Flat Universe Bessel Handling
- At low q: exact hyperspherical functions with HERMITE6 (more points, higher accuracy)
- At high q: flat Bessel with amplitude/argument rescaling (avoids computing expensive
  exact hyperspherical functions for modes where the flat approximation is adequate)
- The transition index `index_q_flat_approximation` is computed once

### 8.6 Parallel Computation
- Main loop over q is OpenMP-parallelized with dynamic scheduling
- Each thread owns its own workspace (sources, tau arrays, Bessel evaluations)
- No inter-thread communication during the q loop

---

## 9. Key Data Flow

```
perturbation module: S(k, tau) on k-grid, tau-grid
    |
    | cubic spline interpolation in k
    v
transfer module: S(q, tau) at each q in q-list
    |
    | construct transfer sources (selection functions, etc.)
    v
sources(tau) for each (q, type)
    |
    |--- For each l: convolve with Phi_l(k*chi)
    |    using trapezoidal quadrature (precomputed weights)
    |    OR Limber approximation at high l
    v
Delta_l(q) stored in transfer table
    |
    | (used by spectra module to compute C_l)
    v
C_l = (2/pi) * integral dq * q^2 * Delta_l(q) * Delta_l(q) * P(k)
```

---

## 10. Key Precision Parameters (Defaults)

| Parameter | Default | Controls |
|-----------|---------|----------|
| `recfast_z_initial` | 10000 | Start of recombination integration |
| `recfast_Nz0` | 20000 | Number of z steps (dz = 0.5) |
| `recfast_fudge_H` | 1.14 | Hydrogen Peebles fudge factor |
| `recfast_fudge_He` | 0.86 | Helium fudge factor |
| `hyper_sampling_flat` | 8.0 | Bessel sampling points per 2pi |
| `hyper_x_min` | 1e-5 | Minimum x for Bessel table |
| `hyper_phi_min_abs` | 1e-10 | Bessel threshold for first point |
| `q_linstep` | 0.45 | Asymptotic linear q step |
| `q_logstep_spline` | 170.0 | Initial log q step parameter |
| `l_switch_limber` | 10 | Limber switch for CMB lensing |
| `l_switch_limber_for_nc_local_over_z` | 100 | Limber switch for number counts |
| `reionization_sampling` | (adaptive) | Relative variation threshold |
| `tol_thermo_integration` | 1e-2 | ODE integration tolerance |
