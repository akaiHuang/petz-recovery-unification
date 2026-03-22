# CLASS ODE Integration Strategy

Extracted from `gdm_class_public/source/perturbations.c` (10,506 lines)
and `gdm_class_public/tools/evolver_ndf15.c`.

---

## 1. High-Level Architecture

The entry point for each (k, IC) pair is `perturb_solve()` (line 2723).
The flow:

```
perturb_solve(k)
  |
  +--> bisection to find tau_ini  (lines 2817-2930)
  |     ensures: tau_c/tau_h < 0.0015 AND k/(aH) < 0.07
  |
  +--> perturb_find_approximation_number()  (line 2938)
  |     counts how many intervals with uniform approximation scheme
  |
  +--> perturb_find_approximation_switches()  (line 2959)
  |     uses bisection to locate each switch time to tolerance tol_tau_approx
  |     result: interval_limit[] and interval_approx[][]
  |
  +--> FOR each interval: (line 3007)
  |     |
  |     +--> set ppw->approx[] to interval_approx[index_interval]
  |     |
  |     +--> perturb_vector_init()
  |     |     - if first interval: set initial conditions
  |     |     - if switch: redistribute variables, seed higher multipoles
  |     |
  |     +--> generic_evolver()  (= evolver_ndf15 or evolver_rk)
  |           integrates from interval_limit[i] to interval_limit[i+1]
  |           calls perturb_derivs() for dy/dtau
  |           calls perturb_timescale() for adaptive stepping
  |           calls perturb_sources() at each tau_sampling point
  |
  +--> fill remaining source array with zeros
```

## 2. Time Intervals and Approximation Scheme

For a typical scalar mode, the time intervals are:

```
tau_ini ---- tau_tca_switch ---- tau_rsa_switch ---- tau_end (= today)
         TCA on, RSA off      TCA off, RSA off      TCA off, RSA on
         (reduced system)      (full hierarchy)      (free-streaming)
```

Additional switch points may exist for:
- **UFA** (ultra-relativistic fluid approximation for massless neutrinos)
- **NCDMFA** (non-cold dark matter fluid approximation)
- **TCA_IDM_DR** (dark matter--dark radiation tight coupling)
- **RSA_IDR** (interacting dark radiation streaming)

The approximation flags are enumerated in `perturbations.h` (lines 29-34):

```c
enum tca_flags   {tca_on, tca_off};         // chronological order
enum rsa_flags   {rsa_off, rsa_on};          // chronological order
enum ufa_flags   {ufa_off, ufa_on};
enum ncdmfa_flags {ncdmfa_off, ncdmfa_on};
```

**Key constraint**: approximations can only switch in one direction (left to right in the enum). The code enforces this and aborts if an approximation tries to go backward.

**Key constraint**: only ONE approximation can switch at any given time. The code checks this (line 3536).

### Workspace indices for approximations (line 631-637 of perturbations.h)

```c
int index_ap_tca;      // tight-coupling approximation
int index_ap_rsa;      // radiation streaming approximation
int index_ap_tca_idm_dr; // dark tight-coupling
int index_ap_rsa_idr;  // dark radiation streaming
int index_ap_ufa;      // ur fluid approximation
int index_ap_ncdmfa;   // ncdm fluid approximation
int ap_size;           // total number of approximation flags
```

## 3. Switch Criteria

### 3.1 Initial time (tau_ini)

Found by bisection (lines 2817-2930). The earliest time satisfying BOTH:

```
tau_c / tau_h < start_small_k_at_tau_c_over_tau_h  (default 0.0015)
k / (aH)     < start_large_k_at_tau_h_over_tau_k   (default 0.07)
```

where:
- `tau_c = 1/kappa'` (photon mean free path in conformal time)
- `tau_h = 1/(aH)` (Hubble time)
- `tau_k = 1/k` (mode time scale)

For ncdm species, also requires `|w - 1/3| < tol_ncdm_initial_w`.

### 3.2 TCA off (tight coupling -> full hierarchy)

Evaluated in `perturb_approximations()` (lines 6064-6072). TCA stays on when BOTH:

```
tau_c / tau_h < tight_coupling_trigger_tau_c_over_tau_h  (default 0.015)
tau_c / tau_k < tight_coupling_trigger_tau_c_over_tau_k  (default 0.01)
```

Equivalently, TCA switches OFF when EITHER:
- `tau_c / tau_h >= 0.015` (photon mean free path becomes comparable to Hubble time)
- `tau_c / tau_k >= 0.01`  (photon mean free path becomes comparable to mode wavelength)

This happens around recombination. For small k (large scales), the tau_h condition triggers first; for large k (small scales), the tau_k condition triggers first.

### 3.3 RSA on (full hierarchy -> radiation streaming)

Evaluated in lines 6104-6112. RSA switches ON when ALL:

```
tau / tau_k > radiation_streaming_trigger_tau_over_tau_k  (default 45.0)
  i.e., k*tau > 45  (mode is deep inside the horizon)

tau > tau_free_streaming  (from thermodynamics module)

radiation_streaming_approximation != rsa_none
```

### 3.4 UFA on (exact ur hierarchy -> fluid approximation)

Line 6147-6154:
```
tau / tau_k > ur_fluid_trigger_tau_over_tau_k  (default 30.0)
  i.e., k*tau > 30
```

### 3.5 NCDMFA on (exact ncdm hierarchy -> fluid approximation)

Line 6159-6166:
```
tau / tau_k > ncdm_fluid_trigger_tau_over_tau_k  (default 31.0)
  i.e., k*tau > 31
```

### 3.6 Summary timeline (typical scalar mode)

```
tau_ini             TCA on, RSA off, UFA off
  |                   (2+2 photon vars: delta_g, theta_g only)
  |                   (2+2 baryon vars: delta_b, theta_b)
  |                   (metric: eta or phi)
  |                   (~10-15 total equations)
  v
tau_tca_switch      TCA OFF -> full photon hierarchy activated
  |                   (l_max_g=12 temperature + l_max_pol_g=10 polarization)
  |                   (~50-70 equations depending on species)
  v
tau_ufa_switch      UFA ON (if ur species present)
  |                   (ur hierarchy truncated to delta, theta, shear only)
  v
tau_rsa_switch      RSA ON -> photon+ur hierarchy removed from ODE
  |                   (delta_g, theta_g computed algebraically)
  |                   (~5-15 equations: metric + matter only)
  v
tau_end = tau_today
```

## 4. TCA Equations (Reduced System)

During TCA, the photon hierarchy is truncated to just `{delta_g, theta_g}`.
Shear and higher multipoles are computed analytically.

### 4.1 Evolved variables during TCA

```
delta_g, theta_g   (photon density and velocity -- integrated)
delta_b, theta_b   (baryon density and velocity -- integrated)
delta_cdm          (CDM density, synchronous gauge)
eta (or phi)       (metric perturbation)
delta_ur, theta_ur, shear_ur, l3_ur, ...  (if ur present, not affected by TCA)
```

### 4.2 Baryon velocity equation during TCA (line 8956-8960)

```
dy[theta_b] = [ -a'/a * theta_b
                + k^2 * (delta_p_b/rho_b + R*(delta_g/4 - s2^2 * sigma_g^TCA))
                + R * slip ] / (1+R)
              + metric_euler
```

where `R = (4/3) * rho_g / rho_b`.

### 4.3 Photon velocity equation during TCA (line 9050-9052)

```
dy[theta_g] = -(dy[theta_b] + a'/a * theta_b - k^2 * delta_p_b/rho_b) / R
              + k^2 * (delta_g/4 - s2^2 * sigma_g^TCA)
              + (1+R)/R * metric_euler
```

### 4.4 TCA slip and shear (`perturb_tca_slip_and_shear`, line 9873)

The function computes two quantities: `slip` and `shear_g`, stored in `ppw->tca_slip` and `ppw->tca_shear_g`.

**Shear at first order in tight coupling** (line 10042):
```
sigma_g = (16/45) * tau_c * (theta_g + metric_shear)
```

**Slip formulas** depend on `ppr->tight_coupling_approximation`:

| Method | Key | Description |
|--------|-----|-------------|
| `first_order_MB` | 0 | Ma & Bertschinger: assumes kappa' ~ a^{-2} |
| `first_order_CAMB` | 1 | Relaxes kappa' ~ a^{-2} assumption |
| `first_order_CLASS` | 2 | Also relaxes c_b^2 ~ a^{-1} |
| `second_order_CRS` | 3 | Full second-order (Cyr-Racine & Sigurdson) |
| `second_order_CLASS` | 4 | Second-order CLASS formulation |
| `compromise_CLASS` | 5 | Leading second-order terms only (default) |

**First-order CAMB/compromise_CLASS slip** (line 10019-10026):
```
slip = (dtau_c/tau_c - 2*a'/a/(1+R)) * (theta_b - theta_g)
     + F * [ -a''/a * theta_b
             + k^2 * (-a'/a * delta_g/2
                      + c_b^2 * (-theta_b - metric_continuity)
                      - (4/3)*(-theta_g - metric_continuity)/4 )
             - a'/a * metric_euler ]
```
where `F = tau_c/(1+R)`.

**Compromise CLASS second-order corrections** (line 10144-10150):
```
slip = (1 - 2*a'/a*F) * slip
     + F*k^2 * (2*a'/a*s2^2*sigma_g + s2^2*sigma_g'
                - (1/3 - c_b^2)*(F*theta' + 2*F'*theta_b))

sigma_g = (1 - 11/6*dtau_c) * sigma_g
        - (11/6)*tau_c*(16/45)*tau_c*(theta' + metric_shear')
```

### 4.5 Seeding at TCA -> Full switch (line 4302-4398)

When TCA switches off, the new variables are seeded:

```c
// Photon shear: from the TCA analytical value
y[shear_g] = ppw->tca_shear_g;

// l=3: second-order TCA
y[l3_g] = (6/7) * k/kappa' * s_l[3] * shear_g;

// Polarization l=0: first-order TCA
y[pol0_g] = 2.5 * shear_g;

// Polarization l=1: second-order TCA
y[pol1_g] = k/kappa' * (5 - 2*s_l[2])/6 * shear_g;

// Polarization l=2: first-order TCA
y[pol2_g] = 0.5 * shear_g;

// Polarization l=3: second-order TCA
y[pol3_g] = k/kappa' * 3*s_l[3]/14 * shear_g;

// Perturbed recombination (if enabled)
y[delta_temp] = (1/3) * delta_b;
y[delta_chi]  = 0;
```

All higher multipoles (l >= 4 for temperature, l >= 4 for polarization) are set to zero.

## 5. RSA Equations (Radiation Streaming Approximation)

When RSA is on, photon (and ur) perturbations are NOT integrated as ODEs.
Instead, they are computed algebraically from the metric at each step via
`perturb_rsa_delta_and_theta()` (line 10181).

### 5.1 Newtonian gauge (rsa_MD)

```
delta_g = -4 * phi
theta_g = 6 * phi'
```

### 5.2 Synchronous gauge (rsa_MD)

```
delta_g = (4/k^2) * (a'/a * h' - k^2 * eta)
theta_g = -(1/2) * h'
```

### 5.3 Reionization correction (rsa_MD_with_reio)

Additional terms proportional to `kappa'` (Thomson scattering rate during reionization):

**Newtonian gauge:**
```
delta_g += -(4/k^2) * kappa' * theta_b

theta_g += (3/k^2) * [ kappa'' * theta_b
                      + kappa' * (-a'/a * theta_b + c_b^2 * k^2 * delta_b + k^2 * phi) ]
```

**Synchronous gauge:**
```
delta_g += -(4/k^2) * kappa' * (theta_b + h'/2)

theta_g += (3/k^2) * [ kappa'' * (theta_b + h'/2)
                      + kappa' * (-a'/a * theta_b + c_b^2*k^2*delta_b - a'/a*h' + k^2*eta) ]
```

### 5.4 Ultra-relativistic neutrinos

Same formulas as photons (without reionization correction):
```
delta_ur = delta_g  (from metric)
theta_ur = theta_g  (from metric)
```

### 5.5 RSA contribution to Einstein equations

When RSA is on, the photon/ur contributions to the total stress-energy
are added via the RSA algebraic values (line 10287-10294):
```
delta_rho      += rho_g * rsa_delta_g
rho_p_theta    += (4/3)*rho_g * rsa_theta_g
(similarly for ur)
```

## 6. Full Hierarchy Equations (TCA off, RSA off)

These are the equations in `perturb_derivs()` starting at line 8919.

### 6.1 Photon temperature hierarchy

**Density** (l=0, line 8925):
```
d(delta_g)/dtau = -(4/3) * (theta_g + metric_continuity)
```

**Velocity** (l=1, line 8976-8979):
```
d(theta_g)/dtau = k^2 * (delta_g/4 - s2^2 * sigma_g)
                + metric_euler
                + kappa' * (theta_b - theta_g)
```

**Shear** (l=2, line 8982-8985):
```
d(sigma_g)/dtau = (1/2) * [ (8/15)*(theta_g + metric_shear)
                           - (3/5)*k*s_l[3]/s_l[2] * F_{g,3}
                           - kappa'*(2*sigma_g - (4/5)/s_l[2] * Pi) ]
```
where `Pi = (G_{g,0} + G_{g,2} + F_{g,2})/8 = (pol0_g + pol2_g + 2*s_l[2]*sigma_g)/8`.

**l=3** (line 8990-8992):
```
d(F_{g,3})/dtau = k/(2l+1) * [l*s_l[l]*2*s_l[2]*sigma_g - (l+1)*s_l[l+1]*F_{g,4}]
                - kappa' * F_{g,3}
```

**l > 3** (line 8995-8999):
```
d(F_{g,l})/dtau = k/(2l+1) * [l*s_l[l]*F_{g,l-1} - (l+1)*s_l[l+1]*F_{g,l+1}]
                - kappa' * F_{g,l}
```

**l = l_max** (line 9003-9006, truncation):
```
d(F_{g,lmax})/dtau = k * [s_l[lmax]*F_{g,lmax-1} - (1+lmax)*cotKgen*F_{g,lmax}]
                   - kappa' * F_{g,lmax}
```

### 6.2 Photon polarization hierarchy

**l=0** (line 9010-9012):
```
d(G_{g,0})/dtau = -k*G_{g,1} - kappa'*(G_{g,0} - 4*Pi)
```

**l=1** (line 9016-9018):
```
d(G_{g,1})/dtau = k/3*(G_{g,0} - 2*s_l[2]*G_{g,2}) - kappa'*G_{g,1}
```

**l=2** (line 9022-9024):
```
d(G_{g,2})/dtau = k/5*(2*s_l[2]*G_{g,1} - 3*s_l[3]*G_{g,3})
                - kappa'*(G_{g,2} - (4/5)*Pi)
```

### 6.3 Baryon equations (TCA off)

**Density** (line 8931):
```
d(delta_b)/dtau = -(theta_b + metric_continuity)
```

**Velocity** (line 8940-8944):
```
d(theta_b)/dtau = -a'/a * theta_b + metric_euler
                + k^2 * delta_p_b/rho_b
                + R * kappa' * (theta_g - theta_b)
```

### 6.4 CDM equations (synchronous gauge)

**Density** (line 9071):
```
d(delta_cdm)/dtau = -metric_continuity  (= -h'/2)
```
Velocity is zero by gauge choice.

### 6.5 GDM equations (gdm_class extension, line 9173-9216)

**Density:**
```
d(delta_gdm)/dtau = -(1+w)*( theta_gdm + metric_continuity )
                  + 3*a'/a * [ (w-ca2)*delta_gdm - Pi_nad ]
```
where `Pi_nad = (cs2 - ca2)*(delta_gdm + 3*a'/a*(1+w)*theta_gdm/k^2)`.

**Velocity:**
```
d(theta_gdm)/dtau = -(1-3*ca2)*a'/a * theta_gdm
                  + k^2/(1+w) * (ca2*delta_gdm + Pi_nad)
                  + metric_euler - s2^2*k^2*sigma_gdm
```

**Shear** (if dynamic):
```
d(sigma_gdm)/dtau = -3*a'/a * sigma_gdm
                  + (8/3)*cv2/(1+w) * (theta_gdm + metric_shear)
```

### 6.6 Ultra-relativistic neutrino hierarchy (RSA off, UFA off)

Same structure as photons but without Thomson scattering:
```
d(delta_ur)/dtau = -(4/3)*(theta_ur + metric_continuity)
d(theta_ur)/dtau = k^2*(delta_ur/4 - s2^2*sigma_ur) + metric_euler
d(sigma_ur)/dtau = (1/2)*(8/15)*(theta_ur + metric_shear) - (3/5)*k*s_l[3]/s_l[2]*F_{ur,3}
...
```

### 6.7 UFA equations (UFA on, RSA off)

When UFA is on, the ur hierarchy is truncated to {delta, theta, shear}.
Shear evolves via a simplified equation (three options):

**CLASS version** (line 9395-9399):
```
d(sigma_ur)/dtau = -(3/tau)*sigma_ur + (2/3)*(theta_ur + metric_ufa_class)
```
where `metric_ufa_class = h'/2` (synchronous) or `-6*phi'` (Newtonian).

## 7. State Vector Layout

The state vector `ppw->pv->y[]` is dynamically constructed at each approximation switch by `perturb_vector_init()` (line 3720+). The indices depend on which approximations are active. The total number of equations is `ppw->pv->pt_size`.

### 7.1 Scalar modes (TCA on)

```
y[index_pt_delta_g]     = delta_g       (photon density)
y[index_pt_theta_g]     = theta_g       (photon velocity)
  -- no shear, no l>=3, no polarization during TCA --
y[index_pt_delta_b]     = delta_b       (baryon density)
y[index_pt_theta_b]     = theta_b       (baryon velocity)
y[index_pt_delta_cdm]   = delta_cdm     (CDM density, if has_cdm)
y[index_pt_delta_gdm]   = delta_gdm     (GDM density, if has_gdm)
y[index_pt_theta_gdm]   = theta_gdm     (GDM velocity, if has_gdm)
y[index_pt_shear_gdm]   = sigma_gdm     (GDM shear, if dynamic_shear_gdm)
y[index_pt_delta_ur]    = delta_ur       (ur density)
y[index_pt_theta_ur]    = theta_ur       (ur velocity)
y[index_pt_shear_ur]    = sigma_ur       (ur shear)
y[index_pt_l3_ur]       = F_{ur,3}      (ur l=3)
y[..delta_ur+l..]       = F_{ur,l}      (ur l=4..l_max_ur, if ufa_off)
y[index_pt_psi0_ncdm1+...] = Psi_l(q)   (ncdm phase-space, if ncdmfa_off)
y[index_pt_eta]         = eta            (synchronous gauge metric)
  -- OR --
y[index_pt_phi]         = phi            (Newtonian gauge metric)
```

Typical count during TCA: ~10-25 equations (depending on species and ur hierarchy truncation).

### 7.2 Scalar modes (TCA off, RSA off) -- full hierarchy

Additional variables compared to TCA:
```
y[index_pt_shear_g]     = sigma_g       (photon shear, l=2)
y[index_pt_l3_g]        = F_{g,3}       (photon l=3)
y[index_pt_l3_g+1..lmax]= F_{g,4..lmax} (photon l=4..l_max_g)
y[index_pt_pol0_g]      = G_{g,0}       (polarization l=0)
y[index_pt_pol1_g]      = G_{g,1}       (polarization l=1)
y[index_pt_pol2_g]      = G_{g,2}       (polarization l=2)
y[index_pt_pol3_g]      = G_{g,3}       (polarization l=3)
y[index_pt_pol3_g+1..]  = G_{g,4..lmax} (polarization l=4..l_max_pol_g)
y[index_pt_perturbed_recombination_delta_temp]  (if perturbed recombination)
y[index_pt_perturbed_recombination_delta_chi]   (if perturbed recombination)
```

Typical count: ~50-70 equations.

### 7.3 Scalar modes (RSA on)

Photon variables {delta_g, theta_g, sigma_g, l>=3, polarization} are REMOVED from the ODE system. Similarly for ur if has_ur. The remaining equations are:
```
delta_b, theta_b, delta_cdm, (GDM if present), eta/phi
```

Plus ncdm phase-space if still integrated. Typical count: ~5-15 equations.

### 7.4 Default Boltzmann hierarchy truncation (precisions.h)

```
l_max_g       = 12   (photon temperature)
l_max_pol_g   = 10   (photon polarization)
l_max_ur      = 17   (ultra-relativistic neutrinos)
l_max_dr      = 17   (decay radiation)
l_max_idr     = 17   (interacting dark radiation)
l_max_ncdm    = 17   (non-cold dark matter)
l_max_g_ten   = 5    (photon temperature, tensor)
l_max_pol_g_ten = 5  (photon polarization, tensor)
```

## 8. NDF15 Evolver

File: `tools/evolver_ndf15.c` by Thomas Tram (2010).

### 8.1 Algorithm

Variable-order (1-5), adaptive-stepsize implicit method based on Numerical Differentiation Formulas (NDF/BDF). Reference: Shampine & Reichelt, "The MATLAB ODE Suite", SIAM J. Sci. Comput. 18(1), 1997.

Key features:
- **Implicit** (suitable for stiff problems like photon-baryon coupling)
- **Variable order k = 1..5** (BDF1 through BDF5)
- **Adaptive step size** with error control
- **Lazy Jacobian**: only recomputed when Newton iterations fail to converge
- **Sparse matrix support**: automatic detection of sparsity pattern
- **Interpolation**: Hermite-like interpolation for output at arbitrary times

### 8.2 Method constants (line 87-88)

```c
G[5] = {1, 3/2, 11/6, 25/12, 137/60}
alpha[5] = {-37/200, -1/9, -8.23e-2, -4.15e-2, 0}
```

### 8.3 Initialization

1. Compute initial Jacobian via numerical differencing (`numjac`)
2. Estimate initial step size from `||f(y)||` and `||J*f||`
3. Start at order k=1 (backward Euler)
4. LU decompose `(I - h*gamma*J)` for the implicit solve

### 8.4 Main stepping loop (line 300+)

Each step:
1. Predict: `y_pred = y + sum(dif[1:k])` using backward differences
2. Correct: simplified Newton iteration (max 4 iterations)
   - Solve: `(I - h*invGa[k]*J) * del = h*invGa[k]*f(y_new) - (psi + difkp1)`
   - Update: `y_new = y_pred + difkp1`
3. Error estimate: `err = ||difkp1 * invwt|| * erconst[k]`
4. If `err > rtol`: reduce step size, possibly reduce order
5. If successful: consider increasing order (k-1, k, k+1 comparison)

### 8.5 Step size selection

```
h_opt = h / max(0.1, 0.833 * (rtol/err)^(1/(k+1)))
h_min = minimum_variation  (from ppr->smallest_allowed_variation)
h_max = (t_final - t_0) / 10
```

Order change heuristic: after k+2 consecutive steps at same h and k, compare error estimates at orders k-1, k, k+1 and choose the one giving the largest allowed step.

### 8.6 Jacobian strategy

- Full Jacobian computed at initialization
- Recomputed ONLY when Newton iteration is "too slow" (converging slowly or diverging)
- Numerical Jacobian via column-wise finite differences
- Optionally stored in compressed sparse column (CSC) format if sufficiently sparse

### 8.7 Timescale callback

The evolver calls `perturb_timescale()` (line 6244) to determine the relevant physical timescale at each point. This is used to set the maximum step size as:

```
h_max_physical = perturb_integration_stepsize * timescale
               = 0.5 * min(tau_k, tau_h, tau_c)
```

where:
- `tau_k = 1/k` (mode timescale, included when RSA off or ncdm present)
- `tau_h = 1/(aH)` (Hubble timescale, always included)
- `tau_c = 1/kappa'` (Thomson scattering timescale, included when TCA off and kappa' != 0)

### 8.8 Integration tolerance

```
tol_perturb_integration = 1e-5  (relative tolerance for error control)
perturb_integration_stepsize = 0.5  (fraction of timescale for max step)
```

### 8.9 Statistics

The evolver tracks (stored in `stepstat[6]`):
```
[0] = successful steps
[1] = failed steps
[2] = total function evaluations
[3] = Jacobian computations
[4] = LU decompositions
[5] = linear system solves
```

## 9. Metric Equations

The metric perturbations are computed in `perturb_einstein()` (called from `perturb_derivs`).

### 9.1 Metric source terms in the ODEs

```
           Synchronous gauge        Newtonian gauge
           -----------------        ---------------
continuity:  h'/2                    -3*phi'
euler:       0                       k^2 * psi
shear:       k^2 * alpha             0
             = (h'+6*eta')/2
```

### 9.2 Gauge choice and CDM

In synchronous gauge, CDM velocity is zero (gauge choice), so CDM has only one equation (density). In Newtonian gauge, CDM has both density and velocity equations.

## 10. Summary: How Many Time Steps?

The number of steps is not fixed; it depends on k and is determined adaptively by the NDF15 evolver. Typical ranges:

- **Small k** (super-horizon modes, k ~ 10^{-4} Mpc^{-1}): ~200-500 steps total
- **Intermediate k** (BAO scales, k ~ 0.01-0.1 Mpc^{-1}): ~500-2000 steps
- **Large k** (small scales, k ~ 1 Mpc^{-1}): ~1000-5000 steps

The number of OUTPUT sample points (tau_sampling) is determined separately by the code and is typically ~5000-10000 points for the full range. The evolver interpolates to provide values at these output times.

## 11. Key Precision Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tol_perturb_integration` | 1e-5 | Relative tolerance for NDF15 |
| `perturb_integration_stepsize` | 0.5 | h_max / timescale ratio |
| `tol_tau_approx` | (small) | Bisection tolerance for switch times |
| `start_small_k_at_tau_c_over_tau_h` | 0.0015 | When to start integration |
| `start_large_k_at_tau_h_over_tau_k` | 0.07 | When to start integration |
| `tight_coupling_trigger_tau_c_over_tau_h` | 0.015 | TCA switch criterion |
| `tight_coupling_trigger_tau_c_over_tau_k` | 0.01 | TCA switch criterion |
| `radiation_streaming_trigger_tau_over_tau_k` | 45.0 | RSA switch criterion (k*tau > 45) |
| `ur_fluid_trigger_tau_over_tau_k` | 30.0 | UFA switch criterion (k*tau > 30) |
| `ncdm_fluid_trigger_tau_over_tau_k` | 31.0 | NCDMFA switch criterion |
| `l_max_g` | 12 | Photon temperature hierarchy depth |
| `l_max_pol_g` | 10 | Photon polarization hierarchy depth |
| `l_max_ur` | 17 | UR neutrino hierarchy depth |
| `tight_coupling_approximation` | compromise_CLASS | TCA implementation variant |
| `radiation_streaming_approximation` | rsa_MD_with_reio | RSA implementation variant |
