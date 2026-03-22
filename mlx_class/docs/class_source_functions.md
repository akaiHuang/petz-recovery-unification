# CLASS Source Function Formulas

Extracted from `gdm_class_public/source/perturbations.c` (10,506 lines).
All line references are to that file unless stated otherwise.

---

## Architecture Overview

CLASS computes **three** scalar temperature source terms (`S_t0`, `S_t1`, `S_t2`), one polarization source (`S_p`), and several metric potential sources. These source terms are then convolved with different radial (Bessel) functions in the line-of-sight (LOS) integration to produce transfer functions:

| Source | Radial Function | Flat-sky limit |
|--------|----------------|----------------|
| `S_t0` (monopole + SW + ISW + Doppler) | `Phi_l(x)` = `j_l(x)` | `j_l(k(tau_0-tau))` |
| `S_t1` (dipole)                         | `dPhi_l/dx` = `j_l'(x)` | `j_l'(k(tau_0-tau))` |
| `S_t2` (quadrupole/polarization)        | `(3 d^2Phi_l/dx^2 + Phi_l)/(2s_2)` | Related to `j_l''` |
| `S_p`  (E-mode polarization)            | `sqrt(3/8 * l(l-1)(l+1)(l+2))/s_2 * Phi_l/(k(tau_0-tau))^2` | `j_l/(x^2)` |

where `s_2 = sqrt(1 - 3K/k^2)` is a curvature factor (=1 for flat space).

The full temperature transfer function is:
```
Delta_T_l(k) = int d(tau) [ S_t0 * Phi_l  +  S_t1 * dPhi_l/dx  +  S_t2 * R_l^(2) ]
```

The `_set_source_` macro (defined in `perturbations.h` line 14):
```c
#define _set_source_(index) ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index][index_tau * ppt->k_size[index_md] + index_k]
```

---

## Newtonian Gauge Sources

### Simple form (commented out, lines 7507-7513)

This is the textbook form but is **not** used because of numerical cancellation:

```c
/* newtonian gauge: simplest form, not efficient numerically */
/*
  if (ppt->gauge == newtonian) {
  _set_source_(ppt->index_tp_t0) = exp(-kappa) * phi'  +  g * delta_g/4;
  _set_source_(ppt->index_tp_t1) = exp(-kappa) * k*psi  +  g * theta_b/k;
  _set_source_(ppt->index_tp_t2) = g * P;
  }
*/
```

In math:
- `S_t0 = e^{-kappa} * phi' + g * delta_g/4`
- `S_t1 = e^{-kappa} * k*psi + g * theta_b/k`
- `S_t2 = g * P`

### Efficient form (ACTUALLY USED, lines 7518-7528)

```c
if (ppt->gauge == newtonian) {
    _set_source_(ppt->index_tp_t0) =
      ppt->switch_sw * g * (delta_g/4. + psi)
      + switch_isw * (g * (phi - psi)
                      + exp_m_kappa * 2. * phi_prime)
      + ppt->switch_dop /k/k * (g * theta_b_dot
                                + g_dot * theta_b);

    _set_source_(ppt->index_tp_t1) =
      switch_isw * exp_m_kappa * k * (psi - phi);

    _set_source_(ppt->index_tp_t2) = ppt->switch_pol * g * P;
}
```

In math notation:

**S_t0** (Newtonian, efficient):
```
S_t0 = [SW]   g(tau) * (delta_g/4 + Psi)
     + [ISW]  g(tau) * (Phi - Psi)  +  e^{-kappa} * 2*Phi'
     + [Dop]  (1/k^2) * [ g * d(theta_b)/d(tau)  +  g' * theta_b ]
```

**S_t1** (Newtonian, efficient):
```
S_t1 = e^{-kappa} * k * (Psi - Phi)
```

**S_t2** (Newtonian, efficient):
```
S_t2 = g(tau) * P
```

**Key insight**: The "efficient" form bundles the ISW contribution differently. The simple form puts `phi'` entirely in `S_t0` and `k*psi` in `S_t1`. The efficient form uses integration by parts to move the Doppler term `g*theta_b/k` from `S_t1` into `S_t0` (as a `(g*theta_b_dot + g_dot*theta_b)/k^2` term), while `S_t1` keeps only the ISW-related `e^{-kappa}*k*(Psi-Phi)` piece. This avoids the large cancellation between the baryon velocity Doppler term and the `k*psi` gravitational redshift in `S_t1`.

---

## Synchronous Gauge Sources

### Simple form (commented out, lines 7532-7538)

```c
/* synchronous gauge: simplest form, not efficient numerically */
/*
  if (ppt->gauge == synchronous) {
  _set_source_(ppt->index_tp_t0) = - exp_m_kappa * h'/6  +  g/4 * delta_g;
  _set_source_(ppt->index_tp_t1) = g * theta_b / k;
  _set_source_(ppt->index_tp_t2) = exp_m_kappa * k^2 * (2/3) * s_l[2] * alpha  +  g * P;
  }
*/
```

In math:
- `S_t0 = -e^{-kappa} * h'/(6)  +  g * delta_g/4`
- `S_t1 = g * theta_b/k`
- `S_t2 = e^{-kappa} * (2/3) * k^2 * s_2 * alpha  +  g * P`

where `alpha = (h' + 6*eta') / (2*k^2)`.

### Efficient form (ACTUALLY USED, lines 7543-7563)

```c
if (ppt->gauge == synchronous) {

    _set_source_(ppt->index_tp_t0) =
      ppt->switch_sw * g * (delta_g/4. + alpha_prime)
      + switch_isw * (g * (eta - alpha_prime - 2*a'/a * alpha)
                      + exp_m_kappa * 2. * (eta_prime
                                             - a''/a * alpha
                                             - a'/a * alpha_prime))
      + ppt->switch_dop * (g * (theta_b_dot/k^2 + alpha_prime)
                           + g_dot * (theta_b/k^2 + alpha));

    _set_source_(ppt->index_tp_t1) =
      switch_isw * exp_m_kappa * k * (alpha_prime
                                      + 2. * a'/a * alpha
                                      - eta);

    _set_source_(ppt->index_tp_t2) =
      ppt->switch_pol * g * P;
}
```

In math notation:

**S_t0** (synchronous, efficient):
```
S_t0 = [SW]   g * (delta_g/4 + alpha')
     + [ISW]  g * (eta - alpha' - 2*H_conf*alpha)
              + e^{-kappa} * 2 * (eta' - H_conf'*alpha - H_conf*alpha')
     + [Dop]  g * (theta_b_dot/k^2 + alpha')
              + g' * (theta_b/k^2 + alpha)
```

where `H_conf = a'/a = aH` and `H_conf' = (a'/a)' = aH' + (aH)^2`.

**S_t1** (synchronous, efficient):
```
S_t1 = e^{-kappa} * k * (alpha' + 2*H_conf*alpha - eta)
```

**S_t2** (synchronous, efficient):
```
S_t2 = g * P
```

**Critical observation**: In the efficient synchronous form, the `alpha`-dependent terms in `S_t2` from the simple form have been **completely absorbed** into `S_t0` and `S_t1` via integration by parts. The efficient `S_t2` is **identical** in both gauges: just `g*P`.

---

## Polarization Source (all gauges)

From line 7567-7576:
```c
/* all gauges. Note that the correct formula for the E source
   should have a minus sign, as shown in Hu & White. We put a
   plus sign to comply with the 'historical convention'
   established in CMBFAST and CAMB. */

_set_source_(ppt->index_tp_p) = sqrt(6.) * g * P;
```

In math:
```
S_p = sqrt(6) * g(tau) * P
```

**Sign convention warning**: The physical sign should be negative (Hu & White), but CLASS uses `+` to match CMBFAST/CAMB historical convention.

---

## The Polarization Combination P

From lines 7478-7491:

**When radiation streaming approximation (RSA) is ON:**
```c
delta_g = ppw->rsa_delta_g;    // RSA approximation
P = 0.;                         // no polarization in RSA
```

**When tight-coupling approximation (TCA) is ON:**
```c
P = 5. * s_l[2] * tca_shear_g / 8.;
```
i.e., `P = (5/8) * s_2 * sigma_g^{TCA}`, where `sigma_g` is the photon shear.

**Full Boltzmann hierarchy (no approximation):**
```c
P = (pol0_g + pol2_g + 2. * s_l[2] * shear_g) / 8.;
```

In math:
```
P = (E_0 + E_2 + 2*s_2*sigma_g) / 8
```

where:
- `pol0_g` = `E_0` = monopole of E-mode polarization
- `pol2_g` = `E_2` = quadrupole of E-mode polarization
- `shear_g` = `sigma_g` = photon quadrupole (shear) `F_gamma,2`
- `s_l[2]` = `s_2 = sqrt(1 - 3K/k^2)` (curvature factor, = 1 for flat)

---

## Tensor (Gravitational Wave) Sources

From lines 7920-7953:

**P for tensors** (lines 7920-7938):
```c
// Full Boltzmann hierarchy:
P = -(1./10.*delta_g_tensor + 2./7.*shear_g_tensor
      + 3./70.*F_g4_tensor
      - 3./5.*pol0_g_tensor + 6./7.*pol2_g_tensor
      - 3./70.*E4_tensor) / sqrt(6.);

// TCA:
P = 2./5.*sqrt(6.) * gwdot / dkappa;
```

**Tensor temperature source** (line 7942):
```c
_set_source_(ppt->index_tp_t2) = - gwdot * exp_m_kappa + g * P;
```

**Tensor polarization source** (line 7953):
```c
_set_source_(ppt->index_tp_p) = sqrt(6.) * g * P;
```

---

## Metric Potential Sources (Bardeen potentials)

From lines 7599-7644:

### Bardeen potential Phi (= phi in Newtonian gauge)
```c
// Newtonian:
source_phi = y[index_pt_phi];
// Synchronous:
source_phi = eta - (a'/a) * alpha;
```

### Phi' (time derivative)
```c
// Newtonian:
source_phi_prime = dy[index_pt_phi];
// Synchronous:
source_phi_prime = eta_dot - (a'/a)' * alpha - (a'/a) * alpha';
```

### Phi + Psi (sum of Bardeen potentials)
```c
// Newtonian:
source_phi_plus_psi = phi + psi;
// Synchronous:
source_phi_plus_psi = eta + alpha';
```

### Psi (= psi in Newtonian gauge)
```c
// Newtonian:
source_psi = psi;
// Synchronous:
source_psi = (a'/a) * alpha + alpha';
```

---

## Gauge Transformation: How CLASS Computes alpha

### YES, CLASS computes alpha explicitly (line 6553-6554):
```c
/* alpha = (h'+6eta')/2k^2 */
ppw->pvecmetric[ppw->index_mt_alpha] =
    (ppw->pvecmetric[ppw->index_mt_h_prime]
     + 6.*ppw->pvecmetric[ppw->index_mt_eta_prime]) / 2. / k2;
```

### alpha' is computed from the fourth Einstein equation (lines 6583-6586):
```c
ppw->pvecmetric[ppw->index_mt_alpha_prime] =
    - 2. * (a'/a) * alpha
    + eta
    - 4.5 * (a^2/k^2) * (rho+p)_shear;
```

In math:
```
alpha' = -2*H_conf*alpha + eta - (9/2)*(a^2/k^2)*pi_tot
```

where `pi_tot` = total anisotropic stress `(rho+p)*sigma`.

### The chain: Einstein equations -> metric quantities

In synchronous gauge, from `perturb_einstein()` (lines 6517-6586):

1. **h'** from energy constraint (line 6520-6521):
```
h' = [k^2 * s_2^2 * eta + (3/2) * a^2 * delta_rho] / [(1/2) * H_conf]
```

2. **eta'** from momentum constraint (line 6545):
```
eta' = [(3/2) * a^2 * (rho+p)*theta + (1/2)*K*h'] / (k^2 * s_2^2)
```

3. **h''** from pressure equation (lines 6548-6551):
```
h'' = -2*H_conf*h' + 2*k^2*s_2^2*eta - 9*a^2*delta_p
```

4. **alpha** = `(h' + 6*eta') / (2*k^2)` (line 6554)

5. **alpha'** from anisotropic stress equation (lines 6583-6586):
```
alpha' = -2*H_conf*alpha + eta - (9/2)*(a^2/k^2)*(rho+p)*sigma
```

### Newtonian gauge metric quantities (lines 6461-6494):

1. **psi** (line 6491):
```
Psi = phi - (9/2) * (a^2/k^2) * (rho+p)*sigma
```

2. **phi'** (line 6494):
```
Phi' = -H_conf * Psi + (3/2) * (a^2/k^2) * (rho+p)*theta
```

(`phi` itself is evolved as a dynamical variable, not from the constraint.)

---

## Gauge-Invariant Combinations Used by CLASS

### The synchronous-to-Newtonian dictionary

From the output/debug section (lines 8205-8214):
```c
// In synchronous gauge:
alpha = pvecmetric[index_mt_alpha];
psi   = H_conf * alpha + alpha';          // Bardeen Psi
phi   = eta - H_conf * alpha;             // Bardeen Phi

// In newtonian gauge:
psi = pvecmetric[index_mt_psi];
phi = y[index_pt_phi];
```

These are the standard gauge transformation formulas:
```
Psi = (a'/a)*alpha + alpha'
Phi = eta - (a'/a)*alpha
```

### Gauge-invariant structure of the efficient sources

In the efficient synchronous sources, we can verify they are gauge-invariant by substituting:

**S_t0 ISW part** (sync):
```
g*(eta - alpha' - 2*H*alpha) + e^{-kappa}*2*(eta' - H'*alpha - H*alpha')
```
Using `Phi = eta - H*alpha` and `Psi = H*alpha + alpha'`:
```
= g*(Phi - Psi + H*alpha - alpha') + e^{-kappa}*2*Phi'
  ... after algebra this equals the Newtonian ISW form:
= g*(Phi - Psi) + e^{-kappa}*2*Phi'
```

**S_t1** (sync):
```
e^{-kappa}*k*(alpha' + 2*H*alpha - eta)
= e^{-kappa}*k*(Psi - Phi)     [exact match to Newtonian form]
```

This confirms the sources are **gauge-invariant** by construction.

### Doppler source handling

In the efficient sync form, the Doppler term is:
```
g*(theta_b_dot/k^2 + alpha') + g'*(theta_b/k^2 + alpha)
```

Note that `theta_b/k^2 + alpha` is NOT `v_b^{Newtonian}/k` exactly, because in synchronous gauge `theta_b^{sync} = theta_b^{Newt} - k^2*alpha`. So `theta_b^{sync}/k^2 + alpha = theta_b^{Newt}/k^2`, confirming gauge invariance.

The integration by parts that CLASS performs:
```
Original:  int g * theta_b/k * j_l'(x) d(tau)
IBP:     = - int [g'*theta_b + g*theta_b'] * j_l(x)/k^2 d(tau)    [boundary terms vanish]
```

This moves the Doppler contribution from the `S_t1 * j_l'` integral into the `S_t0 * j_l` integral, avoiding the cancellation between `g*theta_b/k` and `e^{-kappa}*k*psi` in `S_t1`.

---

## Key Variables Dictionary

### Thermodynamic quantities (`pvecthermo[]`)
| Variable | Symbol | Meaning |
|----------|--------|---------|
| `index_th_dkappa` | `d_kappa/d_tau` | Thomson scattering rate (1/Mpc) |
| `index_th_exp_m_kappa` | `e^{-kappa}` | Optical depth damping factor |
| `index_th_g` | `g(tau)` | Visibility function = `(d_kappa/d_tau) * e^{-kappa}` |
| `index_th_dg` | `g'(tau)` | Visibility function derivative `dg/d_tau` |

### Metric quantities (`pvecmetric[]`)
| Variable | Symbol (sync) | Symbol (Newt) | Meaning |
|----------|--------------|---------------|---------|
| `index_mt_h_prime` | `h'` | N/A | Trace of metric perturbation derivative |
| `index_mt_h_prime_prime` | `h''` | N/A | Second derivative |
| `index_mt_eta_prime` | `eta'` | N/A | Conformal time derivative of eta |
| `index_mt_alpha` | `alpha` | N/A | `(h'+6*eta')/(2*k^2)` |
| `index_mt_alpha_prime` | `alpha'` | N/A | Time derivative of alpha |
| `index_mt_psi` | N/A | `Psi` | Bardeen potential (Newtonian potential) |
| `index_mt_phi_prime` | N/A | `Phi'` | Time derivative of Bardeen Phi |

### Evolved perturbation variables (`y[]`)
| Variable | Meaning |
|----------|---------|
| `index_pt_delta_g` | Photon density contrast `delta_gamma` |
| `index_pt_theta_g` | Photon velocity divergence `theta_gamma` |
| `index_pt_shear_g` | Photon shear (quadrupole) `sigma_gamma = F_{gamma,2}` |
| `index_pt_pol0_g` | E-mode polarization monopole `E_0` |
| `index_pt_pol2_g` | E-mode polarization quadrupole `E_2` |
| `index_pt_delta_b` | Baryon density contrast `delta_b` |
| `index_pt_theta_b` | Baryon velocity divergence `theta_b` |
| `index_pt_eta` | Synchronous gauge metric perturbation `eta` |
| `index_pt_phi` | Newtonian gauge Bardeen potential `Phi` |

### Curvature factor
| Variable | Definition |
|----------|------------|
| `s_l[2]` = `s_2` | `sqrt(1 - 3*K/k^2)` where K = spatial curvature. Equals 1 for flat universe. |

### Control switches
| Variable | Meaning |
|----------|---------|
| `switch_sw` | Include Sachs-Wolfe? (0/1) |
| `switch_eisw` | Include early ISW? (0/1) |
| `switch_lisw` | Include late ISW? (0/1) |
| `switch_isw` | Effective ISW switch (derived from eisw/lisw) |
| `switch_dop` | Include Doppler? (0/1) |
| `switch_pol` | Include polarization contribution? (0/1) |

---

## Answers to Key Questions

### 1. Does CLASS compute alpha = (h'+6eta')/(2k^2) explicitly?

**YES.** Line 6554:
```c
ppw->pvecmetric[ppw->index_mt_alpha] =
    (ppw->pvecmetric[ppw->index_mt_h_prime]
     + 6.*ppw->pvecmetric[ppw->index_mt_eta_prime]) / 2. / k2;
```

CLASS also computes `alpha'` explicitly (line 6583-6586) from the fourth Einstein equation (anisotropic stress). Both `alpha` and `alpha'` are stored in `pvecmetric` and used extensively in the synchronous gauge source functions.

### 2. How does CLASS handle the Doppler source to avoid cancellation?

CLASS uses **integration by parts** to move the Doppler contribution from `S_t1` (convolved with `j_l'`) into `S_t0` (convolved with `j_l`).

**Simple form (UNUSED):**
- `S_t1 = g*theta_b/k + e^{-kappa}*k*psi` (large cancellation between the two terms)

**Efficient form (USED):**
- `S_t1 = e^{-kappa}*k*(Psi-Phi)` (smooth, no cancellation near recombination)
- The Doppler piece moves to `S_t0` as: `(1/k^2) * [g*theta_b' + g'*theta_b]` (Newt.) or `g*(theta_b'/k^2 + alpha') + g'*(theta_b/k^2 + alpha)` (sync.)

This is critical for numerical accuracy: the original `S_t1` has `g*theta_b/k ~ O(1)` and `e^{-kappa}*k*psi ~ O(1)` at recombination, whose difference is `O(0.01)`. Integration by parts eliminates this subtraction.

### 3. What gauge-invariant combinations does CLASS use?

In the **efficient** source functions, all terms are manifestly gauge-invariant. The key combinations are:

- **SW term**: `g * (delta_g/4 + Psi)` (Newtonian) or `g * (delta_g/4 + alpha')` (synchronous)
- **ISW term**: `e^{-kappa} * 2*Phi'` plus correction terms
- **Doppler**: After IBP, `d/d_tau[g*v_b^{Newt}]/k^2` where `v_b^{Newt} = theta_b^{sync}/k + k*alpha`
- **S_t1**: `e^{-kappa}*k*(Psi-Phi)` in both gauges (same expression, different variable names)
- **S_t2**: `g*P` (identical in both gauges, fully gauge-invariant)

The gauge invariance is verified by the substitutions:
```
Phi = eta - H_conf*alpha
Psi = H_conf*alpha + alpha'
```

### 4. Does CLASS apply any corrections or factors that aren't obvious?

**Yes, several:**

1. **The `s_l[2]` curvature factors** (= `sqrt(1-3K/k^2)`): These appear in the polarization combination `P`, in the simple-form `S_t2`, and in the Bessel radial functions for `SCALAR_TEMPERATURE_2` and `SCALAR_POLARISATION_E`. They vanish for flat universes but are critical for open/closed models.

2. **Radiation Streaming Approximation (RSA)**: When active, `delta_g` is replaced by an analytic approximation (lines 10252-10253 for sync gauge):
   ```
   delta_g^{RSA} = (4/k^2) * [H_conf*h' - k^2*eta]
   ```
   and `P = 0` (no polarization in RSA).

3. **Tight-Coupling Approximation (TCA)**: When active, the photon shear is approximated (line 6569):
   ```
   sigma_g^{TCA} = (16/45) * (1/d_kappa) * (theta_g + k^2*alpha)
   ```
   and `P = (5/8)*s_2*sigma_g^{TCA}`.

4. **GDM modifications**: The `gdm_class_public` version adds GDM contributions to `rho_plus_p_shear` (which feeds into `alpha'` via the Einstein equations) and to the algebraic or dynamic shear. This affects the metric perturbations but does NOT modify the CMB source function formulas themselves -- the GDM effects enter only through the changed metric perturbations `(h', eta', alpha, alpha')`.

5. **Reionization correction to RSA**: When `rsa_MD_with_reio` is active, additional Thomson-scattering correction terms are added to `delta_g^{RSA}` and `theta_g^{RSA}` (lines 10257-10270).

6. **Sign convention for E-mode**: CLASS uses a `+` sign for the E-mode source to match CMBFAST/CAMB, while the physical convention (Hu & White) uses `-`.

---

## Summary: What an MLX reimplementation needs

To reproduce CLASS's C_l to high accuracy, an MLX implementation must:

1. **Evolve** the coupled ODE system for `(delta_g, theta_g, shear_g, ..., delta_b, theta_b, eta/phi, ...)` with TCA and RSA approximation switching.

2. **Compute metric perturbations** from Einstein equations at each timestep:
   - Sync: `h' -> eta' -> h'' -> alpha -> (shear corrections) -> alpha'`
   - Newt: `psi -> phi'` (with `phi` evolved dynamically)

3. **Compute the three source functions** using the "efficient" form (with IBP), NOT the simple textbook form.

4. **Convolve** each source with the appropriate radial function:
   - `S_t0` with `j_l(x)`
   - `S_t1` with `j_l'(x)`
   - `S_t2` with `(3*j_l''(x) + j_l(x))/(2*s_2)` (flat: simplified)
   - `S_p` with `sqrt(3/8 * l(l-1)(l+1)(l+2)) * j_l(x) / (s_2 * x^2)`

5. **Include all approximation regimes**: TCA for early times, full Boltzmann for intermediate, RSA for late times. The transitions must be smooth.
