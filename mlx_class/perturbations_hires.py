"""
perturbations_hires.py — High-resolution photon Boltzmann hierarchy solver.

Extends NeutrinoBoltzmannSolver with a configurable photon multipole hierarchy
(l_gamma_max up to 100), replacing the tight-coupling approximation (TCA)
that only tracks Theta_0 and Theta_1.

PHYSICS:
  The photon Boltzmann hierarchy in conformal Newtonian gauge:

  Theta_0' = -k Theta_1 - Phi'
  Theta_1' = k/3 (Theta_0 + Psi) - kappa_dot (Theta_1 - v_b/3)
  Theta_2' = k/5 (2 Theta_1 - 3 Theta_3) - (9/10) kappa_dot Theta_2
  Theta_l' = k/(2l+1) [l Theta_{l-1} - (l+1) Theta_{l+1}]
             - kappa_dot Theta_l                                 (3 <= l < l_max)
  Theta_{l_max}' = k Theta_{l_max-1} - (l_max+1)/tau * Theta_{l_max}

STIFFNESS TREATMENT:
  The Thomson scattering terms are stiff (|kappa_dot| >> k at early times).
  We use Strang operator splitting:

    1. Half-step collision: apply Thomson damping analytically
       Theta_l -> Theta_l * exp(-|kappa_dot| * dt/2)  for l >= 2
       Theta_1, v_b: solve the coupled (Theta_1, v_b) collision system exactly
    2. Full-step streaming: explicit RK4 for the free-streaming + gravity terms
       (all kappa_dot terms removed from the RHS)
    3. Half-step collision: same as step 1

  This is second-order accurate and unconditionally stable for the stiff
  collision terms, regardless of how large |kappa_dot| is.

  The Phi equation remains IMEX (exponential integrator for k^2/calH stiffness).

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time

from .background import (
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    H0_Mpc as _H0_MPC,
)

from .perturbations_neutrino import (
    N_EFF, _f_nu, _OMEGA_GAMMA, _OMEGA_NU, L_NU_MAX,
)


# ============================================================================
# Variable layout
# ============================================================================

IDX_PHI = 0
IDX_DELTA_B = 1
IDX_V_B = 2
IDX_DELTA_C = 3
IDX_V_C = 4
# Theta_0 starts at index 5
IDX_THETA_START = 5


def idx_theta(l):
    """Index of Theta_l in the state vector."""
    return IDX_THETA_START + l


def idx_n_start(l_gamma_max):
    """Starting index of neutrino multipoles N_0."""
    return IDX_THETA_START + l_gamma_max + 1


def n_var_total(l_gamma_max, l_nu_max):
    """Total number of variables per k-mode."""
    return IDX_THETA_START + (l_gamma_max + 1) + (l_nu_max + 1)


# ============================================================================
# Phi equation: IMEX decomposition
# ============================================================================

def _phi_lambda(k_arr, calH):
    """Stiff eigenvalue: lambda = -(k^2/(3 calH) + calH)."""
    return -(k_arr * k_arr / (3.0 * calH) + calH)


def _phi_source(Theta_0, N_0, delta_b, delta_c, k_arr, calH, a,
                Omega_gamma, Omega_nu, Omega_b, Omega_c, H0_Mpc,
                N_2=None, Theta_2=None):
    """
    Non-stiff source F in the IMEX decomposition: Phi' = lambda*Phi + F.

    Includes anisotropic stress correction (Psi != Phi) when N_2/Theta_2
    are provided. This is essential for the radiation driving phase shift
    that places the first acoustic peak at l~220 instead of l~302.

    The full equation: Phi' = lambda*Phi + F_matter + F_aniso
    where F_aniso = -k^2/(3*calH) * (Psi - Phi)
    and   Psi - Phi = -12*H0^2*(Omega_nu*N_2 + Omega_gamma*Theta_2)/(a^2*k^2)
    """
    H02 = H0_Mpc * H0_Mpc
    S = (Omega_gamma / (a * a) * (4.0 * Theta_0)
         + Omega_nu / (a * a) * (4.0 * N_0)
         + Omega_b / a * delta_b
         + Omega_c / a * delta_c)
    F = -0.5 * H02 / calH * S

    # Anisotropic stress correction: F_aniso = -k^2/(3*calH) * (Psi - Phi)
    # Psi - Phi = -12 * H0^2 * (Omega_nu * N_2 + Omega_gamma * Theta_2) / (a^2 * k^2)
    # F_aniso = -k^2/(3*calH) * [-12*H0^2*(Om_nu*N2 + Om_gam*T2)/(a^2*k^2)]
    #         = 4*H0^2/(calH*a^2) * (Omega_nu*N_2 + Omega_gamma*Theta_2)
    if N_2 is not None:
        F_aniso = 4.0 * H02 / (calH * a * a) * Omega_nu * N_2
        if Theta_2 is not None:
            F_aniso = F_aniso + 4.0 * H02 / (calH * a * a) * Omega_gamma * Theta_2
        F = F + F_aniso

    return F


def _phi_int_factor(lam, dt):
    """(exp(lam*dt) - 1) / lam, with Taylor expansion for small |lam*dt|."""
    lam_dt = lam * dt
    exact = (mx.exp(lam_dt) - 1.0) / lam
    taylor = dt * (1.0 + 0.5 * lam_dt + lam_dt * lam_dt / 6.0)
    return mx.where(mx.abs(lam_dt) < 1e-4, taylor, exact)


# ============================================================================
# Thomson collision operator (applied analytically via operator splitting)
# ============================================================================

def _apply_collision_step(y, k_arr, R, kappa_dot_val, l_gamma_max, dt):
    """
    Apply the Thomson collision operator for a time step dt.

    The collision terms are:
      Theta_1' = -kappa_dot * (Theta_1 - v_b/3)
      Theta_2' = -(9/10) * kappa_dot * Theta_2
      Theta_l' = -kappa_dot * Theta_l      (l >= 3)
      v_b'     = kappa_dot / R * (v_b - 3*Theta_1)

    For l >= 2: pure exponential damping Theta_l -> Theta_l * exp(kappa_dot * dt)
    Note: kappa_dot is negative, so exp(kappa_dot * dt) < 1 (damping).

    For (Theta_1, v_b): coupled 2x2 linear system, solved exactly.
      d/dt [Theta_1]   [-kd    kd/3 ] [Theta_1]
           [v_b    ] = [3kd/R -kd/R ] [v_b    ]

    Eigenvalues: lambda = 0 and lambda = -kd(1 + 1/R)
    The zero eigenvalue corresponds to the conserved baryon-photon momentum:
      M = (1+R)*Theta_1 + (R/3)*v_b  (related to total momentum)
    Actually, the conserved quantity of the collision is:
      3*Theta_1 and v_b relax toward each other. In equilibrium, v_b = 3*Theta_1.

    Exact solution:
      Let S = v_b/3 - Theta_1 (slip). Then S' = -kd(1+1/R)*S
      => S(dt) = S(0) * exp(-kd(1+1/R)*dt)
      Theta_1(dt) = Theta_1(0) + [1/(1+R)] * S(0) * [1 - exp(-kd(1+1/R)*dt)]
      v_b(dt)     = v_b(0) - [3R/(1+R)] * S(0)/R * ... wait, let me redo this.

    Actually, note the conserved quantity: (4 rho_gamma) * Theta_1 + rho_b * v_b/3
    In our notation: Theta_1 + R * v_b / (3*(1+R))... Let me just use the matrix exponential.

    The coupled system: x' = A x where x = [Theta_1, v_b]^T
      A = [[-kd,   kd/3  ],
           [3kd/R, -kd/R ]]

    (Note kd < 0, so -kd > 0.)

    Eigenvalues of A: det(A - lambda I) = 0
      (-kd - lambda)(-kd/R - lambda) - (kd/3)(3kd/R) = 0
      lambda^2 + kd(1+1/R)lambda + kd^2/R - kd^2/R = 0
      lambda^2 + kd(1+1/R)lambda = 0
      lambda(lambda + kd(1+1/R)) = 0

    So lambda_1 = 0, lambda_2 = -kd(1+1/R) = |kd|(1+1/R) > 0 (stable decay).

    Eigenvectors:
      lambda_1 = 0: A v1 = 0 => -kd v1_1 + kd/3 v1_2 = 0 => v1_2 = 3 v1_1
        => v1 = [1, 3]^T (equilibrium: v_b = 3*Theta_1)
      lambda_2: v2 = [1, -3/R]^T (the slip mode)

    Solution:
      [Theta_1(t)] = c1 [1 ] + c2 exp(lambda_2 t) [1    ]
      [v_b(t)    ]      [3 ]                       [-3/R ]

    From IC: Theta_1(0) = c1 + c2, v_b(0) = 3c1 - 3c2/R
      c1 = (R*Theta_1(0) + v_b(0)/3) / (R + 1/1) ... let me compute properly.

    c1 + c2 = Theta_1(0)
    3c1 - 3c2/R = v_b(0)

    c1 = Theta_1(0) - c2
    3(Theta_1(0) - c2) - 3c2/R = v_b(0)
    3Theta_1(0) - 3c2 - 3c2/R = v_b(0)
    3Theta_1(0) - 3c2(1 + 1/R) = v_b(0)
    c2 = (3*Theta_1(0) - v_b(0)) / (3*(1 + 1/R))
       = R*(3*Theta_1(0) - v_b(0)) / (3*(R+1))

    c1 = Theta_1(0) - c2

    After time dt:
      Theta_1(dt) = c1 + c2 * exp(lambda_2 * dt)
      v_b(dt)     = 3*c1 - 3/R * c2 * exp(lambda_2 * dt)
    """
    # kappa_dot_val is negative. |kd| = -kappa_dot_val.
    kd = float(kappa_dot_val)  # negative
    _n_start = idx_n_start(l_gamma_max)

    # --- Theta_l for l >= 2: pure damping ---
    # Theta_2: damping rate is (9/10)*|kd|
    # Theta_l (l>=3): damping rate is |kd|
    # Note: kd < 0, so exp(kd * dt) = exp(-|kd|*dt) < 1 (damping)
    if l_gamma_max >= 2:
        exp_kd_dt_theta2 = mx.array(float(np.exp(0.9 * kd * dt)))  # (9/10)*kd < 0
        exp_kd_dt_higher = mx.array(float(np.exp(kd * dt)))  # kd < 0

        # Build the damping vector for all columns
        # We need to multiply each Theta_l column by its damping factor
        parts = []
        # Columns before Theta_2: no collision damping (Phi, delta_b, v_b, delta_c, v_c, Theta_0, Theta_1)
        parts.append(y[:, :idx_theta(2)])

        # Theta_2
        parts.append(y[:, idx_theta(2):idx_theta(2)+1] * exp_kd_dt_theta2)

        # Theta_3 ... Theta_{l_gamma_max}
        if l_gamma_max >= 3:
            parts.append(y[:, idx_theta(3):idx_theta(l_gamma_max)+1] * exp_kd_dt_higher)

        # Neutrinos: no collision
        parts.append(y[:, _n_start:])

        y = mx.concatenate(parts, axis=1)

    # --- (Theta_1, v_b): coupled system ---
    # lambda_2 = -kd * (1 + 1/R) [remember kd < 0, so lambda_2 > 0 => exp decays]
    # Actually lambda_2 = -kd*(1+1/R). Since kd<0, lambda_2 = |kd|*(1+1/R) > 0.
    # But in the exponent we have exp(lambda_2 * dt) which would GROW. That's wrong.
    # Let me recheck: the eigenvalue of the COLLISION ONLY system A is lambda_2 = -kd(1+1/R).
    # With kd < 0: lambda_2 = |kd|(1+1/R) > 0. This means the slip mode GROWS? No!
    #
    # Wait, the collision operator for Theta_1 is: Theta_1' = -kd*(Theta_1 - v_b/3)
    # Since kd < 0: -kd > 0, so if Theta_1 > v_b/3 (positive slip), then Theta_1 decreases.
    # That's DAMPING. The matrix element A[0,0] = -kd > 0... that means the diagonal is positive,
    # which gives growth. Something is wrong with my sign convention.
    #
    # Let me rewrite with |kd| = abs(kappa_dot):
    # Theta_1' = |kd| * (Theta_1 - v_b/3)  ... WAIT this would be growth!
    #
    # Actually kappa_dot is the scattering rate with sign convention:
    # In the equations: Theta_1' = k/3 (Theta_0 + Psi) + kappa_dot * (Theta_1 - v_b/3)
    # kappa_dot < 0 (our convention), so kappa_dot * (Theta_1 - v_b/3) < 0 when Theta_1 > v_b/3.
    # This IS damping. Good.
    #
    # So the collision-only system is:
    # Theta_1' = kd * (Theta_1 - v_b/3)     [kd < 0 => damping of slip]
    # v_b'     = -kd/R * (v_b - 3*Theta_1)  [kd < 0, -kd > 0, damping of slip]
    #
    # Wait, the v_b equation from the task description:
    # v_b' = -calH v_b + k Psi + kappa_dot / R * (v_b - 3 Theta_1)
    # Collision part: kappa_dot / R * (v_b - 3 Theta_1)
    # = kd / R * (v_b - 3*Theta_1)  [kd < 0]
    # If v_b > 3*Theta_1: kd/R * (positive) = negative => damping. Good.
    #
    # Matrix A (collision only):
    # A = [[kd,    -kd/3  ],   (from Theta_1' = kd*Theta_1 - kd*v_b/3)
    #      [-3kd/R, kd/R  ]]   (from v_b' = -3kd/R*Theta_1 + kd/R*v_b)
    #
    # Wait let me be more careful:
    # Theta_1' = kd*(Theta_1 - v_b/3) = kd*Theta_1 - kd/3*v_b
    # v_b' = kd/R*(v_b - 3*Theta_1) = -3kd/R*Theta_1 + kd/R*v_b
    #
    # A = [[kd,     -kd/3],
    #      [-3kd/R,  kd/R]]
    #
    # Eigenvalues: (kd - lambda)(kd/R - lambda) - (-kd/3)(-3kd/R) = 0
    # (kd-lambda)(kd/R-lambda) - kd^2/R = 0
    # kd^2/R - kd*lambda/R - kd*lambda + lambda^2 - kd^2/R = 0
    # lambda^2 - kd*lambda*(1+1/R) = 0
    # lambda*(lambda - kd*(1+1/R)) = 0
    # lambda_1 = 0, lambda_2 = kd*(1+1/R)
    #
    # Since kd < 0 and R > 0: lambda_2 = kd*(1+1/R) < 0 => DECAYING. Good!
    #
    # Eigenvectors:
    # lambda_1 = 0: (kd)*v1 + (-kd/3)*v2 = 0 => v2 = 3*v1 => [1, 3]
    # lambda_2: (kd - kd(1+1/R))*v1 + (-kd/3)*v2 = 0
    #           (-kd/R)*v1 - kd/3*v2 = 0 => v2 = -3/(R)*v1 => [1, -3/R]
    #
    # General solution:
    # [Theta_1(t)] = c1*[1] + c2*exp(lambda_2*t)*[1   ]
    # [v_b(t)    ]     [3]                        [-3/R]
    #
    # At t=0: c1 + c2 = Theta_1(0), 3c1 - 3c2/R = v_b(0)
    # => c2 = R*(3*Theta_1(0) - v_b(0)) / (3*(R+1))
    # => c1 = Theta_1(0) - c2

    R_val = float(R)
    R_val_safe = max(R_val, 1e-10)
    lambda_2 = kd * (1.0 + 1.0 / R_val_safe)  # negative (decaying)
    exp_lam2_dt = float(np.exp(lambda_2 * dt))

    Theta_1 = y[:, idx_theta(1)]
    v_b = y[:, IDX_V_B]

    slip = 3.0 * Theta_1 - v_b  # = 3*(Theta_1 - v_b/3)
    c2 = mx.array(float(R_val_safe / (3.0 * (R_val_safe + 1.0)))) * slip
    c1 = Theta_1 - c2

    exp_factor = mx.array(exp_lam2_dt)
    Theta_1_new = c1 + c2 * exp_factor
    v_b_new = 3.0 * c1 - 3.0 / mx.array(float(R_val_safe)) * c2 * exp_factor

    # Replace Theta_1 and v_b in state vector
    y = mx.concatenate([
        y[:, :IDX_V_B],
        v_b_new[:, None],
        y[:, IDX_V_B+1:idx_theta(1)],
        Theta_1_new[:, None],
        y[:, idx_theta(1)+1:]
    ], axis=1)

    return y


# ============================================================================
# Streaming + gravity derivative (NO collision terms)
# ============================================================================

def deriv_streaming(y, k_arr, calH, R, tau,
                    Omega_gamma, Omega_nu, Omega_b, Omega_c,
                    a, H0_Mpc, l_gamma_max, l_nu_max, tight_coupling):
    """
    RHS for the streaming + gravity part (collision terms removed).

    In the Strang splitting, this is the "non-stiff" part that we integrate
    with explicit RK4. The stiff Thomson collision terms are handled
    analytically in _apply_collision_step().
    """
    _n_start = idx_n_start(l_gamma_max)

    Phi = y[:, IDX_PHI]
    delta_b = y[:, IDX_DELTA_B]
    v_b = y[:, IDX_V_B]
    delta_c = y[:, IDX_DELTA_C]
    v_c = y[:, IDX_V_C]

    Theta_0 = y[:, idx_theta(0)]
    Theta_1 = y[:, idx_theta(1)]

    # Neutrino multipoles
    N_0 = y[:, _n_start]
    N_1 = y[:, _n_start + 1]
    N_2 = y[:, _n_start + 2] if l_nu_max >= 2 else mx.zeros_like(Phi)

    # Diagnose Psi from anisotropic stress
    H02 = H0_Mpc * H0_Mpc
    Psi = Phi - 12.0 * H02 * Omega_nu * N_2 / (a * a * k_arr * k_arr)
    Theta_2_val = mx.zeros_like(Phi)
    if l_gamma_max >= 2 and not tight_coupling:
        Theta_2_val = y[:, idx_theta(2)]
        Psi = Psi - 12.0 * H02 * Omega_gamma * Theta_2_val / (a * a * k_arr * k_arr)

    # Phi' with anisotropic stress correction (Psi != Phi)
    # Full: Phi' = -calH*Phi - k^2/(3*calH)*Psi - 0.5*H0^2/calH * S
    #            = lambda*Phi + F_matter + F_aniso
    # where F_aniso = -k^2/(3*calH)*(Psi - Phi)
    S = (Omega_gamma / (a * a) * (4.0 * Theta_0)
         + Omega_nu / (a * a) * (4.0 * N_0)
         + Omega_b / a * delta_b
         + Omega_c / a * delta_c)
    Phi_dot = (-calH * Phi
               - k_arr * k_arr * Psi / (3.0 * calH)
               - 0.5 * H02 / calH * S)

    # ---- Photon hierarchy (streaming only, no collision) ----
    dTheta = []

    # l = 0: Theta_0' = -k Theta_1 - Phi'  (no collision for monopole)
    dTheta.append(-k_arr * Theta_1 - Phi_dot)

    if tight_coupling:
        # TCA mode: Theta_1 via combined baryon-photon
        dTheta_1 = (k_arr / 3.0 * (Theta_0 + Psi) - calH * R * Theta_1) / (1.0 + R)
        dTheta.append(dTheta_1)
        for l in range(2, l_gamma_max + 1):
            dTheta.append(mx.zeros_like(Phi))
    else:
        # l = 1: streaming part only = k/3 (Theta_0 + Psi)
        # collision part (kd*(Theta_1 - v_b/3)) is handled by operator splitting
        dTheta.append(k_arr / 3.0 * (Theta_0 + Psi))

        if l_gamma_max >= 2:
            # l = 2: streaming part = k/5 (2 Theta_1 - 3 Theta_3)
            Theta_2 = y[:, idx_theta(2)]
            Theta_3 = y[:, idx_theta(3)] if l_gamma_max >= 3 else mx.zeros_like(Phi)
            dTheta.append(k_arr / 5.0 * (2.0 * Theta_1 - 3.0 * Theta_3))

            # l = 3 to l_max - 1: streaming only
            for l in range(3, l_gamma_max):
                Theta_lm1 = y[:, idx_theta(l - 1)]
                Theta_lp1 = y[:, idx_theta(l + 1)]
                dTheta.append(k_arr / (2.0 * l + 1.0) * (l * Theta_lm1 - (l + 1.0) * Theta_lp1))

            # l = l_max: truncation
            tau_safe = mx.maximum(mx.array(float(tau)), mx.array(1e-10))
            Theta_lm1 = y[:, idx_theta(l_gamma_max - 1)]
            Theta_lmax = y[:, idx_theta(l_gamma_max)]
            dTheta.append(k_arr * Theta_lm1 * l_gamma_max / (2.0 * l_gamma_max + 1.0)
                          - (l_gamma_max + 1.0) / tau_safe * Theta_lmax)

    # ---- Baryon equations ----
    if tight_coupling:
        d_v_b = 3.0 * dTheta[1]
        d_delta_b = -k_arr * v_b - 3.0 * Phi_dot
    else:
        # Streaming part of v_b: -calH v_b + k Psi
        # Collision part (kd/R*(v_b - 3 Theta_1)) handled by operator splitting
        d_v_b = -calH * v_b + k_arr * Psi
        d_delta_b = -k_arr * v_b - 3.0 * Phi_dot

    # ---- CDM ----
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot)
    d_v_c = -calH * v_c + k_arr * Psi

    # ---- Neutrino hierarchy (no collision) ----
    dN = []
    dN.append(-k_arr * N_1 - Phi_dot)
    dN.append(k_arr / 3.0 * (N_0 - 2.0 * N_2 + Psi))

    for l in range(2, l_nu_max):
        N_lm1 = y[:, _n_start + l - 1]
        N_lp1 = y[:, _n_start + l + 1]
        dN.append(k_arr / (2.0 * l + 1.0) * (l * N_lm1 - (l + 1.0) * N_lp1))

    tau_safe_nu = mx.maximum(mx.array(float(tau)), mx.array(1e-10))
    N_last = y[:, _n_start + l_nu_max]
    N_prev = y[:, _n_start + l_nu_max - 1]
    dN.append(k_arr * N_prev * l_nu_max / (2.0 * l_nu_max + 1.0)
              - (l_nu_max + 1.0) / tau_safe_nu * N_last)

    # ---- Stack ----
    parts = [
        Phi_dot[:, None],
        d_delta_b[:, None],
        d_v_b[:, None],
        d_delta_c[:, None],
        d_v_c[:, None],
    ]
    for l in range(l_gamma_max + 1):
        parts.append(dTheta[l][:, None])
    for l in range(l_nu_max + 1):
        parts.append(dN[l][:, None])

    return mx.concatenate(parts, axis=1)


# ============================================================================
# Column replacement helpers
# ============================================================================

def _set_phi(y, Phi_new):
    """Replace the Phi column (index 0) in the state vector."""
    return mx.concatenate([Phi_new[:, None], y[:, 1:]], axis=1)


# ============================================================================
# IMEX RK4 step (streaming only, no collision)
# ============================================================================

def _imex_rk4_streaming_step(y, k_arr, calH, R, tau_val,
                              Omega_gamma, Omega_nu, Omega_b, Omega_c,
                              a, H0_Mpc, l_gamma_max, l_nu_max, dtau,
                              tight_coupling):
    """
    IMEX RK4 step for the STREAMING part only.
    Phi uses exponential integrator; everything else uses explicit RK4.
    Thomson collision terms are NOT included (handled by operator splitting).
    """
    _n_start = idx_n_start(l_gamma_max)
    lam = _phi_lambda(k_arr, calH)
    exp_lam_dt = mx.exp(lam * dtau)
    exp_lam_half = mx.exp(lam * 0.5 * dtau)

    def get_F(state):
        # Extract anisotropic stress quadrupoles for Psi != Phi correction
        _N_2 = state[:, _n_start + 2] if l_nu_max >= 2 else None
        _Theta_2 = (state[:, idx_theta(2)]
                    if l_gamma_max >= 2 and not tight_coupling else None)
        return _phi_source(
            state[:, idx_theta(0)], state[:, _n_start],
            state[:, IDX_DELTA_B], state[:, IDX_DELTA_C],
            k_arr, calH, a, Omega_gamma, Omega_nu, Omega_b, Omega_c, H0_Mpc,
            N_2=_N_2, Theta_2=_Theta_2)

    def rhs(state):
        return deriv_streaming(state, k_arr, calH, R, tau_val,
                               Omega_gamma, Omega_nu, Omega_b, Omega_c,
                               a, H0_Mpc, l_gamma_max, l_nu_max, tight_coupling)

    # Stage 1
    Phi_1 = y[:, IDX_PHI]
    F1 = get_F(y)
    k1 = rhs(y)

    # Stage 2
    Phi_2 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F1
    y2 = y + 0.5 * dtau * k1
    y2 = _set_phi(y2, Phi_2)
    F2 = get_F(y2)
    k2 = rhs(y2)

    # Stage 3
    Phi_3 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F2
    y3 = y + 0.5 * dtau * k2
    y3 = _set_phi(y3, Phi_3)
    F3 = get_F(y3)
    k3 = rhs(y3)

    # Stage 4
    Phi_4 = exp_lam_dt * Phi_1 + _phi_int_factor(lam, dtau) * F3
    y4 = y + dtau * k3
    y4 = _set_phi(y4, Phi_4)
    F4 = get_F(y4)
    k4 = rhs(y4)

    # Combine RK4
    y_new = y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # Phi: Simpson's rule exponential integrator
    F_mid = 0.5 * (F2 + F3)
    Phi_new = (exp_lam_dt * Phi_1
               + (dtau / 6.0) * (exp_lam_dt * F1
                                  + 4.0 * exp_lam_half * F_mid
                                  + F4))
    y_new = _set_phi(y_new, Phi_new)

    return y_new


# ============================================================================
# Full Strang splitting step
# ============================================================================

def strang_split_step(y, k_arr, calH, R, tau_val, kappa_dot_val,
                      Omega_gamma, Omega_nu, Omega_b, Omega_c,
                      a, H0_Mpc, l_gamma_max, l_nu_max, dtau,
                      tight_coupling):
    """
    Strang operator splitting step (second-order):
      1. Half-step collision (dt/2)
      2. Full-step streaming (dt)
      3. Half-step collision (dt/2)

    In TCA mode, collision is not applied separately (it's embedded in the
    TCA equations for Theta_1 via the (1+R) denominator).
    """
    if tight_coupling:
        # In TCA, collision is already part of the equation
        y = _imex_rk4_streaming_step(
            y, k_arr, calH, R, tau_val,
            Omega_gamma, Omega_nu, Omega_b, Omega_c,
            a, H0_Mpc, l_gamma_max, l_nu_max, dtau,
            tight_coupling=True)
        # Enforce v_b = 3 * Theta_1
        Theta_1_new = y[:, idx_theta(1)]
        v_b_locked = 3.0 * Theta_1_new
        y = mx.concatenate([
            y[:, :IDX_V_B],
            v_b_locked[:, None],
            y[:, IDX_V_B + 1:]
        ], axis=1)
    else:
        # Step 1: half collision
        y = _apply_collision_step(y, k_arr, R, kappa_dot_val, l_gamma_max, 0.5 * dtau)
        # Step 2: full streaming
        y = _imex_rk4_streaming_step(
            y, k_arr, calH, R, tau_val,
            Omega_gamma, Omega_nu, Omega_b, Omega_c,
            a, H0_Mpc, l_gamma_max, l_nu_max, dtau,
            tight_coupling=False)
        # Step 3: half collision
        y = _apply_collision_step(y, k_arr, R, kappa_dot_val, l_gamma_max, 0.5 * dtau)

    return y


# ============================================================================
# Tau grid
# ============================================================================

def build_tau_grid_hires(bg, k_max, N_early=300, N_late=800):
    """Build integration grid extending to tau_rec."""
    tau_rec = bg.tau_rec
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * tau_rec)

    tau_early = np.geomspace(tau_init, tau_early_end, N_early)

    cs_approx = 1.0 / np.sqrt(3.0)
    dtau_acoustic = 2.0 * np.pi / (k_max * cs_approx) / 25.0
    N_needed = int((tau_rec - tau_early_end) / dtau_acoustic) + 10
    N_late_actual = max(N_late, N_needed)

    tau_late = np.linspace(tau_early_end, tau_rec, N_late_actual)
    tau_grid = np.unique(np.concatenate([tau_early, tau_late]))
    return tau_grid


# ============================================================================
# Initial conditions
# ============================================================================

def adiabatic_ic_hires(k_arr_np, bg, l_gamma_max, l_nu_max=L_NU_MAX):
    """
    Adiabatic initial conditions for the high-resolution system.

    In the tight-coupling limit (early times):
      Phi     = 1.0
      delta_b = -3/2 Phi
      v_b     = k tau / 6 (= 3 * Theta_1)
      delta_c = -3/2 Phi
      v_c     = k tau / 6
      Theta_0 = -Phi/2
      Theta_1 = k tau / 18
      Theta_l = 0  for l >= 2
      N_0     = -Phi/2
      N_1     = k tau / 18
      N_2     = (k tau)^2 / 60
    """
    N_k = len(k_arr_np)
    tau_init = bg.tau_grid[1]
    nvar = n_var_total(l_gamma_max, l_nu_max)
    _n_start = idx_n_start(l_gamma_max)

    y0 = np.zeros((N_k, nvar), dtype=np.float32)

    y0[:, IDX_PHI] = 1.0
    y0[:, IDX_DELTA_B] = -1.5
    y0[:, IDX_V_B] = k_arr_np * tau_init / 6.0
    y0[:, IDX_DELTA_C] = -1.5
    y0[:, IDX_V_C] = k_arr_np * tau_init / 6.0

    y0[:, idx_theta(0)] = -0.5
    y0[:, idx_theta(1)] = k_arr_np * tau_init / 18.0

    # Neutrinos
    y0[:, _n_start + 0] = -0.5
    y0[:, _n_start + 1] = k_arr_np * tau_init / 18.0
    if l_nu_max >= 2:
        y0[:, _n_start + 2] = (k_arr_np * tau_init) ** 2 / 60.0
    for l in range(3, min(l_nu_max + 1, 6)):
        prod_val = 1.0
        for j in range(1, l + 1):
            prod_val *= (2 * j + 1)
        y0[:, _n_start + l] = (k_arr_np * tau_init) ** l / prod_val

    return mx.array(y0)


# ============================================================================
# Solver class
# ============================================================================

class HiResBoltzmannSolver:
    """
    High-resolution photon Boltzmann hierarchy solver with Strang splitting.

    Key features:
    1. Full photon hierarchy with l_gamma_max up to 100
    2. Strang operator splitting: collision (implicit/analytic) + streaming (RK4)
    3. TCA -> full hierarchy switch when |kappa_dot| drops below threshold
    4. Photon anisotropic stress (Theta_2) contribution to Psi
    5. Separate v_b with exact baryon-photon momentum exchange
    """

    def __init__(self, bg, k_arr_Mpc, l_gamma_max=25, l_nu_max=L_NU_MAX,
                 tca_threshold=50.0):
        """
        Parameters
        ----------
        bg : Background
            Solved background cosmology.
        k_arr_Mpc : array
            Wavenumber array in Mpc^-1.
        l_gamma_max : int
            Maximum photon multipole (default 25).
        l_nu_max : int
            Maximum neutrino multipole (default 20).
        tca_threshold : float
            Switch from TCA to full hierarchy when |kappa_dot| / k_max < threshold.
        """
        self.bg = bg
        self.k_arr_np = np.asarray(k_arr_Mpc, dtype=np.float32)
        self.N_k = len(self.k_arr_np)
        self.k_arr = mx.array(self.k_arr_np)
        self.l_gamma_max = l_gamma_max
        self.l_nu_max = l_nu_max
        self.tca_threshold = tca_threshold
        self.n_var = n_var_total(l_gamma_max, l_nu_max)

        k_max = float(self.k_arr_np[-1])
        self.tau_grid = build_tau_grid_hires(bg, k_max)

        self.Omega_gamma = _OMEGA_GAMMA
        self.Omega_nu = _OMEGA_NU
        self.Omega_b = _OMEGA_B
        self.Omega_c = float(bg.Omega_cdm)

        mem_kb = self.N_k * self.n_var * 4 / 1024
        print(f"[HiRes] N_k={self.N_k}, l_gamma_max={l_gamma_max}, "
              f"l_nu_max={l_nu_max}, N_var={self.n_var}")
        print(f"[HiRes] Memory: {mem_kb:.1f} KB ({self.N_k} x {self.n_var} x 4 bytes)")
        print(f"[HiRes] Grid: {len(self.tau_grid)} steps "
              f"[{self.tau_grid[0]:.2e}, {self.tau_grid[-1]:.1f}] Mpc "
              f"(tau_rec={bg.tau_rec:.1f})")
        print(f"[HiRes] TCA threshold: |kappa_dot|/k_max < {tca_threshold}")

    def solve(self):
        t0 = time.time()
        bg = self.bg
        k_arr = self.k_arr
        tau_grid = self.tau_grid

        y = adiabatic_ic_hires(self.k_arr_np, bg, self.l_gamma_max, self.l_nu_max)

        # Precompute background quantities
        calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
        R_np = bg.R_at_tau(tau_grid).astype(np.float32)
        a_np = bg.a_at_tau(tau_grid).astype(np.float32)
        kappa_dot_np = bg.kappa_dot_at_tau(tau_grid).astype(np.float32)

        calH_arr = mx.array(calH_np)
        R_arr = mx.array(R_np)
        a_arr = mx.array(a_np)

        _Og = self.Omega_gamma
        _On = self.Omega_nu
        _Ob = self.Omega_b
        _Oc = self.Omega_c
        _H0 = _H0_MPC
        _lgmax = self.l_gamma_max
        _lnmax = self.l_nu_max
        k_max = float(self.k_arr_np[-1])

        # Determine TCA -> full hierarchy switch point
        abs_kappa_dot = np.abs(kappa_dot_np)
        tca_switch_idx = None
        for i in range(len(tau_grid)):
            if abs_kappa_dot[i] / k_max < self.tca_threshold:
                tca_switch_idx = i
                break

        if tca_switch_idx is None:
            tca_switch_idx = len(tau_grid) - 1
            print(f"[HiRes] WARNING: TCA never switches off!")

        tau_switch = tau_grid[tca_switch_idx]
        z_switch = 1.0 / float(bg.a_at_tau(np.array([tau_switch]))[0]) - 1.0
        print(f"[HiRes] TCA -> full hierarchy at tau = {tau_switch:.1f} Mpc "
              f"(z = {z_switch:.0f}, step {tca_switch_idx}/{len(tau_grid)})")

        # ---- TCA phase ----
        print(f"[HiRes] Integrating (TCA phase)...")
        for i in range(tca_switch_idx):
            dt = float(tau_grid[i + 1] - tau_grid[i])
            y = strang_split_step(
                y, k_arr,
                calH_arr[i], R_arr[i], tau_grid[i], float(kappa_dot_np[i]),
                _Og, _On, _Ob, _Oc,
                a_arr[i], _H0, _lgmax, _lnmax, dt,
                tight_coupling=True)
            if i % 200 == 0:
                mx.eval(y)

        mx.eval(y)
        t_tca = time.time() - t0
        print(f"[HiRes] TCA phase done in {t_tca:.2f}s ({tca_switch_idx} steps)")

        # Seed Theta_2 from quasi-static equilibrium:
        # At equilibrium: Theta_2 ~ (2k/(9|kd|/10)) * Theta_1 * (5/(2*5+1)) ... simplified:
        # Theta_2 ~ (2/9) * (k * Theta_1) / |kd| * (5/(2*2+1)) = (2k*Theta_1)/(9*|kd|) * (5/5)
        # Actually from the equation: 0 = k/5 (2 Theta_1) - (9/10) |kd| Theta_2
        # => Theta_2 = (10/9) * (k/(5*|kd|)) * 2 * Theta_1 = (4/9) * k * Theta_1 / |kd|
        kd_switch = abs_kappa_dot[tca_switch_idx]
        if kd_switch > 1e-10:
            Theta_1_switch = y[:, idx_theta(1)]
            Theta_2_seed = (4.0 / 9.0) * k_arr * Theta_1_switch / kd_switch
            _idx_t2 = idx_theta(2)
            y = mx.concatenate([
                y[:, :_idx_t2],
                Theta_2_seed[:, None],
                y[:, _idx_t2 + 1:]
            ], axis=1)
            mx.eval(y)

        # ---- Full hierarchy phase ----
        print(f"[HiRes] Integrating (full hierarchy phase)...")
        t_full_start = time.time()
        for i in range(tca_switch_idx, len(tau_grid) - 1):
            dt = float(tau_grid[i + 1] - tau_grid[i])
            y = strang_split_step(
                y, k_arr,
                calH_arr[i], R_arr[i], tau_grid[i], float(kappa_dot_np[i]),
                _Og, _On, _Ob, _Oc,
                a_arr[i], _H0, _lgmax, _lnmax, dt,
                tight_coupling=False)
            if i % 200 == 0:
                mx.eval(y)

        mx.eval(y)
        t_full = time.time() - t_full_start
        t_total = time.time() - t0
        print(f"[HiRes] Full hierarchy phase done in {t_full:.2f}s "
              f"({len(tau_grid) - 1 - tca_switch_idx} steps)")
        print(f"[HiRes] Total integration: {t_total:.2f}s")

        return HiResResult(y=y, k_arr=self.k_arr_np, bg=bg,
                           l_gamma_max=self.l_gamma_max,
                           l_nu_max=self.l_nu_max,
                           Omega_gamma=_Og, Omega_nu=_On, H0=_H0)


class HiResResult:
    """Result container for the high-resolution solver."""

    def __init__(self, y, k_arr, bg, l_gamma_max, l_nu_max,
                 Omega_gamma, Omega_nu, H0):
        self.y = y
        self.k_arr = k_arr
        self.bg = bg
        self.l_gamma_max = l_gamma_max
        self.l_nu_max = l_nu_max
        self.Omega_gamma = Omega_gamma
        self.Omega_nu = Omega_nu
        self.H0 = H0
        self.N_k = len(k_arr)

    def source_at_recombination(self):
        """
        Extract source functions at tau_rec with Silk damping.

        Returns (Theta_0, Phi, Psi, v_b, Theta_2).
        """
        y_np = np.array(self.y)
        _n_start = idx_n_start(self.l_gamma_max)

        Phi = y_np[:, IDX_PHI]
        Theta_0 = y_np[:, idx_theta(0)]
        v_b = y_np[:, IDX_V_B]
        Theta_2 = y_np[:, idx_theta(2)] if self.l_gamma_max >= 2 else np.zeros_like(Phi)

        N_2 = y_np[:, _n_start + 2] if self.l_nu_max >= 2 else np.zeros_like(Phi)

        # Diagnose Psi
        H02 = self.H0 ** 2
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])
        Psi = (Phi
               - 12.0 * H02 * self.Omega_nu * N_2 / (a_rec ** 2 * self.k_arr ** 2)
               - 12.0 * H02 * self.Omega_gamma * Theta_2 / (a_rec ** 2 * self.k_arr ** 2))

        silk = np.exp(-(self.k_arr / self.bg.k_D) ** 2)
        return (Theta_0 * silk, Phi * silk, Psi * silk,
                v_b * silk, Theta_2 * silk)

    def sw_source(self):
        """Sachs-Wolfe source: Theta_0 + Psi (silk-damped)."""
        Theta_0, Phi, Psi, v_b, Theta_2 = self.source_at_recombination()
        return Theta_0 + Psi

    def photon_multipoles(self):
        """Return all photon multipoles at tau_rec (no damping)."""
        y_np = np.array(self.y)
        multipoles = {}
        for l in range(self.l_gamma_max + 1):
            multipoles[l] = y_np[:, idx_theta(l)]
        return multipoles

    def phi_psi_ratio(self):
        """Psi/Phi at tau_rec (no Silk damping)."""
        y_np = np.array(self.y)
        _n_start = idx_n_start(self.l_gamma_max)
        Phi = y_np[:, IDX_PHI]
        N_2 = y_np[:, _n_start + 2] if self.l_nu_max >= 2 else np.zeros_like(Phi)
        Theta_2 = y_np[:, idx_theta(2)] if self.l_gamma_max >= 2 else np.zeros_like(Phi)

        H02 = self.H0 ** 2
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])
        Psi = (Phi
               - 12.0 * H02 * self.Omega_nu * N_2 / (a_rec ** 2 * self.k_arr ** 2)
               - 12.0 * H02 * self.Omega_gamma * Theta_2 / (a_rec ** 2 * self.k_arr ** 2))
        mask = np.abs(Phi) > 1e-10
        return np.where(mask, Psi / Phi, 1.0)


# ============================================================================
# Convergence test
# ============================================================================

def _test_hires_convergence():
    """
    Test the high-resolution photon hierarchy solver:
    1. Stability for l_gamma_max = 25, 50, 100
    2. C_l convergence
    3. Compare with TCA-only solver
    4. Performance benchmark
    """
    print("=" * 70)
    print("TEST: High-resolution photon Boltzmann hierarchy (Strang splitting)")
    print("=" * 70)

    from .background import Background
    from .spectra import compute_cl

    # Background
    print("\n--- Step 1: Background ---")
    bg = Background(khronon=False)
    bg.solve()

    # k-grid
    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

    # ell values for C_l
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1), np.arange(30, 100, 2),
        np.arange(100, 500, 4), np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)

    # ================================================================
    # Run for each l_gamma_max
    # ================================================================
    results = {}
    Dl_dict = {}
    timings = {}

    for lgmax in [25, 50, 100]:
        print(f"\n{'='*60}")
        print(f"  l_gamma_max = {lgmax}")
        print(f"{'='*60}")

        solver = HiResBoltzmannSolver(bg, k_arr, l_gamma_max=lgmax)
        t0 = time.time()
        result = solver.solve()
        t_solve = time.time() - t0
        timings[lgmax] = t_solve

        # Stability check
        Theta_0, Phi, Psi, v_b, Theta_2 = result.source_at_recombination()
        has_nan = any(np.any(np.isnan(x)) for x in [Theta_0, Phi, Psi])
        has_inf = any(np.any(np.isinf(x)) for x in [Theta_0, Phi, Psi])
        print(f"  Stability: NaN={has_nan}, Inf={has_inf}")

        if has_nan or has_inf:
            print(f"  FAIL: Unstable for l_gamma_max={lgmax}!")
            for name, arr in [("Theta_0", Theta_0), ("Phi", Phi), ("Psi", Psi)]:
                bad = np.where(np.isnan(arr) | np.isinf(arr))[0]
                if len(bad) > 0:
                    print(f"    {name}: first bad at k={k_arr[bad[0]]:.4f}")
            results[lgmax] = None
            continue

        results[lgmax] = result

        # Compute C_l (with Doppler term)
        sw = result.sw_source()
        print(f"  SW source range: [{np.min(sw):.4f}, {np.max(sw):.4f}]")
        _, _, Dl = compute_cl(sw, k_arr, ell_values, bg.D_A, source_Dop=v_b)
        Dl_dict[lgmax] = Dl

        # Multipole spectrum
        multipoles = result.photon_multipoles()
        rms_spectrum = []
        for l in [0, 1, 2, 5, 10, 20, min(50, lgmax), lgmax]:
            if l <= lgmax:
                rms = np.sqrt(np.mean(multipoles[l]**2))
                rms_spectrum.append((l, rms))
        print(f"  Multipole RMS at tau_rec:")
        for l, rms in rms_spectrum:
            print(f"    l={l:3d}: {rms:.6e}")

    # ================================================================
    # TCA reference
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  Reference: TCA (perturbations_neutrino.py)")
    print(f"{'='*60}")
    from .perturbations_neutrino import NeutrinoBoltzmannSolver
    solver_tca = NeutrinoBoltzmannSolver(bg, k_arr, l_nu_max=L_NU_MAX)
    t0 = time.time()
    result_tca = solver_tca.solve()
    t_tca = time.time() - t0
    timings['TCA'] = t_tca

    Theta_0_tca, Phi_tca, Psi_tca, v_b_tca, _, _ = result_tca.source_at_recombination()
    sw_tca = Theta_0_tca + Psi_tca
    _, _, Dl_tca = compute_cl(sw_tca, k_arr, ell_values, bg.D_A, source_Dop=v_b_tca)
    Dl_dict['TCA'] = Dl_tca

    # ================================================================
    # Convergence analysis
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  CONVERGENCE ANALYSIS")
    print(f"{'='*60}")

    from scipy.interpolate import interp1d
    from scipy.signal import find_peaks as _find_peaks

    l_full = np.arange(2, 2501)
    peak_data = {}
    for label, Dl in Dl_dict.items():
        if Dl is None:
            continue
        f = interp1d(ell_values, Dl, kind='cubic', fill_value='extrapolate')
        Dl_interp = np.maximum(f(l_full), 0.0)
        peaks, _ = _find_peaks(Dl_interp, distance=80, prominence=20)
        if len(peaks) >= 3:
            peak_ls = l_full[peaks[:7]]
            peak_Ds = Dl_interp[peaks[:7]]
            peak_data[label] = (peak_ls, peak_Ds)
            print(f"\n  {label}:")
            print(f"    Peaks at l = {peak_ls}")
            print(f"    Heights = {[f'{d:.0f}' for d in peak_Ds]} uK^2")
            if len(peaks) >= 2:
                print(f"    1st/2nd ratio: {peak_Ds[0]/peak_Ds[1]:.3f}")

    # C_l differences
    if 25 in Dl_dict and 100 in Dl_dict and Dl_dict[25] is not None and Dl_dict[100] is not None:
        Dl_25 = Dl_dict[25]
        Dl_100 = Dl_dict[100]
        mask = Dl_100 > 10
        rel_diff = np.abs(Dl_25[mask] - Dl_100[mask]) / Dl_100[mask]
        print(f"\n  C_l relative difference (l_max=25 vs l_max=100):")
        print(f"    Max:  {np.max(rel_diff)*100:.2f}%")
        print(f"    Mean: {np.mean(rel_diff)*100:.2f}%")
        print(f"    RMS:  {np.sqrt(np.mean(rel_diff**2))*100:.2f}%")

        for l_lo, l_hi, name in [(2, 50, "SW plateau"),
                                  (100, 350, "1st peak"),
                                  (350, 700, "2nd peak"),
                                  (700, 1200, "3rd-4th peaks"),
                                  (1200, 2500, "damping tail")]:
            band_mask = (ell_values >= l_lo) & (ell_values < l_hi) & (Dl_100 > 10)
            if np.sum(band_mask) > 0:
                band_diff = np.mean(np.abs(Dl_25[band_mask] - Dl_100[band_mask])
                                    / Dl_100[band_mask]) * 100
                print(f"    {name} (l={l_lo}-{l_hi}): {band_diff:.2f}%")

    if 50 in Dl_dict and 100 in Dl_dict and Dl_dict[50] is not None and Dl_dict[100] is not None:
        Dl_50 = Dl_dict[50]
        Dl_100 = Dl_dict[100]
        mask = Dl_100 > 10
        rel_diff_50 = np.abs(Dl_50[mask] - Dl_100[mask]) / Dl_100[mask]
        print(f"\n  C_l relative difference (l_max=50 vs l_max=100):")
        print(f"    Max:  {np.max(rel_diff_50)*100:.2f}%")
        print(f"    Mean: {np.mean(rel_diff_50)*100:.2f}%")
        print(f"    RMS:  {np.sqrt(np.mean(rel_diff_50**2))*100:.2f}%")

    if 'TCA' in Dl_dict and 100 in Dl_dict and Dl_dict[100] is not None:
        mask = Dl_dict[100] > 10
        rel_diff_tca = np.abs(Dl_dict['TCA'][mask] - Dl_dict[100][mask]) / Dl_dict[100][mask]
        print(f"\n  C_l relative difference (TCA vs l_max=100):")
        print(f"    Max:  {np.max(rel_diff_tca)*100:.2f}%")
        print(f"    Mean: {np.mean(rel_diff_tca)*100:.2f}%")

    # ================================================================
    # Performance
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  PERFORMANCE")
    print(f"{'='*60}")
    for label, t in timings.items():
        nvar = n_var_total(label, L_NU_MAX) if isinstance(label, int) else 'N/A'
        print(f"  {label}: {t:.2f}s (N_var={nvar})")

    # ================================================================
    # Plot
    # ================================================================
    print(f"\n--- Generating plot ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        colors = {'TCA': 'gray', 25: 'blue', 50: 'green', 100: 'red'}
        styles = {'TCA': '--', 25: '-', 50: '-', 100: '-'}

        # Panel 1: C_l
        ax = axes[0, 0]
        for label in ['TCA', 25, 50, 100]:
            if label in Dl_dict and Dl_dict[label] is not None:
                lbl = f'$l_{{\\gamma,max}}$={label}' if isinstance(label, int) else 'TCA'
                ax.plot(ell_values, Dl_dict[label],
                        color=colors[label], ls=styles[label], lw=1.2, label=lbl)
        ax.set_xlabel('Multipole $\\ell$')
        ax.set_ylabel(r'$D_\ell$ [$\mu K^2$]')
        ax.set_title('CMB TT: Photon Hierarchy Convergence')
        ax.legend(fontsize=9)
        ax.set_xlim(2, 2500)

        # Panel 2: Ratio to l_max=100
        ax = axes[0, 1]
        if 100 in Dl_dict and Dl_dict[100] is not None:
            Dl_ref = Dl_dict[100]
            for label in ['TCA', 25, 50]:
                if label in Dl_dict and Dl_dict[label] is not None:
                    mask_p = Dl_ref > 10
                    ratio = np.ones_like(Dl_ref)
                    ratio[mask_p] = Dl_dict[label][mask_p] / Dl_ref[mask_p]
                    lbl = f'{label}' if isinstance(label, int) else 'TCA'
                    ax.plot(ell_values, ratio, color=colors[label],
                            ls=styles[label], lw=1, label=lbl)
        ax.axhline(1.0, color='gray', ls=':', alpha=0.5)
        ax.set_xlabel('Multipole $\\ell$')
        ax.set_ylabel(r'$D_\ell / D_\ell^{l_{max}=100}$')
        ax.set_title('Ratio to $l_{\\gamma,max}=100$')
        ax.set_ylim(0.8, 1.2)
        ax.legend(fontsize=9)
        ax.set_xlim(2, 2500)

        # Panel 3: Multipole spectrum
        ax = axes[1, 0]
        has_data = False
        for lgmax in [25, 50, 100]:
            if lgmax in results and results[lgmax] is not None:
                mp = results[lgmax].photon_multipoles()
                ls_list = list(range(lgmax + 1))
                rms = [np.sqrt(np.mean(mp[l]**2)) for l in ls_list]
                ax.semilogy(ls_list, rms, '-o', markersize=2,
                            color=colors[lgmax], label=f'$l_{{max}}$={lgmax}')
                has_data = True
        ax.set_xlabel('Multipole $l$')
        ax.set_ylabel(r'$\sqrt{\langle \Theta_l^2 \rangle}$')
        ax.set_title(r'Photon multipole spectrum at $\tau_{rec}$')
        if has_data:
            ax.legend(fontsize=9)

        # Panel 4: Source functions
        ax = axes[1, 1]
        has_data = False
        for lgmax in [25, 100]:
            if lgmax in results and results[lgmax] is not None:
                sw = results[lgmax].sw_source()
                ax.semilogx(k_arr, sw, color=colors[lgmax], lw=1,
                            label=f'$l_{{max}}$={lgmax}')
                has_data = True
        if 'TCA' in Dl_dict:
            ax.semilogx(k_arr, sw_tca, color='gray', ls='--', lw=1, label='TCA')
            has_data = True
        ax.set_xlabel(r'$k$ [Mpc$^{-1}$]')
        ax.set_ylabel(r'$\Theta_0 + \Psi$')
        ax.set_title('SW Source at Recombination')
        if has_data:
            ax.legend(fontsize=9)

        plt.tight_layout()
        import os
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'hires_convergence.png')
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"  Plot saved: {out_path}")
    except ImportError:
        print("  matplotlib not available")

    # ================================================================
    # Summary
    # ================================================================
    print(f"\n{'='*70}")
    print("SUMMARY: High-resolution photon Boltzmann hierarchy")
    print(f"{'='*70}")
    all_stable = all(lgmax in results and results[lgmax] is not None
                     for lgmax in [25, 50, 100])
    print(f"  Stability: {'ALL PASS' if all_stable else 'SOME FAILED'}")
    for label, t in timings.items():
        print(f"  {label}: {t:.2f}s")
    if 25 in Dl_dict and 100 in Dl_dict and Dl_dict[25] is not None and Dl_dict[100] is not None:
        mask = Dl_dict[100] > 10
        rms_25_100 = np.sqrt(np.mean(((Dl_dict[25][mask] - Dl_dict[100][mask])
                                       / Dl_dict[100][mask])**2)) * 100
        print(f"  C_l RMS difference (l_max=25 vs 100): {rms_25_100:.2f}%")
    if 50 in Dl_dict and 100 in Dl_dict and Dl_dict[50] is not None and Dl_dict[100] is not None:
        mask = Dl_dict[100] > 10
        rms_50_100 = np.sqrt(np.mean(((Dl_dict[50][mask] - Dl_dict[100][mask])
                                       / Dl_dict[100][mask])**2)) * 100
        print(f"  C_l RMS difference (l_max=50 vs 100): {rms_50_100:.2f}%")
    print(f"{'='*70}")

    return all_stable


if __name__ == '__main__':
    _test_hires_convergence()
