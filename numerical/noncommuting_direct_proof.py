#!/usr/bin/env python3
"""
NON-PERTURBATIVE PROOF attempt:

We want to show F^2 > exp(-DeltaD) for [omega, tau] != 0.

F^2 = (1/d) Tr[omega tau^{-1/2} omega tau^{-1/2}]
G   = (1/d) exp(D(omega||tau))

Key identities in the tau-eigenbasis {|e_j>}, tau = sum_j t_j |e_j><e_j|:

F^2 = (1/d) sum_{j,k} |omega_{jk}|^2 / sqrt(t_j t_k)

This is EXACT. Now for G:

G = (1/d) exp(Tr[omega log omega - omega log tau])
  = (1/d) exp(Tr[omega log omega] - sum_j omega_{jj} log(t_j))

Note: Tr[omega log tau] = sum_{j,k} omega_{jk} (log tau)_{kj}
Since tau is diagonal: (log tau)_{kj} = log(t_j) * delta_{jk}
So Tr[omega log tau] = sum_j omega_{jj} log(t_j). This is determined
solely by the diagonal elements of omega.

KEY APPROACH: Use Schur's majorization + convexity.

Let lambda_1, ..., lambda_d be eigenvalues of omega (in decreasing order).
Let w_j = omega_{jj} (diagonal elements in tau-eigenbasis).
By Schur's inequality: (w_1, ..., w_d) is majorized by (lambda_1, ..., lambda_d).

Tr[omega log omega] = sum_k lambda_k log lambda_k = -S(omega)

By the Schur concavity of entropy (or equivalently, convexity of x log x):
sum_k lambda_k log lambda_k >= sum_j w_j log w_j

(Since x*log(x) is convex, and lambdas majorize w's, the sum of convex
function applied to lambdas >= sum applied to w's. Wait, it's the other
way: if lambda majorizes w, then sum f(lambda_i) >= sum f(w_i) for convex f.)

Actually, Schur's theorem says: if lambda majorizes w, then for CONVEX f,
sum f(lambda_i) >= sum f(w_i). Since x*log(x) is convex for x > 0:

sum lambda_k log lambda_k >= sum w_j log w_j

with equality iff lambda = w (up to permutation), which means omega
is diagonal in the tau-eigenbasis.

Therefore:
D(omega||tau) = sum lambda_k log lambda_k - sum w_j log t_j
             >= sum w_j log w_j - sum w_j log t_j
              = sum w_j log(w_j/t_j)
              = D(omega_diag || tau)

With equality iff omega is diagonal in the tau-eigenbasis (i.e., [omega,tau]=0).

So: D(omega||tau) >= D(omega_diag||tau), strict when [omega,tau] != 0.

Now:
F^2 = (1/d)[sum_j w_j^2/t_j + sum_{j!=k} |omega_{jk}|^2/sqrt(t_j*t_k)]
    >= (1/d) sum_j w_j^2/t_j    (since off-diagonal terms are non-negative)

And:
G = (1/d) exp(D(omega||tau))
  >= (1/d) exp(D(omega_diag||tau))    (since D(omega||tau) >= D(omega_diag||tau))

Hmm, the inequality for G goes the wrong way! D(omega||tau) >= D(omega_diag||tau)
means exp(D(omega||tau)) >= exp(D(omega_diag||tau)), so G >= G_diag.
This makes the gap F^2 - G <= F^2 - G_diag.

So we can't just use the diagonal bound on G. We need a tighter upper bound on G.

Let me think about this differently.

Actually, the correct approach is:

F^2 = F^2_diag + F^2_off  where both are non-negative.
G = (1/d) exp(D(omega||tau))

At the saturation condition (commuting case with constant ratio):
F^2_diag = G_diag and F^2_off = 0 and D(omega||tau) = D(omega_diag||tau).

When [omega,tau] != 0:
F^2_off > 0 (at least one off-diagonal element with t_j != t_k)
D(omega||tau) > D(omega_diag||tau) (from Schur)

The first effect INCREASES F^2. The second effect INCREASES G.
These two effects compete. We need to show that the F^2 increase wins.

This is essentially what the perturbative analysis shows at second order.
But for a non-perturbative proof, we need a different approach.

ALTERNATIVE APPROACH: Use the Peierls-Bogoliubov inequality.

For any Hermitian H and density matrix rho:
Tr[rho exp(H)] >= exp(Tr[rho H])

Equivalently (dual): for PSD omega with Tr(omega) = 1:
Tr[omega H] >= log(Tr[exp(log omega + H)])  ... hmm, not quite.

Let me try yet another angle.

APPROACH VIA THE MATRIX POWER MEAN:

F^2 = (1/d) Tr[omega tau^{-1/2} omega tau^{-1/2}]
    = (1/d) Tr[(tau^{-1/4} omega tau^{-1/4})^2]  [by cyclic trace]

This is the SQUARE OF THE HILBERT-SCHMIDT NORM of tau^{-1/4} omega tau^{-1/4},
divided by d.

Let Z = tau^{-1/4} omega tau^{-1/4}. This is PSD with
Tr[Z] = Tr[tau^{-1/2} omega] and
Tr[Z^2] = d * F^2.

Now we can also write Z = tau^{-1/4} omega tau^{-1/4}, so
Z has the same eigenvalues as tau^{-1/2} omega (which are real positive).

The eigenvalues of Z are {z_i}. Then:
- F^2 = (1/d) sum z_i^2  (arithmetic mean of z_i^2)
- Tr[Z] = sum z_i

For the relative entropy:
D(omega||tau) = Tr[omega log omega] - Tr[omega log tau]

Let me use the spectral decomposition omega = U diag(lambda) U^dagger
where U is unitary. Then:
Tr[omega log omega] = sum_k lambda_k log lambda_k

And in the tau-eigenbasis:
Tr[omega log tau] = sum_j omega_{jj} log t_j = sum_j (U diag(lambda) U^dagger)_{jj} log t_j
= sum_j sum_k |U_{jk}|^2 lambda_k log t_j

So D(omega||tau) = sum_k lambda_k log lambda_k - sum_{j,k} |U_{jk}|^2 lambda_k log t_j
= sum_k lambda_k [log lambda_k - sum_j |U_{jk}|^2 log t_j]

Hmm, this is getting complicated. Let me try numerical exploration of the
non-perturbative regime to see if there's a pattern.
"""

import numpy as np
from scipy.linalg import logm, fractional_matrix_power

def compute_F2(omega, tau, d):
    tau_inv_half = fractional_matrix_power(tau, -0.5)
    X = tau_inv_half @ omega @ tau_inv_half
    return np.real(np.trace(omega @ X)) / d

def compute_G(omega, tau, d):
    D_rel = np.real(np.trace(omega @ (logm(omega) - logm(tau))))
    return np.exp(D_rel) / d

def compute_F2_decomposed(omega, tau, d):
    """Decompose F^2 into diagonal and off-diagonal parts in tau-eigenbasis."""
    # Get tau eigenbasis
    t_vals, V = np.linalg.eigh(tau)
    # Transform omega to tau eigenbasis
    omega_rot = V.conj().T @ omega @ V
    F2_diag = 0
    F2_off = 0
    for j in range(d):
        F2_diag += np.abs(omega_rot[j,j])**2 / t_vals[j]
        for k in range(d):
            if k != j:
                F2_off += np.abs(omega_rot[j,k])**2 / np.sqrt(t_vals[j] * t_vals[k])
    return np.real(F2_diag / d), np.real(F2_off / d)

def random_density_matrix(d):
    G = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    rho = G @ G.conj().T
    return rho / np.trace(rho)

# Test the decomposition
np.random.seed(42)
print("=" * 72)
print("EXPLORING THE NON-PERTURBATIVE REGIME")
print("=" * 72)

print("\nDecomposition of F^2 and G for random (omega, tau) pairs:")
print(f"{'d':>3} {'F2':>12} {'F2_diag':>12} {'F2_off':>12} {'G':>12} {'G_diag':>12} {'gap':>12} {'Schur_gap':>12}")
print("-" * 100)

for trial in range(20):
    d = np.random.choice([2, 3, 4])
    omega = random_density_matrix(d)
    tau = random_density_matrix(d)
    tau = 0.8 * tau + 0.2 * np.eye(d) / d
    tau /= np.trace(tau)

    # Compute F^2
    F2 = compute_F2(omega, tau, d)
    F2_diag, F2_off = compute_F2_decomposed(omega, tau, d)

    # Compute G
    G = compute_G(omega, tau, d)

    # Compute G_diag (exp(-DeltaD) for diagonal part)
    t_vals, V = np.linalg.eigh(tau)
    omega_rot = V.conj().T @ omega @ V
    omega_diag_rot = np.diag(np.real(np.diag(omega_rot)))
    tau_diag = np.diag(t_vals)
    G_diag = compute_G(omega_diag_rot, tau_diag, d)

    # Schur gap: D(omega||tau) - D(omega_diag||tau)
    D_full = np.real(np.trace(omega @ (logm(omega) - logm(tau))))
    D_diag = np.real(np.trace(omega_diag_rot @ (logm(omega_diag_rot) - logm(tau_diag))))
    schur_gap = D_full - D_diag

    gap = F2 - G
    print(f"{d:3d} {F2:12.6f} {F2_diag:12.6f} {F2_off:12.6f} {G:12.6f} {G_diag:12.6f} {gap:12.6f} {schur_gap:12.6f}")

print("\n" + "=" * 72)
print("KEY INSIGHT: Can we bound G in terms of G_diag + corrections?")
print("=" * 72)

# The gap is F^2 - G = (F2_diag - G) + F2_off
# = (F2_diag - G_diag) + (G_diag - G) + F2_off
# The first term: F2_diag >= G_diag by AM-GM (commuting case)
# The second term: G_diag - G = (1/d)(exp(D_diag) - exp(D_full))
#   Since D_full >= D_diag (Schur), we have exp(D_full) >= exp(D_diag),
#   so G_diag - G <= 0. This is the "bad" term.
# The third term: F2_off >= 0. This is the "good" term.
#
# So gap = [F2_diag - G_diag] + F2_off - [G - G_diag]
#        = [AM-GM gap, >= 0] + [off-diag F^2, >= 0] - [Schur increase, >= 0]
#
# We need: [AM-GM gap] + F2_off > [G - G_diag]
#
# The perturbative analysis shows that at order eps^2:
# F2_off = O(eps^2) with coefficient 1/sqrt(t_j*t_k)
# G - G_diag = O(eps^2) with coefficient log(t_j/t_k)/(t_j-t_k)
# And 1/sqrt(t_j*t_k) > log(t_j/t_k)/(t_j-t_k) by logarithmic mean inequality.

# But non-perturbatively, G - G_diag can be larger because exp is convex.
# Let's check: is it always true that F2_off >= G - G_diag?

print("\nIs F2_off >= G - G_diag always?")
print(f"{'d':>3} {'F2_off':>12} {'G-G_diag':>12} {'diff':>12} {'OK?':>5}")
print("-" * 50)

violations = 0
for trial in range(10000):
    d = np.random.choice([2, 3, 4, 5])
    omega = random_density_matrix(d)
    tau = random_density_matrix(d)
    tau = 0.8 * tau + 0.2 * np.eye(d) / d
    tau /= np.trace(tau)

    F2 = compute_F2(omega, tau, d)
    F2_diag, F2_off = compute_F2_decomposed(omega, tau, d)
    G = compute_G(omega, tau, d)

    t_vals, V = np.linalg.eigh(tau)
    omega_rot = V.conj().T @ omega @ V
    omega_diag_rot = np.diag(np.real(np.diag(omega_rot)))
    tau_diag = np.diag(t_vals)

    try:
        G_diag = compute_G(omega_diag_rot, tau_diag, d)
    except:
        continue

    diff = F2_off - (G - G_diag)
    if diff < -1e-10:
        violations += 1
        if violations <= 5:
            print(f"{d:3d} {F2_off:12.6f} {G-G_diag:12.6f} {diff:12.6f} {'NO':>5}")

if violations == 0:
    print("  No violations in 10000 tests!")
else:
    print(f"  {violations} violations found!")

print("\n" + "=" * 72)
print("CHECKING: F2_off >= G - G_diag even without AM-GM gap?")
print("=" * 72)

# For the saturation case (all ratios w_j/t_j equal), AM-GM gap = 0.
# Then we need F2_off >= G - G_diag alone.
# This is exactly what the perturbative proof shows at second order.
# But does it hold non-perturbatively?

# Create cases where AM-GM gap is exactly zero (constant ratio)
print("\nCases with constant diagonal ratio w_j/t_j = 1 (AM-GM gap = 0):")
print(f"{'d':>3} {'eps':>8} {'F2_off':>12} {'G-G_diag':>12} {'diff':>12} {'F2-G':>12}")
print("-" * 65)

for d in [2, 3, 4]:
    for eps in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4]:
        # tau with distinct eigenvalues
        t_vals = np.random.dirichlet(np.ones(d))
        # Make sure they're not too small
        t_vals = 0.8 * t_vals + 0.2 / d
        t_vals /= np.sum(t_vals)
        tau = np.diag(t_vals)

        # omega = tau + eps * off-diagonal perturbation
        C = np.random.randn(d, d) + 1j * np.random.randn(d, d)
        C = (C + C.conj().T) / 2  # Hermitian
        C = C - np.diag(np.diag(C))  # off-diagonal only
        C = C / np.linalg.norm(C, 'fro')  # normalize

        omega = tau + eps * C
        # Ensure PSD
        eigvals = np.linalg.eigvalsh(omega)
        if eigvals[0] < 1e-10:
            omega = omega + (abs(eigvals[0]) + 0.01) * np.eye(d)
            omega = omega / np.trace(omega)
            # This changes the diagonal, breaking constant ratio. Skip.
            continue
        omega = omega / np.trace(omega)
        # Check that diagonal ratio is not constant anymore
        # (adding off-diagonal perturbation and renormalizing changes diagonal)
        # Actually, since we only add off-diagonal, Tr is preserved, so omega_jj = t_j.
        # Ratios are w_j/t_j = 1 for all j. AM-GM gap = 0.

        F2 = compute_F2(omega, tau, d)
        F2_diag, F2_off = compute_F2_decomposed(omega, tau, d)
        G = compute_G(omega, tau, d)

        omega_diag = np.diag(np.real(np.diag(omega)))
        G_diag = compute_G(omega_diag, tau, d)

        diff = F2_off - (G - G_diag)
        print(f"{d:3d} {eps:8.4f} {F2_off:12.6f} {G-G_diag:12.6f} {diff:12.6f} {F2-G:12.6f}")

print("\n" + "=" * 72)
print("NON-PERTURBATIVE PROOF STRATEGY")
print("=" * 72)

print("""
The full non-perturbative proof requires showing:

  F^2_off >= G - G_diag  whenever omega_diag = tau (constant ratio = 1)

where:
  F^2_off = (1/d) sum_{j!=k} |omega_{jk}|^2 / sqrt(t_j*t_k)
  G - G_diag = (1/d)[exp(D(omega||tau)) - exp(D(omega_diag||tau))]
             = (1/d)[exp(D(omega||tau)) - 1]   (since D(tau||tau) = 0)

Since D(omega||tau) >= 0 (Gibbs' inequality), G >= 1/d.

The key is: exp(D(omega||tau)) - 1 vs sum_{j!=k} |omega_{jk}|^2/sqrt(t_j*t_k)

For omega = tau + C (C off-diagonal, Tr C = 0):

D(omega||tau) = Tr[(tau+C)log(tau+C)] - Tr[(tau+C)log tau]
              = Tr[(tau+C)log(tau+C)] - Tr[tau log tau] - 0
              = -S(tau+C) + S(tau) + Tr[(tau+C)log(tau+C)] + S(tau+C) - Tr[tau log tau] - S(tau)

Hmm, let me use a cleaner notation. Actually for the constant-ratio case
(omega_diag = tau), D(omega_diag||tau) = D(tau||tau) = 0, so:

  G = (1/d) exp(D(omega||tau))
  G_diag = 1/d

The question becomes: is F^2_off >= (1/d)[exp(D(omega||tau)) - 1]?

i.e., sum_{j!=k} |omega_{jk}|^2/sqrt(t_j*t_k) >= exp(D(omega||tau)) - 1?

Note that exp(x) - 1 <= x*exp(x) for x >= 0 is NOT helpful directly.
But exp(x) - 1 >= x for x >= 0 (by convexity of exp).

So we need: sum |omega_{jk}|^2/sqrt(t_j*t_k) >= D(omega||tau).

And by the Schur analysis:
D(omega||tau) = -S(omega) - Tr[omega log tau]
              = -S(omega) - sum_j omega_{jj} log t_j
              = -S(omega) - sum_j t_j log t_j   (since omega_jj = t_j)
              = -S(omega) + S(tau)  ... wait, S(tau) = -sum t_j log t_j.
              = S(tau) - S(omega)

YES! For omega_diag = tau:
D(omega||tau) = Tr[omega log omega - omega log tau]
= Tr[omega log omega] - Tr[omega log tau]
= -S(omega) - (-S(tau))   [since Tr[omega log tau] = sum t_j log t_j = -S(tau)]
Wait no: -S(omega) = Tr[omega log omega], and
Tr[omega log tau] = sum_j omega_{jj} log t_j = sum_j t_j log t_j = -S(tau)
So D(omega||tau) = -S(omega) - (-S(tau)) = S(tau) - S(omega).

Since omega = tau + C with C off-diagonal Hermitian, the eigenvalues of omega
majorize the diagonal elements (which equal the eigenvalues of tau).
By Schur concavity of entropy: S(omega) <= S(omega_diag) = S(tau).
So D(omega||tau) = S(tau) - S(omega) >= 0. Good, this confirms Gibbs' inequality.

Now we need: sum_{j!=k} |omega_{jk}|^2/sqrt(t_j*t_k) >= exp(S(tau) - S(omega)) - 1

Hmm, this is still hard to prove non-perturbatively in general.
But at second order in the off-diagonal elements:
S(tau) - S(omega) ~ sum_{j!=k} |omega_{jk}|^2 * log(t_j/t_k)/(t_j - t_k) + higher order

And exp(S(tau) - S(omega)) - 1 ~ S(tau) - S(omega) + O((S(tau)-S(omega))^2)
~ sum |omega_{jk}|^2 * log(t_j/t_k)/(t_j-t_k) + higher order

So the question reduces to showing that the "higher order" terms don't win.

Actually, I realize there might be a cleaner approach using a known matrix inequality.
""")

# Let me test a specific conjecture: is D(omega||tau) a CONCAVE function
# of the off-diagonal elements? If so, then exp(D) - 1 is bounded by the
# tangent at zero, giving exactly the second-order result.

print("Testing concavity of D(omega||tau) in off-diagonal perturbation:")
print()

for d in [2, 3]:
    for trial in range(3):
        t_vals = np.random.dirichlet(np.ones(d))
        t_vals = 0.8 * t_vals + 0.2 / d
        t_vals /= np.sum(t_vals)
        tau = np.diag(t_vals)

        # Random off-diagonal direction
        C = np.random.randn(d, d) + 1j * np.random.randn(d, d)
        C = (C + C.conj().T) / 2
        C = C - np.diag(np.diag(C))
        C = C / np.linalg.norm(C, 'fro')

        print(f"  d={d}, tau eigenvalues = {t_vals.round(4)}")
        eps_values = np.linspace(0, 0.3, 30)
        D_values = []
        F2off_values = []
        for eps in eps_values:
            omega = tau + eps * C
            eigvals = np.linalg.eigvalsh(omega)
            if eigvals[0] < 1e-10:
                break
            omega = omega / np.trace(omega)  # renormalize

            D_val = np.real(np.trace(omega @ (logm(omega) - logm(tau))))
            D_values.append(D_val)

            _, F2_off = compute_F2_decomposed(omega, tau, d)
            F2off_values.append(F2_off)

        # Check if D(eps) is concave in eps
        D_values = np.array(D_values)
        eps_valid = eps_values[:len(D_values)]

        if len(D_values) > 5:
            # Second derivative estimate
            d2D = np.diff(D_values, 2) / np.diff(eps_valid, 2)[:-1]**2
            # Wait, this isn't right. Let me use finite differences properly.
            h = eps_valid[1] - eps_valid[0]
            d2D = np.diff(D_values, 2) / h**2
            is_concave = np.all(d2D <= 1e-6)
            print(f"    D(eps) concave? {is_concave}, max(d2D) = {max(d2D):.6f}")
        print()

print("D(omega||tau) is NOT concave in general (it can be convex).")
print("So the simple concavity argument fails.")
print()

# BUT: The exponential exp(D) also grows, and the competition
# between F^2_off and exp(D)-1 is the key. Let's check directly.

print("=" * 72)
print("DIRECT CHECK: F^2_off vs exp(D(omega||tau)) - 1")
print("=" * 72)

# For constant diagonal ratio (omega_jj = t_j):
np.random.seed(123)
violations = 0
n_tested = 0

for trial in range(50000):
    d = np.random.choice([2, 3, 4, 5])
    t_vals = np.random.dirichlet(np.ones(d))
    t_vals = 0.8 * t_vals + 0.2 / d
    t_vals /= np.sum(t_vals)
    tau = np.diag(t_vals)

    # Random off-diagonal perturbation (varying magnitude)
    eps = np.random.uniform(0.001, 0.5)
    C = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    C = (C + C.conj().T) / 2
    C = C - np.diag(np.diag(C))
    C = C / np.linalg.norm(C, 'fro') * eps

    omega = tau + C
    eigvals = np.linalg.eigvalsh(omega)
    if eigvals[0] < 1e-10:
        continue

    # Don't renormalize - keep omega_jj = t_j exactly
    # Actually Tr(omega) = Tr(tau) + Tr(C) = 1 + 0 = 1 since C is traceless.
    # And omega_jj = t_j since C is off-diagonal.

    try:
        F2 = compute_F2(omega, tau, d)
        G = compute_G(omega, tau, d)
        gap = np.real(F2 - G)
        n_tested += 1

        if gap < -1e-10:
            violations += 1
    except:
        continue

print(f"Tested {n_tested} cases with constant diagonal ratio, {violations} violations")
if violations == 0:
    print("CONFIRMED: F^2 > exp(-DeltaD) even for large off-diagonal perturbations.")
