# Optimal Warp Bubble Shape — Minimizing Fisher-Information Energy

**Date**: 2026-03-28
**Context**: E_warp = (c⁴/8G) ∫|∇Σ|² d³x from Paper 2 (Σ-field Fisher information)

---

## Setup

The warp bubble is spherically symmetric with Σ(r) = Σ₀ · h(r), where:
- h(r) = 1 for r ≤ R - δ  (interior, flat space for passenger)
- h(r) = 0 for r ≥ R + δ  (exterior, normal space)
- h(r) smooth, monotonically decreasing in the wall region r ∈ [R-δ, R+δ]

The energy functional:

$$E = \frac{c^4 \Sigma_0^2}{8G} \int |\nabla h|^2 \, d^3x = \frac{c^4 \Sigma_0^2}{8G} \cdot 4\pi \int_0^\infty (h')^2 \, r^2 \, dr$$

Since h is constant outside the wall, the integral reduces to:

$$F[h] \equiv \int_{R-\delta}^{R+\delta} (h')^2 \, r^2 \, dr$$

with boundary conditions h(R-δ) = 1, h(R+δ) = 0.

---

## Task 1: Variational Problem — The Optimal Wall Shape

### Euler-Lagrange Equation

Minimize F[h] = ∫(h')² r² dr subject to the boundary conditions.

The integrand is L(r, h, h') = (h')² r². Since L does not depend on h explicitly (only on h'), the Euler-Lagrange equation is:

$$\frac{d}{dr}\left(\frac{\partial L}{\partial h'}\right) = \frac{\partial L}{\partial h}$$

$$\frac{d}{dr}\left(2 h' r^2\right) = 0$$

Therefore:

$$h'(r) \cdot r^2 = C \quad \text{(constant)}$$

$$\boxed{h'(r) = \frac{C}{r^2}}$$

Integrating:

$$h(r) = A - \frac{C}{r}$$

This is the **Coulomb/Newton solution**: h(r) = A + B/r where B = -C.

### Applying Boundary Conditions

h(R - δ) = 1:  A + B/(R-δ) = 1
h(R + δ) = 0:  A + B/(R+δ) = 0

From the second equation: A = -B/(R+δ)

Substituting into the first:

$$-\frac{B}{R+\delta} + \frac{B}{R-\delta} = 1$$

$$B \left(\frac{1}{R-\delta} - \frac{1}{R+\delta}\right) = 1$$

$$B \cdot \frac{2\delta}{R^2 - \delta^2} = 1$$

$$\boxed{B = \frac{R^2 - \delta^2}{2\delta}}$$

$$A = -\frac{B}{R+\delta} = -\frac{R-\delta}{2\delta}$$

### The Optimal Shape Function

$$\boxed{h_{\rm opt}(r) = \frac{R^2 - \delta^2}{2\delta}\left(\frac{1}{r} - \frac{1}{R+\delta}\right) \quad \text{for } r \in [R-\delta, R+\delta]}$$

Equivalently:

$$h_{\rm opt}(r) = \frac{(R+\delta)(R-\delta)}{2\delta} \cdot \frac{R+\delta - r}{r(R+\delta)} = \frac{R-\delta}{2\delta} \cdot \frac{R+\delta - r}{r}$$

**Verification**:
- h(R-δ) = (R-δ)/(2δ) · (2δ)/(R-δ) = 1 ✓
- h(R+δ) = (R-δ)/(2δ) · 0/(R+δ) = 0 ✓

### Physical Interpretation

The optimal shape is **not a linear ramp or sigmoid** — it is a **1/r (Coulomb) profile** within the wall. This makes deep physical sense:

> The 1/r profile is the **harmonic function** on R³ (satisfies ∇²h = 0). Minimizing ∫|∇h|²d³x with fixed boundary values is equivalent to finding the harmonic function — this is the **Dirichlet principle**!

The Euler-Lagrange equation h'r² = const is exactly ∇²h = 0 in spherical symmetry:

$$\nabla^2 h = \frac{1}{r^2}\frac{d}{dr}\left(r^2 \frac{dh}{dr}\right) = 0$$

**This is a classic result from potential theory**: the minimum-energy field configuration with given boundary data is the harmonic (Laplacian-free) solution.

### Optimal Energy

For the optimal shape, h' = -B/r² = -(R²-δ²)/(2δr²):

$$F_{\rm opt} = \int_{R-\delta}^{R+\delta} \frac{B^2}{r^4} \cdot r^2 \, dr = B^2 \int_{R-\delta}^{R+\delta} \frac{dr}{r^2}$$

$$= B^2 \left[\frac{1}{R-\delta} - \frac{1}{R+\delta}\right] = B^2 \cdot \frac{2\delta}{R^2-\delta^2}$$

Using B = (R²-δ²)/(2δ):

$$F_{\rm opt} = \frac{(R^2-\delta^2)^2}{4\delta^2} \cdot \frac{2\delta}{R^2-\delta^2} = \frac{R^2-\delta^2}{2\delta}$$

$$\boxed{F_{\rm opt} = \frac{R^2 - \delta^2}{2\delta}}$$

Full energy:

$$\boxed{E_{\rm opt} = \frac{\pi c^4 \Sigma_0^2}{2G} \cdot \frac{R^2 - \delta^2}{2\delta} = \frac{\pi c^4 \Sigma_0^2 (R^2-\delta^2)}{4G\delta}}$$

---

## Task 2: Compare Bubble Shapes

We compute F[h] = 4π ∫(h')² r² dr for each shape. Define the wall region as r ∈ [R-δ, R+δ].

### (a) Alcubierre Sigmoid: h = 1/(1 + exp((r-R)/σ))

For the sigmoid with width parameter σ ≈ δ/4 (so that the transition effectively occurs within [R-δ, R+δ]):

h'(r) = -1/σ · h(1-h) = -(1/σ) · exp((r-R)/σ) / (1+exp((r-R)/σ))²

Substituting u = (r-R)/σ:

$$F_{\rm sig} = \int_{R-\delta}^{R+\delta} \frac{1}{\sigma^2} \cdot \frac{e^{2u}}{(1+e^u)^4} \cdot r^2 \, dr$$

For R >> δ, we can approximate r² ≈ R² and extend the integral:

$$F_{\rm sig} \approx \frac{R^2}{\sigma^2} \cdot \sigma \int_{-\infty}^{\infty} \frac{e^{2u}}{(1+e^u)^4} du = \frac{R^2}{\sigma} \cdot \frac{1}{6}$$

The integral ∫ sech⁴(u/2)/16 du = 1/6 (standard result: ∫₋∞^∞ e²ᵘ/(1+eᵘ)⁴ du = 1/6).

With σ = δ/4 (chosen so that h(R±δ) ≈ 0.98/0.02, approximately matching the boundary conditions):

$$\boxed{F_{\rm sig} \approx \frac{2R^2}{3\delta}}$$

**Note**: The sigmoid never exactly reaches 0 or 1, so it doesn't strictly satisfy the same BCs. With σ = δ/4, h(R-δ) = 0.982 and h(R+δ) = 0.018, which is close enough for comparison. A wider sigmoid (larger σ) would have less energy but would violate the wall-thickness constraint — it "cheats" by spreading the transition beyond [R-δ, R+δ]. The comparison at σ = δ/4 is the fair one.

### (b) Linear Ramp: h = (R+δ-r)/(2δ) in the wall

h'(r) = -1/(2δ) everywhere in the wall.

$$F_{\rm lin} = \int_{R-\delta}^{R+\delta} \frac{1}{4\delta^2} \cdot r^2 \, dr = \frac{1}{4\delta^2}\left[\frac{r^3}{3}\right]_{R-\delta}^{R+\delta}$$

$$= \frac{1}{4\delta^2} \cdot \frac{(R+\delta)^3 - (R-\delta)^3}{3}$$

Expanding: (R+δ)³ - (R-δ)³ = 6R²δ + 2δ³

$$\boxed{F_{\rm lin} = \frac{6R^2\delta + 2\delta^3}{12\delta^2} = \frac{R^2}{2\delta} + \frac{\delta}{6}}$$

### (c) Optimal (Coulomb/Harmonic): already computed

$$\boxed{F_{\rm opt} = \frac{R^2 - \delta^2}{2\delta} = \frac{R^2}{2\delta} - \frac{\delta}{2}}$$

### (d) Gaussian: h = exp(-(r-R)²/(2δ²))

This doesn't satisfy h → 1 for r << R in the same way. We interpret it as a wall profile that transitions around r = R. More precisely, define h so that it equals 1 for r ≤ R and falls as a Gaussian for r > R (one-sided Gaussian wall of thickness δ). But for fair comparison, let's use a symmetric version:

h(r) = ½[1 + erf((R-r)/(√2 δ_eff))] where δ_eff = δ/2.

Then h' = -1/(√(2π) δ_eff) · exp(-(r-R)²/(2δ_eff²))

For R >> δ, approximate r² ≈ R²:

$$F_{\rm gauss} \approx \frac{R^2}{2\pi \delta_{\rm eff}^2} \cdot \delta_{\rm eff}\sqrt{\pi} = \frac{R^2}{\sqrt{2\pi}\,\delta_{\rm eff}} = \frac{R^2\sqrt{2}}{\sqrt{\pi}\,\delta}$$

$$\boxed{F_{\rm gauss} \approx \frac{R^2}{\sqrt{2\pi}\,(\delta/2)} = \frac{R^2\sqrt{2/\pi}}{\delta} \approx \frac{0.798 \, R^2}{\delta}}$$

### Comparison Table (leading order, R >> δ)

| Shape | F[h] / (R²/δ) | Ratio to Optimal |
|-------|----------------|------------------|
| **(c) Optimal (Coulomb)** | **0.500** | **1.000** |
| (b) Linear ramp | 0.500 + δ²/(6R²) | 1.000 + O(δ²/R²) |
| (a) Alcubierre sigmoid | 0.667 | 1.333 |
| (d) Gaussian (erf) | 0.798 | 1.597 (≈1.60) |

**Key result**: The optimal and linear shapes have the **same leading-order** energy! The linear ramp is near-optimal because in the thin-wall limit (δ << R), the r² weight is approximately constant across the wall, and a linear ramp minimizes ∫(h')² dr with no r² weight. The correction is only O(δ²/R²).

The Alcubierre sigmoid wastes **33% more energy** than optimal.
The Gaussian wastes **~60% more energy** than optimal.

### Exact Ratio (beyond leading order)

For finite δ/R, compare optimal vs linear:

$$\frac{F_{\rm lin}}{F_{\rm opt}} = \frac{R^2/(2\delta) + \delta/6}{R^2/(2\delta) - \delta/2} = \frac{3R^2 + \delta^2}{3R^2 - 3\delta^2} = \frac{3R^2 + \delta^2}{3(R^2 - \delta^2)}$$

For δ/R = 0.1: ratio = 1.014 (linear is 1.4% worse)
For δ/R = 0.3: ratio = 1.132 (13.2% worse)
For δ/R = 0.5: ratio = 1.444 (44.4% worse)

**Bottom line**: For thin walls (δ/R < 0.1), the linear ramp is near-optimal. For thick walls (δ/R > 0.3), the 1/r Coulomb profile wins significantly — the linear ramp wastes >13% energy.

---

## Task 3: Scaling Laws

### General Form

$$E = \frac{c^4 \Sigma_0^2}{8G} \cdot 4\pi \cdot F[h] \equiv \frac{\pi c^4 \Sigma_0^2}{2G} \cdot \mathcal{G}(R,\delta)$$

where G(R,δ) = F[h] is the geometric factor.

### Geometric Factors

| Shape | G(R,δ) |
|-------|--------|
| Optimal | (R² - δ²)/(2δ) |
| Linear | R²/(2δ) + δ/6 |
| Sigmoid | ~2R²/(3δ) |
| Gaussian | ~0.798 R²/δ |

### Scaling with R (bubble size)

$$\boxed{E \propto R^2} \quad \text{(for all shapes, leading order)}$$

Energy scales as the **surface area** of the bubble, not volume. This is because the gradient energy lives on the wall, which is a 2D surface of area ~4πR².

### Scaling with δ (wall thickness)

$$\boxed{E \propto 1/\delta} \quad \text{(for all shapes)}$$

**Thicker walls = less energy**. This is the gradient penalty: spreading the transition over a wider region reduces |∇Σ|². Intuitively, |∇h| ~ 1/(2δ), so |∇h|² ~ 1/(4δ²), times volume ~ R²·2δ, giving R²/δ.

### Scaling with v (velocity, via Σ₀)

From the Σ interpretation: Σ₀ quantifies the metric deformation needed for velocity v. In the weak-field/low-velocity regime:

$$\Sigma_0 \sim \frac{v^2}{c^2}$$

Therefore:

$$\boxed{E \propto \Sigma_0^2 \propto v^4/c^4}$$

**Fourth-power scaling** — not quadratic! This is extremely favorable for subluminal bubbles.

### Optimal δ for Given R

There is **no finite optimal δ from energy minimization alone** — E decreases monotonically as δ increases. The constraint must come from elsewhere:

1. **Bubble usability**: The passenger region has radius R - δ. Requiring this to be positive: δ < R.
2. **Passenger volume constraint**: If we need interior radius r_int = R - δ ≥ L (some minimum size), then δ ≤ R - L.
3. **Fixed interior + exterior**: If we fix the interior radius a = R - δ and exterior radius b = R + δ, then R = (a+b)/2, δ = (b-a)/2, and:

$$E_{\rm opt} = \frac{\pi c^4 \Sigma_0^2}{2G} \cdot \frac{ab}{b-a}$$

This is minimized by making b → ∞ (infinite wall), giving E → πc⁴Σ₀²a/(2G). But that defeats the purpose.

**Practical optimum**: For fixed interior radius a (passenger space), minimize over R:

$$E_{\rm opt}(R) = \frac{\pi c^4 \Sigma_0^2}{4G} \cdot \frac{R^2 - (R-a)^2}{R - a} = \frac{\pi c^4 \Sigma_0^2}{4G} \cdot \frac{2aR - a^2}{R-a}$$

With δ = R - a (so a = R - δ, a fixed), we can express G purely in terms of R:

$$G = \frac{R^2 - \delta^2}{2\delta} = \frac{R^2 - (R-a)^2}{2(R-a)} = \frac{a(2R-a)}{2(R-a)}$$

$$\frac{dG}{dR} = \frac{2a \cdot 2(R-a) - a(2R-a)\cdot 2}{4(R-a)^2} = \frac{2a[2(R-a) - (2R-a)]}{4(R-a)^2} = \frac{2a(2R-2a-2R+a)}{4(R-a)^2} = \frac{-2a^2}{4(R-a)^2}$$

This is **always negative** — G decreases as R increases (for fixed a). So larger R (thicker wall) always wins.

$$\boxed{\text{Optimal strategy: make } \delta \text{ as large as practical (R as large as possible for fixed interior } a\text{)}}$$

**Physical limit**: The wall cannot extend to infinity. External constraints (interaction with ambient spacetime, energy supply geometry, etc.) set a maximum R.

If we impose R ≤ R_max, then optimal is δ = R_max - a:

$$G_{\rm practical} = \frac{a(2R_{\max}-a)}{2(R_{\max}-a)}$$

For R_max >> a: G → a (constant), so E → πc⁴Σ₀²a/(2G). The wall energy becomes negligible and you just pay for the interior.

---

## Task 4: Subluminal Scaling

### The v⁴ advantage

With Σ₀ = v²/c²:

Using the optimal shape in the thin-wall limit (F ≈ R²/(2δ)):

$$E(v) = \frac{\pi c^4 \Sigma_0^2}{2G} \cdot \frac{R^2}{2\delta} = \frac{\pi c^4}{4G\delta} \cdot \frac{v^4}{c^4} \cdot R^2 = \frac{\pi v^4 R^2}{4G\delta}$$

### Reference energy at v = c

$$E(c) = \frac{\pi c^4 R^2}{4G\delta}$$

Numerically for R = 1 m, δ = 0.1 m:

$$E(c) = \frac{\pi (3\times10^8)^4 \times 1}{4 \times 6.674\times10^{-11} \times 0.1} \approx \frac{\pi \times 8.1\times10^{33}}{2.67\times10^{-11}} \approx 9.5 \times 10^{44} \text{ J}$$

Converting to mass: E/c² ≈ 9.5×10⁴⁴ / (9×10¹⁶) ≈ 1.06×10²⁸ kg

That's about **5 Jupiter masses** for a 1-meter light-speed bubble (with δ = 0.1 m). This matches order-of-magnitude estimates in the literature.

### Scaling table

| v | v/c | Σ₀ = v²/c² | E/E(c) = v⁴/c⁴ | E (mass-equiv) |
|---|-----|-------------|-----------------|-----------------|
| c | 1 | 1 | 1 | 1.05×10²⁸ kg ≈ 5.5 M_Jupiter |
| 0.1c | 10⁻¹ | 10⁻² | 10⁻⁴ | 1.05×10²⁴ kg ≈ 0.18 M_Earth |
| 0.01c | 10⁻² | 10⁻⁴ | 10⁻⁸ | 1.05×10²⁰ kg ≈ 0.001 M_Moon |
| 0.001c | 10⁻³ | 10⁻⁶ | 10⁻¹² | 1.05×10¹⁶ kg ≈ small asteroid |
| 10⁴ m/s | 3.3×10⁻⁵ | 1.1×10⁻⁹ | 1.2×10⁻¹⁸ | 1.29×10¹⁰ kg ≈ 13 billion kg |
| 10³ m/s | 3.3×10⁻⁶ | 1.1×10⁻¹¹ | 1.2×10⁻²² | 1.29×10⁶ kg ≈ 1290 tons |
| 10² m/s | 3.3×10⁻⁷ | 1.1×10⁻¹³ | 1.2×10⁻²⁶ | 129 kg |
| 10 m/s | 3.3×10⁻⁸ | 1.1×10⁻¹⁵ | 1.2×10⁻³⁰ | 1.3×10⁻² g |

### Finding the "achievable" velocity threshold

We want E < 10⁶ kg × c² for R = 1 m, δ = 0.1 m:

$$\frac{v^4}{c^4} < \frac{10^6}{1.06 \times 10^{28}} = 9.4 \times 10^{-23}$$

$$\frac{v}{c} < (9.4 \times 10^{-23})^{1/4} = 3.1 \times 10^{-6}$$

$$\boxed{v < 938 \text{ m/s} \approx 1 \text{ km/s} \quad \text{(for } E < 10^6 \text{ kg, R=1m, } \delta\text{=0.1m)}}$$

That's about **Mach 2.7** — supersonic but hardly useful for interstellar travel.

### Optimizing with larger δ

If we increase δ to 1 m (with R = 2 m to keep 1 m interior):

$$E(c) = \frac{\pi c^4 \times 4}{4G \times 1} \approx \frac{\pi \times 8.1\times10^{33} \times 4}{4 \times 6.674\times10^{-11}} = 3.8 \times 10^{44} \text{ J} \implies 4.2 \times 10^{27} \text{ kg}$$

This is only ~4× better. Wall thickness helps as 1/δ but R² hurts.

For truly large δ (R = 100 m, interior a = 1 m, δ = 99 m):

$$G = \frac{R^2 - \delta^2}{2\delta} = \frac{10000 - 9801}{198} = \frac{199}{198} \approx 1.005 \text{ m}$$

$$E(c) = \frac{\pi c^4 \times 1.005}{2G} \approx \frac{\pi \times 8.1\times10^{33}}{2 \times 6.674\times10^{-11}} \approx 1.9 \times 10^{44} \text{ J} \implies 2.1 \times 10^{27} \text{ kg}$$

The geometric factor G → a = 1 m in the large-R limit, so:

$$E_{\min}(c) = \frac{\pi c^4 a}{2G} \approx 2 \times 10^{27} \text{ kg}$$

This is the **irreducible minimum** for a 1-meter interior light-speed bubble.

For this optimized case, E < 10⁶ kg requires:

$$v < c \times \left(\frac{10^6}{2.1\times10^{27}}\right)^{1/4} \approx 4.7 \times 10^{-6} c \approx 1400 \text{ m/s}$$

Only marginally better (~50% improvement over the thin-wall case).

---

## Summary of Key Results

### 1. Optimal shape is the Coulomb (harmonic) profile

$$h_{\rm opt}(r) = \frac{R^2 - \delta^2}{2\delta}\left(\frac{1}{r} - \frac{1}{R+\delta}\right)$$

This follows from the **Dirichlet principle**: minimizing ∫|∇h|²d³x ⟺ solving ∇²h = 0.

### 2. Energy hierarchy

$$E_{\rm optimal} \leq E_{\rm linear} < E_{\rm sigmoid} < E_{\rm Gaussian}$$

The linear ramp is near-optimal (within ~1% for thin walls). Alcubierre's sigmoid wastes 33% energy.

### 3. Universal scaling

$$\boxed{E \propto \frac{R^2}{\delta} \cdot \frac{v^4}{c^4}}$$

- Surface area scaling (R²), not volume
- Inverse wall thickness (1/δ) — thicker is better
- **Fourth power in velocity** — the dominant handle

### 4. Achievability threshold

For R = 1 m bubble, E < 10⁶ kg requires v < ~1 km/s. Not useful for transport.

For E < 1 kg (laboratory demonstration):

$$v < c \times \left(\frac{1}{2.1\times10^{27}}\right)^{1/4} \approx 1.5 \times 10^{-7} c \approx 44 \text{ m/s}$$

A ~44 m/s (160 km/h) warp bubble requiring ~1 kg mass-energy equivalent. Absurd but not physically impossible — it's "merely" an engineering problem of converting rest mass entirely to Σ-field gradient energy.

### 5. Connection to Σ-framework

The energy E = (c⁴/8G)∫|∇Σ|²d³x is the **Fisher information** of the Σ field. The optimal bubble minimizes Fisher information for given boundary conditions. This is the gravitational analog of:

- Electrostatics: minimum ∫|∇φ|²d³x ⟹ Laplace equation
- Elasticity: minimum strain energy ⟹ biharmonic equation
- Information geometry: minimum Fisher information ⟹ maximum ignorance (Σ changes as slowly as possible)

**The optimal warp bubble is the one that "knows" as little as possible about where the wall is — it spreads the information content maximally.**

---

## Appendix: Verification of Coulomb Optimality

**Claim**: Among all h: [R-δ, R+δ] → [0,1] with h(R-δ)=1, h(R+δ)=0, the function h_opt = A + B/r minimizes ∫(h')²r²dr.

**Proof**: Let h = h_opt + ε·η where η(R-δ) = η(R+δ) = 0.

$$F[h_{\rm opt} + \varepsilon\eta] = F[h_{\rm opt}] + 2\varepsilon\int_{R-\delta}^{R+\delta} h_{\rm opt}' \eta' r^2 dr + \varepsilon^2 \int (\eta')^2 r^2 dr$$

The cross term, integrating by parts:

$$\int h_{\rm opt}' \eta' r^2 dr = \left[h_{\rm opt}' \eta r^2\right]_{R-\delta}^{R+\delta} - \int \frac{d}{dr}(r^2 h_{\rm opt}') \eta \, dr = 0 - 0 = 0$$

The boundary term vanishes because η = 0 at endpoints. The integral vanishes because d(r²h'_opt)/dr = 0 (the Euler-Lagrange equation). Therefore:

$$F[h_{\rm opt} + \varepsilon\eta] = F[h_{\rm opt}] + \varepsilon^2 \int (\eta')^2 r^2 dr \geq F[h_{\rm opt}]$$

with equality iff η' = 0, i.e., η = const = 0 (by boundary conditions). □

**The Coulomb profile is the unique global minimizer.**
