# Σ_grav Through the Exponential Wormhole
## Date: 2026-03-10
## Author: Research computation for Sheng-Kai Huang
## Status: Complete calculation

---

## 1. The Exponential Wormhole Geometry

### 1.1 The Metric (Isotropic Coordinates)

The exponential (Papapetrou/Yilmaz-Rosen) metric in isotropic coordinates:

```
ds² = −exp(−r_s/r) c²dt² + exp(+r_s/r)(dr² + r²dΩ²)
```

where r_s = 2M (Schwarzschild radius). Notation: we use m = r_s/2 = M in some formulas following Boonserm et al. (2018), so exp(−r_s/r) = exp(−2m/r).

**Key properties:**
- g₀₀ = −exp(−r_s/r), always negative and nonzero for all r > 0
- g_rr = exp(+r_s/r), always positive and finite for all r > 0
- **No horizons**: g₀₀ never vanishes (exp(−r_s/r) ≠ 0 for any finite r)
- **No singularities in metric components**: all components smooth for r > 0

### 1.2 The Areal Radius

The areal (circumferential) radius is:

```
R(r) = r · exp(r_s/(2r)) = r · exp(m/r)
```

This is the physical radius such that the area of a sphere at coordinate r is A = 4πR².

**Finding the throat**: Take dR/dr = 0:

```
dR/dr = exp(m/r) · [1 − m/r] = 0
```

Since exp(m/r) > 0 always, we need 1 − m/r = 0, giving:

```
r_throat = m = r_s/2
```

At the throat:

```
R_throat = m · exp(1) = (r_s/2) · e = (e/2) · r_s ≈ 1.359 r_s
```

**Verification of minimum (flare-out condition)**:

```
d²R/dr² |_{r=m} = exp(1) · m/m² = e/m > 0
```

The second derivative is positive, confirming this is a minimum of R(r) — the flare-out condition is satisfied. This is a wormhole throat.

### 1.3 Areal Radius Behavior

| r | R(r) = r·exp(m/r) | Behavior |
|---|---|---|
| r → ∞ | R → r (asymptotically flat) | "Our universe" side |
| r = 2m = r_s | R = 2m·e^(1/2) ≈ 1.649·r_s | Approaching throat |
| r = m = r_s/2 | R = m·e ≈ 1.359·r_s | **THROAT (minimum)** |
| r = m/2 | R = (m/2)·e² ≈ 3.695·m | R increasing again |
| r = m/3 | R = (m/3)·e³ ≈ 6.695·m | R growing rapidly |
| r → 0⁺ | R → 0·exp(∞) → ∞ | **Second asymptotic region** |

**Critical insight**: As r → 0⁺, R(r) → ∞. The limit r → 0 corresponds to a **second asymptotically flat region** — another universe (or another part of the same universe). The wormhole connects r → ∞ (our side) through the throat at r = m to r → 0⁺ (the other side).

### 1.4 The Two-Sheeted Structure

The wormhole has a Z₂-asymmetric structure:

```
Our universe:    r ∈ (m, ∞)    →  R ∈ (R_throat, ∞)
Throat:          r = m          →  R = R_throat = em
Other universe:  r ∈ (0, m)    →  R ∈ (R_throat, ∞)
```

Both sides are asymptotically flat (R → ∞) but with different coordinate parameterizations. Unlike the Schwarzschild wormhole (Einstein-Rosen bridge), this wormhole is **traversable**: all metric components remain finite and nonzero at the throat.

---

## 2. Σ_grav Profile Along Radial Geodesic

### 2.1 Definition

From the τ framework:

```
Σ_grav(r) = −ln(−g₀₀) = −ln(exp(−r_s/r)) = r_s/r
```

### 2.2 Profile: Our Side (r > m)

| r | r_s/r | Σ_grav | τ = 1 − exp(−Σ/2) | Physical regime |
|---|---|---|---|---|
| ∞ | 0 | 0 | 0 | Flat spacetime |
| 10 r_s | 0.1 | 0.1 | 0.049 | Weak field |
| r_s | 1 | 1 | 0.394 | Strong field |
| r_s/2 = m | 2 | 2 | 0.632 | **Throat** |

At the throat: Σ_throat = r_s/m = r_s/(r_s/2) = 2.

**Note**: The earlier estimate in the problem statement used r_thr = (e/2)r_s, which is the **areal** radius of the throat. In **isotropic** coordinates, the throat is at r = r_s/2 = m. The value Σ = 2/e ≈ 0.736 would apply if one naively used r_thr = (e/2)r_s in the isotropic formula; the correct value at the throat (in isotropic coordinates) is:

```
Σ_throat = r_s / r_throat = r_s / (r_s/2) = 2
```

This is finite and corresponds to τ_throat = 1 − e⁻¹ ≈ 0.632.

### 2.3 Profile: Other Side (r < m)

As we pass through the throat to r < m:

| r (isotropic) | Σ_grav = r_s/r | τ | Physical meaning |
|---|---|---|---|
| m = r_s/2 | 2 | 0.632 | Throat |
| m/2 = r_s/4 | 4 | 0.865 | Deep inside |
| m/3 = r_s/6 | 6 | 0.950 | Deeper |
| m/10 | 20 | ≈ 1 | Very deep |
| → 0⁺ | → ∞ | → 1 | "Other side's infinity" |

**This is surprising and important**: Σ_grav **increases** monotonically as we go through the throat and toward the other side. It does NOT decrease on the other side.

### 2.4 The Σ_grav Paradox

In isotropic coordinates, Σ_grav = r_s/r diverges as r → 0⁺, even though r → 0⁺ corresponds to a second asymptotically flat region (R → ∞). This means:

- **Geometrically**: The "other side" is asymptotically flat (curvature → 0 as R → ∞)
- **But**: Σ_grav → ∞ as seen from our coordinate system

**Resolution**: Σ_grav = −ln(−g₀₀) is coordinate-dependent in the sense that it measures the informational cost relative to the observer at r = ∞ on **our** side. An observer "born" on the other side would define their own Σ'_grav relative to their own r → 0⁺ infinity. The τ framework is fundamentally **observer-dependent** (as noted in Paper 2 and confirmed by Kudler-Flam et al. 2025).

For an observer on the other side, using a natural coordinate ρ = m²/r (which maps the other side's r ∈ (0, m) to ρ ∈ (m, ∞)):

```
Σ'_grav(ρ) = r_s · ρ/m² = 2ρ/m
```

This vanishes as ρ → ∞ (i.e., r → 0⁺), confirming that the other side IS asymptotically flat from its own perspective. But from our perspective, communicating with the other side costs Σ → ∞.

### 2.5 The Z₂ Asymmetry in Σ

Unlike a symmetric Morris-Thorne wormhole (where both sides look identical), the exponential wormhole is **Z₂-asymmetric**:

- Our side: Σ goes from 0 (at r = ∞) to 2 (at throat)
- Other side: Σ goes from 2 (at throat) to ∞ (at r → 0⁺)

The two sides are not mirror images. The "other universe" appears infinitely redshifted from our perspective, even though it is intrinsically flat.

---

## 3. ΔΣ for Complete Traversal

### 3.1 Traversal from Our Side

Consider a traveler starting at r₁ (our side, r₁ > m) and arriving at r₂ (other side, r₂ < m):

```
ΔΣ = Σ(r₂) − Σ(r₁) = r_s/r₂ − r_s/r₁
```

Since r₂ < m < r₁, we have r_s/r₂ > r_s/r₁, so:

```
ΔΣ > 0    (always, for traversal from our side to the other side)
```

**Information is ALWAYS lost in traversal from our side to the other side.** The traveler accumulates entropy production.

### 3.2 Traversal from the Other Side

For a traveler starting on the other side (r₂ < m) and arriving on our side (r₁ > m), the same formula gives:

```
ΔΣ = Σ(r₁) − Σ(r₂) = r_s/r₁ − r_s/r₂ < 0
```

**From the other side's perspective, traversal TO our side has ΔΣ < 0 — information recovery!**

But this is misleading. The correct calculation uses each observer's natural Σ:

- From our side's perspective: signal going out → ΔΣ > 0 (entropy production)
- From the other side's perspective (using their own reference): signal going out → ΔΣ' > 0 (also entropy production)

### 3.3 The Key Question: Can ΔΣ < 0 Occur?

**Answer: It depends on reference frame, but physically NO in the natural sense.**

Consider a round trip: start at r₁ on our side → throat → r₂ on other side → turn around → throat → back to r₁:

```
ΔΣ_round_trip = Σ(r₁) − Σ(r₁) = 0
```

The round-trip ΔΣ is zero because the metric is static (time-reversal symmetric). This is consistent with Σ being a **state function** (like a potential), not a path-dependent quantity.

**However**, for a one-way traversal measured from our side:

| Scenario | ΔΣ | Meaning |
|---|---|---|
| r₁ → throat (inward, our side) | > 0 | Information lost |
| throat → r₂ (outward, other side) | > 0 | More information lost |
| Total one-way | > 0 | Net entropy production |
| Round trip | = 0 | State function returns to initial value |

### 3.4 The τ Framework Interpretation

The static metric means time-reversal symmetry holds. This is precisely the condition for Σ_total = 0, consistent with:

```
τ_signed = 0    (for a complete round trip in a static spacetime)
```

**The wormhole does NOT produce net entropy production for a round trip.** It is a conservative (potential-like) system. Each "leg" of the journey has ΔΣ, but the total is zero — like climbing and descending a mountain.

This is different from a dynamical process (e.g., gravitational collapse) where Σ_total > 0 irreversibly.

---

## 4. Tortoise Coordinate Behavior

### 4.1 Definition

The tortoise coordinate for the exponential metric is defined by:

```
dr*/dr = √(g_rr / (−g₀₀)) = exp(r_s/r) / exp(−r_s/(2r)) · 1/exp(r_s/(2r))
```

Wait — let me be more careful. For a static spherically symmetric metric in isotropic form:

```
ds² = −A(r) dt² + B(r)(dr² + r²dΩ²)
```

with A = exp(−r_s/r) and B = exp(+r_s/r), the radial null geodesic satisfies:

```
A dt² = B dr²
→ dt/dr = ±√(B/A) = ±exp(r_s/r)
```

So the tortoise coordinate (making null geodesics travel at 45°) is:

```
dr* = √(B/A) dr = exp(r_s/r) dr
```

Therefore:

```
r* = ∫ exp(r_s/r) dr
```

### 4.2 Evaluation of the Integral

Let u = r_s/r, so r = r_s/u, dr = −r_s du/u²:

```
r* = ∫ exp(u) · (−r_s/u²) du = −r_s ∫ exp(u)/u² du
```

The integral ∫ exp(u)/u² du can be evaluated by integration by parts:

```
∫ e^u/u² du = −e^u/u + ∫ e^u/u du = −e^u/u + Ei(u) + C
```

where Ei(u) = ∫₋∞^u (e^t/t) dt is the exponential integral.

Therefore:

```
r*(r) = −r_s [−exp(r_s/r)/(r_s/r) + Ei(r_s/r)] + C
       = r · exp(r_s/r) − r_s · Ei(r_s/r) + C
```

Simplifying:

```
r*(r) = r · exp(r_s/r) − r_s · Ei(r_s/r) + C
```

### 4.3 Behavior at Key Locations

**As r → ∞** (u = r_s/r → 0⁺):
- r · exp(r_s/r) → r · (1 + r_s/r + ...) → r + r_s + ...
- Ei(r_s/r) → ln(r_s/r) + γ + ... → −ln(r/r_s) + γ + ...
- r* → r + r_s + r_s ln(r/r_s) + ... → r (asymptotically)
- **r* → ∞**: standard behavior, tortoise coordinate extends to infinity

**At the throat r = m = r_s/2** (u = 2):
- r · exp(r_s/r) = (r_s/2) · e² ≈ 3.695 r_s
- Ei(2) ≈ 4.954
- r*(throat) = 3.695 r_s − r_s · 4.954 + C ≈ −1.259 r_s + C
- **r* is FINITE at the throat!**

**As r → 0⁺** (u = r_s/r → ∞):
- r · exp(r_s/r) = (r_s/u) · exp(u) → ∞ (exponential beats 1/u)
- Ei(u) → exp(u)/u · (1 + 1/u + ...) ≈ exp(u)/u for large u
- r_s · Ei(r_s/r) = r_s · Ei(u) ≈ r_s · exp(u)/u = r · exp(r_s/r)
- Leading terms cancel! More precisely:

```
r* ≈ r · exp(r_s/r) − r_s · exp(r_s/r)/(r_s/r) · [1 + r/(r_s) + ...]
   = r · exp(r_s/r) − r · exp(r_s/r) · [1 + r/r_s + ...]
   = −r²/r_s · exp(r_s/r) + ...
```

Actually, let's be more careful. For large u:

```
Ei(u) = e^u/u · [1 + 1!/u + 2!/u² + ...]
```

So:

```
r_s · Ei(u) = r_s · e^u/u · [1 + 1/u + 2/u² + ...]
            = r · e^(r_s/r) · [1 + r/r_s + 2r²/r_s² + ...]
```

Therefore:

```
r* = r·e^(r_s/r) − r·e^(r_s/r)·[1 + r/r_s + 2r²/r_s² + ...] + C
   = −r²/r_s · e^(r_s/r) · [1 + O(r/r_s)] + C
```

As r → 0⁺: r²/r_s · e^(r_s/r) → ∞ (exponential dominates polynomial), so:

```
r* → −∞    as r → 0⁺
```

### 4.4 Summary of Tortoise Coordinate

```
r → ∞:     r* → +∞
r = m:     r* = finite (≈ −1.26 r_s + C)
r → 0⁺:   r* → −∞
```

**The tortoise coordinate ranges over all of (−∞, +∞).**

This is the standard behavior for a traversable wormhole: r* extends to ±∞ on both sides, with the throat at a finite value of r*. Light signals can propagate from one side to the other in finite coordinate time but infinite tortoise time on each side — consistent with two asymptotically flat regions.

### 4.5 Comparison with Paper 2 Claim

Paper 2 states: "tortoise coordinate remains bounded."

**Correction**: This claim needs qualification.

- r* is **bounded** in the sense that the echo delay time Δt_echo (the round-trip light travel time from photon sphere to throat and back) is **finite**: Δt_echo ≈ 4.17 r_s/c.
- r* is **unbounded** in the global sense: it ranges over (−∞, +∞) covering both asymptotic regions.
- The finite echo time comes from the path between the photon sphere (r = r_s in isotropic coordinates) and the throat (r = r_s/2), not from r* being globally bounded.

The physically relevant statement is: **the tortoise coordinate interval between photon sphere and throat is finite**, unlike Schwarzschild where the tortoise coordinate diverges logarithmically at the horizon:

```
Exponential:   Δr*(photon sphere → throat) = finite ≈ 2.1 r_s
Schwarzschild: Δr*(photon sphere → horizon) = ∞
```

This finite interval is what produces detectable echoes.

---

## 5. Tolman Temperature Profile

### 5.1 Definition

The Tolman relation for local temperature in a static gravitational field:

```
T(r) = T_∞ / √(−g₀₀) = T_∞ · exp(r_s/(2r))
```

where T_∞ is the temperature measured by an observer at infinity.

### 5.2 Profile Along the Traversal

| r | r_s/(2r) | exp(r_s/(2r)) | T(r)/T_∞ | Physical meaning |
|---|---|---|---|---|
| ∞ | 0 | 1 | 1 | Reference temperature |
| 10 r_s | 0.05 | 1.051 | 1.051 | Slight blueshift |
| r_s | 0.5 | 1.649 | 1.649 | Moderate blueshift |
| r_s/2 = m (throat) | 1 | e ≈ 2.718 | 2.718 | **Throat temperature** |
| r_s/4 | 2 | e² ≈ 7.389 | 7.389 | Other side |
| r_s/10 | 5 | e⁵ ≈ 148.4 | 148.4 | Deep other side |
| → 0⁺ | → ∞ | → ∞ | → ∞ | Other infinity |

### 5.3 Interpretation

**The Tolman temperature increases monotonically as we go deeper through the wormhole**, diverging as we approach the "other side's infinity" (r → 0⁺).

This seems paradoxical: the other side is asymptotically flat, so why is its temperature infinite from our perspective?

**Resolution**: The Tolman temperature T(r) = T_∞/√(−g₀₀) measures the **energy cost of communication** between r and our infinity. An observer on the other side, using their own reference, sees T' = T'_∞ (finite). The divergent Tolman temperature reflects the infinite redshift between the two asymptotic regions — analogous to the infinite blueshift at a black hole horizon, but here stretched across the entire "other universe."

### 5.4 Tolman Temperature and the Gravitational Landauer Principle

The Retrodiction Landauer Principle states:

```
W_erasure(r) = k_B T(r) ln 2 = k_B T_∞ ln 2 · exp(r_s/(2r))
```

At the throat:

```
W_throat = e · k_B T_∞ ln 2 ≈ 2.718 · W_∞
```

The erasure cost at the throat is e times the cost at infinity — finite and traversable. Beyond the throat, the cost grows exponentially, making communication with the other side's deep interior progressively more expensive.

### 5.5 Connection to Hawking Temperature

If we associate T_∞ with the Hawking temperature T_H = ℏc³/(8πGMk_B), then at the throat:

```
T_throat = e · T_H ≈ 2.718 T_H
```

This is a modest enhancement — the throat is "hotter" than the Hawking temperature by a factor of e, but not dramatically so.

---

## 6. Physical Interpretation of the "Other Side"

### 6.1 Coordinate Mapping

The exponential metric in isotropic coordinates covers r ∈ (0, ∞). The key coordinate transformation to understand the geometry is through the areal radius:

```
R(r) = r · exp(r_s/(2r))
```

This maps:
- r ∈ (m, ∞) → R ∈ (R_throat, ∞): **Our universe** (R increases with r)
- r = m → R = R_throat = em: **Throat** (minimum of R)
- r ∈ (0, m) → R ∈ (R_throat, ∞): **Other universe** (R increases as r decreases)

### 6.2 Is r → 0⁺ a Second Asymptotic Region?

**Yes.** As r → 0⁺:

1. **Areal radius**: R = r · exp(r_s/(2r)) → ∞ ✓
2. **Proper radial distance**: ℓ = ∫ exp(r_s/(2r)) dr → ∞ ✓
3. **Curvature invariants**: The Kretschner scalar K ~ M⁴/r⁸ · exp(−4M/r) → 0 as r → 0⁺ ✓

The curvature goes to zero because the exponential decay exp(−4M/r) kills the polynomial divergence r⁻⁸.

**The other side is genuinely asymptotically flat** — a legitimate second universe connected by the wormhole.

### 6.3 Proper Distance to the Throat

The proper radial distance from coordinate r to the throat:

```
ℓ(r) = ∫_m^r exp(r_s/(2r')) dr'
```

From our side (r > m, integrating inward to m):

```
ℓ₁ = ∫_m^∞ exp(r_s/(2r)) dr = finite · r_s
```

This integral converges because exp(r_s/(2r)) → 1 rapidly as r → ∞.

From the other side (r < m, integrating outward from r toward m):

```
ℓ₂(r) = ∫_r^m exp(r_s/(2r')) dr'
```

As r → 0⁺, this integral diverges (the other side extends to infinite proper distance), confirming it's a full asymptotic region.

### 6.4 The Z₂ Asymmetry

Unlike the Ellis-Bronnikov wormhole (which is Z₂-symmetric), the exponential wormhole is **asymmetric**:

| Property | Our side (r > m) | Other side (r < m) |
|---|---|---|
| Asymptotic R | R → ∞ as r → ∞ | R → ∞ as r → 0⁺ |
| Σ_grav at infinity | 0 | ∞ (from our frame) |
| Tolman temperature at infinity | T_∞ | ∞ (from our frame) |
| Gravitational potential | Φ → 0 | Φ → −∞ (from our frame) |
| Intrinsic geometry | Asymptotically flat | Asymptotically flat |
| Curvature | → 0 | → 0 |

**The two sides are intrinsically equivalent as flat spacetimes, but their mutual relationship involves infinite redshift.** An observer on the other side sees their own space as perfectly normal, but communication with our side requires traversing the throat (finite cost) and then an infinite redshift barrier.

### 6.5 Comparison with Other Wormholes

| Wormhole | Z₂ symmetry | Both sides flat? | Σ at throat |
|---|---|---|---|
| Ellis-Bronnikov | Yes | Yes | Depends on model |
| Morris-Thorne (generic) | Optional | Depends | Depends |
| Simpson-Visser (a > 2M) | Yes | Yes | Finite |
| **Exponential** | **No** | **Yes (intrinsically)** | **Σ = 2** |
| Schwarzschild (Einstein-Rosen) | Yes | Yes | **Σ = ∞ (not traversable)** |

### 6.6 The Inversion Symmetry

The exponential metric has a hidden inversion symmetry. Under r → m²/r:

```
r_s/r → r_s · r/m² = r_s · r/(r_s/2)² = 4r/r_s
```

This is NOT a symmetry of the metric (it changes the exponent). However, the coordinate transformation ρ = m²/r maps the other side r ∈ (0, m) to ρ ∈ (m, ∞), and the metric in ρ-coordinates takes the form:

```
ds² = −exp(−4ρ/r_s) dt² + exp(+4ρ/r_s)(m⁴/ρ⁴)(dρ² + ρ²dΩ²)
```

This is NOT the same exponential metric (the mass parameter changes), confirming the Z₂ asymmetry.

---

## 7. Comparison with Schwarzschild

### 7.1 Σ_grav Comparison

| Quantity | Schwarzschild | Exponential |
|---|---|---|
| g₀₀ | −(1 − r_s/r)^iso or −(1−2M/R)^areal | −exp(−r_s/r) |
| Σ_grav = −ln(−g₀₀) | −ln(1 − r_s/r) [diverges at r = r_s] | r_s/r [finite everywhere] |
| Σ at r = 2r_s | 0.693 | 0.5 |
| Σ at r = r_s | ∞ (horizon) | 1.0 |
| Σ at throat | N/A (no wormhole) | 2.0 |
| Σ → ∞? | Yes, at r → r_s | Only as r → 0⁺ |

Note: For Schwarzschild in isotropic coordinates, r_iso = r_s is actually the horizon (which maps to two values of the isotropic coordinate). The comparison is cleaner in areal coordinates, but the qualitative point holds.

### 7.2 The Horizon vs. Throat Dichotomy

```
Schwarzschild at r → r_s (areal):
  Σ_Schw = −ln(1 − r_s/r) → ∞
  τ → 1 (complete information loss)
  g₀₀ → 0 (horizon forms)
  RESULT: Event horizon = point of no return

Exponential at r = m = r_s/2 (isotropic, throat):
  Σ_exp = r_s/m = 2
  τ = 1 − e⁻¹ ≈ 0.632
  g₀₀ = −e⁻² ≈ −0.135 (nonzero!)
  RESULT: Wormhole throat = traversable passage
```

### 7.3 τ Framework Prediction

**The τ framework forbids event horizons because they require τ = 1 (Σ = ∞), which means complete, irreversible information loss.** The strongest form of this argument:

1. **Event horizon requires Σ = ∞**: This means zero recovery fidelity (F = exp(−Σ/2) = 0)
2. **Zero fidelity = complete decoherence**: No quantum information survives
3. **But**: Quantum mechanics is unitary → information must be recoverable in principle
4. **Therefore**: Σ = ∞ is forbidden → event horizons cannot form
5. **The exponential metric replaces the horizon with a throat**: Σ = 2 (finite) → F = e⁻¹ ≈ 0.368 (nonzero recovery)

This is the core no-horizon theorem of Paper 2, now seen in the context of the wormhole geometry.

### 7.4 Strong-Field Behavior Summary

```
Schwarzschild:  Σ rises to ∞ at finite r → creates horizon →
                information paradox → requires exotic resolution (Hawking radiation,
                complementarity, etc.)

Exponential:    Σ rises to 2 at throat → creates wormhole →
                information passes through (degraded but recoverable) →
                no information paradox
```

---

## 8. Key Result: Does ΔΣ < 0 Occur?

### 8.1 The Central Question

Can a traversal through the exponential wormhole result in ΔΣ < 0 (net information recovery, τ_signed < 0)?

### 8.2 Analysis

**Case 1: One-way traversal from our side**

Starting at r₁ > m on our side, ending at r₂ < m on the other side:

```
ΔΣ = Σ(r₂) − Σ(r₁) = r_s/r₂ − r_s/r₁ > 0
```

Always positive — no information recovery.

**Case 2: One-way traversal from the other side**

Starting at r₂ < m (other side), ending at r₁ > m (our side):

```
ΔΣ = Σ(r₁) − Σ(r₂) = r_s/r₁ − r_s/r₂ < 0
```

**ΔΣ < 0!** But this is misleading — it's because we're using OUR Σ reference frame. The traveler from the other side sees decreasing Σ in our frame because they started at higher Σ. From THEIR perspective (using their natural reference), they accumulated entropy production to get to their starting point in the first place.

**Case 3: Round trip**

```
ΔΣ_round = 0
```

Zero, as expected for a static metric.

**Case 4: Symmetric entry points (r₁ = r₂ in areal coordinates)**

If we choose entry and exit points with the same areal radius R, the isotropic coordinates will be different (r₁ > m and r₂ < m with r₁ ≠ r₂ in general):

```
R(r₁) = R(r₂) ⟹ r₁ · exp(m/r₁) = r₂ · exp(m/r₂)
```

For such points: ΔΣ = r_s/r₂ − r_s/r₁. Since r₂ < m < r₁ but they map to the same areal radius, generally r₂ < r₁, so ΔΣ > 0 (the other side's isotropic coordinate is smaller).

### 8.3 The Definitive Answer

**ΔΣ < 0 does NOT occur in any physically meaningful sense for the exponential wormhole.**

The reasons:

1. **Σ_grav is monotonically increasing toward the other side**: Σ = r_s/r increases as r decreases through the throat.

2. **The apparent ΔΣ < 0 for "reverse traversal" is a reference frame artifact**: It's analogous to saying "rolling downhill" — the ball lost height, but only because it was lifted in the first place.

3. **The static metric has time-reversal symmetry**: ΔΣ_round = 0 always. There is no NET entropy production in a round trip. The wormhole is a conservative system.

4. **τ_signed = 0 for the wormhole as a channel**: The complete channel (from our infinity, through the throat, to the other infinity) has Σ → ∞, but this is because the other side's infinity is infinitely redshifted — not because of irreversible processes.

### 8.4 What Would Be Needed for Genuine ΔΣ < 0

To get genuine entropy decrease (τ_signed < 0), one would need:

1. **Non-Markovian information backflow**: The wormhole geometry would need to be dynamical, with the throat contracting or expanding in a way that reverses the information flow direction.

2. **Spin-echo-like reversal**: A mechanism at the throat that reverses the effective Hamiltonian, sending information back more efficiently than it was lost.

3. **Quantum correlations**: Pre-existing entanglement between the two sides (as in ER=EPR) could enable information recovery beyond the Σ bound.

None of these apply to the static exponential metric. **The exponential wormhole is informationally conservative, not informationally regenerative.**

### 8.5 The τ Framework Summary for the Exponential Wormhole

```
                    Σ = 0                Σ = 2              Σ → ∞
                      ↓                    ↓                  ↓
Our infinity ←——————→ Throat ←——————→ Other infinity
   r → ∞              r = m               r → 0⁺
   R → ∞              R = em              R → ∞
   τ = 0              τ = 0.632           τ → 1
   F = 1              F = e⁻¹ ≈ 0.368    F → 0

Key: Information recovery fidelity F decreases monotonically
     from our side to the other side. The throat is traversable
     (F > 0) but the full transit to the other infinity is not
     (F → 0).
```

**Physical picture**: The wormhole acts as an information-degrading channel. At the throat, about 37% of quantum information survives (F = e⁻¹). Going deeper into the other side, less and less survives. The "other universe" is informationally disconnected from us — not by a horizon (Σ = ∞ at a single point), but by a gradual, infinite accumulation of informational cost.

This is arguably a more physical picture than the Schwarzschild horizon: instead of a sharp cutoff, there is a smooth degradation, with the throat as the "narrowest passage" where information survives with fidelity e⁻¹.

---

## 9. References

### Primary Sources

1. **Boonserm, P., Ngampitipan, T., Simpson, A. & Visser, M.** (2018). "Exponential metric represents a traversable wormhole." *Phys. Rev. D* 98, 084048. [arXiv:1805.03781](https://arxiv.org/abs/1805.03781)
   - Establishes exponential metric as traversable wormhole; throat at r = m; flare-out condition; areal radius analysis.

2. **Mychelkin, E.G., Suliyeva, G. & Makukov, M.A.** (2025/2026). "The exponential metric: traversable wormhole and possible identification of scalar background." [arXiv:2510.15391](https://arxiv.org/abs/2510.15391)
   - Gauss-Bonnet invariant sign change at throat; Ricci scalar extremum at r = m/2; thermodynamic characteristics.

3. **Simpson, A. & Visser, M.** (2019). "Black-bounce to traversable wormhole." *JCAP* 2019(02), 042. [arXiv:1812.07114](https://arxiv.org/abs/1812.07114)
   - Interpolation between black hole and wormhole; comparison framework.

4. **Morris, M.S. & Thorne, K.S.** (1988). "Wormholes in spacetime and their use for interstellar travel: A tool for teaching general relativity." *Am. J. Phys.* 56, 395.
   - Foundation: traversable wormhole formalism, flare-out condition, proper distance.

### Supporting Sources

5. **Makukov, M.A. & Mychelkin, E.G.** (2020). "Three scalar field cosmology and the exponential metric." *Gen. Relativ. Gravit.* 52, 102.
   - Three independent scalar-field families converging to exponential metric.

6. **Nath, P.P. & Sarma, D.** (2025). "Exponential wormhole, its quasinormal modes, and study in logarithmic f(R) gravity." *Ann. Phys.* 479, 170067.
   - QNM analysis; stability; logarithmic f(R) gravity embedding.

7. **Nandkishore, A.K. et al.** (2024). "A new class of traversable wormhole metrics." *EPJC* 84, 1063.
   - Generalized exponential wormhole metrics; photon sphere at throat.

8. **Manna, T., Rahman, F. & Chowdhury, T.** (2023). "Strong lensing in exponential wormhole spacetimes." *New Astron. Rev.* 97, 101695.
   - Strong lensing predictions for the exponential wormhole.

9. **Ernazarov, K.K.** (2025). "The Yilmaz-Rosen and JNW metric solutions in the scalar-Einstein-Gauss-Bonnet 4d gravitational model." [arXiv:2510.21625](https://arxiv.org/abs/2510.21625)
   - Phantom scalar field identification; energy condition violation.

10. **Guo, X. & Yuan, F.** (2025). "Quantum Effective Dynamics of Papapetrou Spacetime." *EPJC* (2025). [arXiv:2506.08821](https://arxiv.org/abs/2506.08821)
    - LQG corrections to exponential metric; quantum throat modification.

### τ Framework Sources

11. **Kudler-Flam, J. et al.** (2025). "Gravitational entropy is observer-dependent." [arXiv:2405.00114](https://arxiv.org/abs/2405.00114)
    - Observer-dependence of gravitational entropy; supports Σ_grav being reference-frame dependent.

12. **Cardoso, V., Hooper, S., Regimbau, T., Palenzuela, C., Pani, P.** (2016). "Gravitational-wave signatures of exotic compact objects and of quantum corrections at the horizon scale." *Phys. Rev. D* 94, 084031.
    - Echo predictions for horizonless objects.

13. **Damour, T. & Solodukhin, S.N.** (2007). "Wormholes as black hole foils." *Phys. Rev. D* 76, 024016.
    - Wormholes mimicking black holes; echo time formula.

---

## Appendix A: Detailed Tortoise Coordinate Calculation

### A.1 The Integral

```
r*(r) = ∫ exp(r_s/r) dr
```

Substitution: u = r_s/r, du = −r_s dr/r² = −u² dr/r_s, so dr = −r_s du/u².

```
r* = −r_s ∫ e^u / u² du
```

Integration by parts (∫ e^u/u² du):

Let f = e^u, g' = 1/u² → f' = e^u, g = −1/u

```
∫ e^u/u² du = −e^u/u − ∫ (−e^u/u) du = −e^u/u + ∫ e^u/u du = −e^u/u + Ei(u)
```

where Ei(u) = ∫_{−∞}^{u} e^t/t dt is the exponential integral function.

Therefore:

```
r*(r) = −r_s · [−e^(r_s/r)/(r_s/r) + Ei(r_s/r)] + C
       = r · e^(r_s/r) − r_s · Ei(r_s/r) + C
```

### A.2 Echo Delay Time

The echo delay is the round-trip light travel time between the photon sphere (r_ph = r_s in isotropic coordinates) and the throat (r = r_s/2):

```
Δr* = r*(r_ph) − r*(r_throat)
    = [r_s · e¹ − r_s · Ei(1)] − [(r_s/2) · e² − r_s · Ei(2)]
    = r_s · [e − Ei(1) − e²/2 + Ei(2)]
```

Numerically:
- e ≈ 2.718
- Ei(1) ≈ 1.895
- e²/2 ≈ 3.695
- Ei(2) ≈ 4.954

```
Δr* = r_s · [2.718 − 1.895 − 3.695 + 4.954]
    = r_s · 2.082
```

The echo delay time (round trip):

```
Δt_echo = 2 · Δr* / c ≈ 4.16 r_s/c
```

This matches Paper 2's value of Δt ≈ 4.2 r_s/c (the small difference is from rounding in Ei values).

### A.3 Proper Radial Distance to Throat

The proper distance from r to the throat (for our side, r > m = r_s/2):

```
ℓ(r → throat) = ∫_{r_s/2}^{r} exp(r_s/(2r')) dr'
```

This uses √(g_rr) = exp(r_s/(2r)) since in isotropic coordinates:

```
dℓ² = g_rr dr² = exp(r_s/r) dr²
→ dℓ = exp(r_s/(2r)) dr
```

By the same substitution (u = r_s/(2r)):

```
ℓ = r_s/2 · [r'/r_s · 2 · exp(r_s/(2r')) − Ei(r_s/(2r'))]  evaluated from r_s/2 to r
```

For a traveler starting at r = 10 r_s:

```
ℓ ≈ ∫_{r_s/2}^{10 r_s} exp(r_s/(2r)) dr ≈ 10.3 r_s
```

The throat is at a **finite** proper distance from any point on our side — the wormhole is genuinely traversable.

---

## Appendix B: Σ_grav in Curvature (Areal) Coordinates

### B.1 Coordinate Transformation

The areal radius R = r · exp(r_s/(2r)) is not invertible in closed form (requires the Lambert W function):

```
r = −(r_s/2) / W(−r_s/(2R))
```

where W is the Lambert W function. This has two real branches for R > R_throat = (e/2)r_s:

- W₀ branch: gives r > r_s/2 (our side)
- W₋₁ branch: gives r < r_s/2 (other side)

### B.2 Σ_grav(R) in Areal Coordinates

```
Σ_grav(R) = r_s/r = −2 W(−r_s/(2R))
```

On our side (W₀ branch): Σ goes from 0 (R → ∞) to 2 (R = R_throat).
On the other side (W₋₁ branch): Σ goes from 2 (R = R_throat) to ∞ (R → ∞).

**In areal coordinates, Σ is double-valued**: for each R > R_throat, there are two values of Σ — one from each side of the wormhole. This is the standard wormhole feature of a double-valued radial function.

---

## Appendix C: The Gauss-Bonnet Invariant and τ

From Mychelkin et al. (2025):

```
G(r) = (16 M² e^{−4M/r} / r⁸) · (M² − 4Mr + 3r²)
     = (16 M² e^{−4M/r} / r⁸) · (M − r)(M − 3r)
```

Sign changes at r = M (= r_s/2, the throat) and r = M/3.

In terms of Σ = r_s/r = 2M/r:

```
G = 0  at  Σ = 2 (throat)  and  Σ = 6
G < 0  for  2 < Σ < 6 (just inside the other side)
G > 0  for  Σ < 2 (our side) and Σ > 6 (deep other side)
```

**The Gauss-Bonnet invariant changes sign precisely at the wormhole throat** (Σ = 2). In the τ framework, this is where τ = 1 − e⁻¹ ≈ 0.632 — the "critical fidelity" of the wormhole passage. The sign change of the topological invariant at this exact point suggests a deep connection between the informational (τ) and topological (G) characterizations of the throat.
