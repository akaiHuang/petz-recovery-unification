# Route 5: Unified Sigma Formula -- Running G, MOND, and Galaxy Rotation Curves

**Date**: 2026-03-07
**Agent task**: Investigate whether a single QRE formula Sigma = D(rho_spacetime || rho_matter) can produce different gravitational behaviors at different scales, and how this compares to MOND

## Executive Summary

A single QRE formula CAN in principle produce different gravitational behaviors at different scales, via the DPI-constrained gravitational channel N_grav(k). The synthesis of Casini-Huerta (RG = DPI), Dorau-Much (QRE = Einstein), and Kumar (running G = flat rotation curves) forms a coherent chain. However, critical gaps remain: deriving eta=1 from information theory, explaining CMB acoustic peaks, and making cluster-scale predictions. The framework is fundamentally different from MOND (modifies information processing, not acceleration law) and avoids MOND's relativistic completion problem.

## Key Findings

### 1. Running G from QRE/DPI: The Casini-Huerta Bridge

Casini, Teste & Torroba (2017, arXiv:1611.00016) proved that RG coarse-graining IS a quantum channel, and QRE monotonicity under RG flow IS the DPI. The logical chain:

1. RG flow = quantum channel (coarse-graining) -- **established**
2. QRE monotonicity under this channel = DPI -- **established** (Casini-Huerta 2017)
3. QRE --> Einstein equations at each scale -- **established** (Dorau-Much 2025, PRL)
4. Therefore: scale-dependent QRE --> scale-dependent Einstein equations --> running G

**Step 4 is genuinely novel. No one has done this synthesis.**

### 2. Kumar (2025) -- Marginal IR Running of G

From arXiv:2509.05246 (Physics Letters B 871, 2025):

```
Modified potential:  Phi(r) = -G_N M / r + (2 G_N M k_*) / pi * ln(r/r_0)
Force law:           F(r) = -G_N M / r^2 - (2 G_N M k_*) / (pi r)
Rotation curve:      v^2(r) = G_N M_bar / r + 2 G_N M_bar k_* / pi
```

At r >> r_c: v^2 --> constant (flat rotation curves).
Crossover scale: r_c = 1/k_* ~ 36-38 kpc (universal across 3 galaxies tested).

**Important caveat**: eta=1 (marginal running) comes from dimensional/scaling argument, NOT a first-principles derivation. More principled than MOND's a_0 but not derived from a complete theory.

### 3. Key Differences from MOND

| Feature | MOND | Running G | QRE Framework |
|---------|------|-----------|---------------|
| What is modified | Acceleration law | Gravitational coupling G(r) | Gravitational quantum channel N_grav(k) |
| Origin of scale | a_0 phenomenological | k_* from RG | Emerges from de Sitter IR cutoff |
| External Field Effect | Yes (unique prediction) | No | Potentially (channel depends on environment) |
| Relativistic completion | TeVeS (problematic, ruled out by GW170817) | Natural within GR | Built on GR from start |
| Gravitational waves | v_GW != c in many extensions | v_GW = c | v_GW = c |

### 4. Where MOND Fails and Whether We Avoid Those Failures

| Failure Mode | MOND Status | Running G / Our Status |
|-------------|-------------|----------------------|
| Galaxy clusters | Factor ~2 missing mass remains | Untested -- open calculation needed |
| CMB acoustic peaks | Cannot explain 3rd peak without DM | Running G contributes only 0.05% at z~1100 -- **severe problem** |
| Bullet Cluster | Cannot explain mass/baryon offset without DM | Moffat's MOG claims success but needs extra vector field |
| Gravitational lensing | Requires TeVeS (ruled out) | Automatic in GR + running G |

### 5. The Three-Scale Picture

| Scale | Dominant Sigma Term | Physics | Source |
|-------|-------------------|---------|--------|
| Strong field (r ~ r_s) | Sigma ~ r_s/r (area-law, classical G) | Schwarzschild geometry | Paper 2 |
| Galactic (r ~ 10 kpc) | Sigma ~ r_s/r + (2k_*r_s/pi)*ln(r/r_0) | Running G from IR RG flow | Paper 3 |
| Cosmological (r ~ L_dS) | Sigma includes volume-law de Sitter contribution | Entropy displacement (Verlinde) | Paper 3 |

All from: Sigma = D(rho_spacetime || rho_matter) - D(N_grav(k)(rho_spacetime) || N_grav(k)(rho_matter))

### 6. Verlinde (2016) Connection

Verlinde's emergent gravity (arXiv:1611.02269) provides the conceptual foundation for the cosmological-scale contribution:
- De Sitter volume-law entanglement displaced by matter produces additional gravitational force
- Predicts Omega_D ~ 0.26 from Omega_B ~ 0.05
- Ghari & Haghi (2026, arXiv:2601.01715): Verlinde beats MOND at 5.2 sigma for 23 dwarf spheroidals
- BUT: fails at cluster cores (>5 sigma rejection), internal consistency questioned

## Key Equations

```
Running G:     G(k) ~ G_N * (k_*/k)^eta,  eta = 1 (marginal)

Modified force: F(r) = -G_N M/r^2 - (2G_N M k_*)/(pi r)

Flat v^2:      v^2(r >> r_c) = 2 G_N M_bar k_* / pi = const

Sigma unified: Sigma = D(rho_spacetime || rho_matter) with scale-dependent N_grav(k)

Verlinde:      g_D = sqrt(a_0 * g_B / 6),  a_0 = c*H_0
```

## Key References

### Primary
- Casini, Teste & Torroba (2017): arXiv:1611.00016 -- Relative entropy and RG flow
- Kumar (2025): arXiv:2509.05246 -- Marginal IR running of Gravity (PLB 871)
- Gubitosi et al. (2024): arXiv:2403.00531 -- RGGR phenomenology with SPARC
- Dorau & Much (2025 PRL): arXiv:2510.24491 -- QRE to semiclassical Einstein equations
- Verlinde (2016): arXiv:1611.02269 -- Emergent gravity and the dark universe
- Bianconi (2025 PRD): arXiv:2408.14391 -- Gravity from entropy
- Ghari & Haghi (2026): arXiv:2601.01715 -- Verlinde vs MOND in dwarf spheroidals
- Reuter & Weyer (2004): arXiv:hep-th/0410117 -- Running Newton constant and rotation curves

### Cluster and CMB Tests
- Brouwer et al. (2017): arXiv:1612.03034 -- First test of Verlinde with weak lensing
- Testing emergent gravity on clusters (2019): arXiv:1901.05505 -- Verlinde fails at cluster cores
- Logarithmic corrections and LSS (2021): arXiv:2102.09602

### Observational Discriminators
- Hooks & Bends in the RAR (MNRAS 2024) -- Non-monotonic RAR tracks challenge MOND
- CMB constraints on G (2025): A&A -- G consistent with lab value to 1.8%
- Dai & Stojkovic (2017): arXiv:1710.00946 -- Verlinde internal inconsistencies

## Assessment

- **Likelihood of unified Sigma producing all scales**: MEDIUM-HIGH (chain is logically coherent but key derivations missing)
- **What works**: Casini-Huerta + Dorau-Much + Kumar chain is genuinely novel; running G avoids MOND's relativistic problems; flat rotation curves emerge naturally
- **What doesn't**: eta=1 not derived from QRE; CMB is severe unsolved problem; cluster predictions absent
- **Genuine novelty**: The synthesis connecting QRE scale-dependence to running G to rotation curves has not been done by anyone

## Honest Gaps

1. **eta = 1 from QRE**: Central missing derivation -- why the gravitational channel's information loss runs marginally
2. **CMB acoustic peaks**: Running G contributes only ~0.05% at recombination -- cannot explain 3rd peak without additional physics
3. **Cluster-scale predictions**: No quantitative predictions exist
4. **N_grav(k) definition**: Gravitational channel at arbitrary scale not rigorously defined in AQFT
5. **Verlinde consistency**: Internal self-consistency still debated

## Next Steps

1. **For Paper 3**: Present the three-scale picture with honest gap acknowledgment
2. **Derive eta=1**: Try to show marginal running follows from DPI saturation or information-theoretic channel capacity bounds
3. **CMB strategy**: Investigate whether Verlinde's volume-law entropy at z~1100 can substitute for CDM in driving acoustic peaks
4. **Cluster calculations**: Compute running G predictions for cluster mass profiles
5. **Test against SPARC**: Compare QRE-derived rotation curves vs MOND vs CDM on the full 175-galaxy SPARC dataset (Gubitosi et al. did 5 galaxies)
6. **External Field Effect**: Calculate whether N_grav(k) has environment dependence that mimics MOND's EFE
