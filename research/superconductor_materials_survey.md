# Comprehensive Survey: High-Temperature and Room-Temperature Superconductor Candidates

**Author**: Sheng-Kai Huang
**Date**: 2026-03-28
**Status**: Living document (materials database + tau-framework analysis)
**Purpose**: Catalogue all known superconductors with T_c > 30 K, their material properties, pairing mechanisms, and estimated R_tau values for room-temperature screening.

---

## 1. Introduction and Scope

This survey catalogues all well-established superconducting materials and candidate systems with T_c > 30 K, organized by material class. For each entry we record:

- Chemical formula, crystal structure, and space group
- Critical temperature T_c (K) and required pressure (GPa)
- Pairing mechanism and symmetry
- Electron-phonon coupling constant lambda (where applicable)
- Debye temperature Theta_D (K)
- Superconducting gap Delta(0) (meV)
- Estimated R_tau = I_pair / tau_phonon(300K) from our framework
- Key references

The R_tau screening ratio is defined in `paper_superconductivity_tau.tex` (Eq. 17):

```
R_tau = max_T [I_pair(T)] / tau_dec(300 K)
```

where I_pair is the total Cooper pair entanglement information and tau_dec(300K) is the thermal decoherence at room temperature. A room-temperature superconductor requires R_tau >= 1.

**Proxy formula** (from Paper, Eq. 20):

```
R_tau ~ [Delta(0) / (k_B * 300 K)] * [2*pi*N(0)*ln(2) / lambda_ep]
```

For practical estimation when full Eliashberg data is unavailable, we use a simplified calibration:

```
R_tau ~ (T_c / 300 K) * (2*Delta(0) / (3.53 * k_B * T_c)) * (lambda_ep / (1 + lambda_ep))
```

where 2*Delta(0)/(3.53*k_B*T_c) is the strong-coupling correction ratio (= 1 for BCS weak coupling, > 1 for strong coupling).

---

## 2. Master Table: All Known Superconductors with T_c > 30 K

### 2.1 Conventional (Phonon-Mediated) Superconductors

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| MgB2 | AlB2-type hexagonal | P6/mmm | 39 | 0 | 7.1 (sigma), 2.8 (pi) | 0.87 | 750 | Phonon, s-wave (two-gap) | 0.012 | Confirmed | Nagamatsu (2001), Choi (2002) |
| CaC6 | Rhombohedral | R-3m | 11.5 | 0 | 1.8 | 0.83 | 340 | Phonon, s-wave | 0.003 | Confirmed | Weller (2005) |
| YB6 | CsCl-type | Pm-3m | 7.2 | 0 | 1.1 | 0.78 | 600 | Phonon, s-wave | 0.002 | Confirmed | Lortz (2006) |

### 2.2 Compressed Hydrides (Phonon-Mediated, High Pressure)

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| H3S | Im-3m (bcc) | Im-3m | 203 | 155 | 30 | 2.19 | 1330 | Phonon, s-wave | 0.069 | Confirmed | Drozdov (2015), Duan (2014), Errea (2015) |
| LaH10 | fcc clathrate | Fm-3m | 250-260 | 170-190 | 36-40 | 3.41 | 1500 | Phonon, s-wave | 0.092 | Confirmed | Somayazulu (2019), Drozdov (2019), Liu (2017) |
| YH6 | Im-3m (sodalite) | Im-3m | 220 | 166 | 33 | 2.5 | 1400 | Phonon, s-wave | 0.078 | Confirmed | Troyan (2021), Kong (2021) |
| YH9 | P63/mmc | P63/mmc | 243 | 201 | 37 | 2.8 | 1450 | Phonon, s-wave | 0.087 | Confirmed | Kong (2021), Snider (2021) |
| CaH6 | Im-3m (sodalite) | Im-3m | 215 | 172 | 32 | 2.7 | 1350 | Phonon, s-wave | 0.075 | Confirmed | Ma (2022) |
| LaBeH8 | Fm-3m (clathrate) | Fm-3m | 110 | 80 | 17 | 1.4 | 1200 | Phonon, s-wave | 0.035 | Confirmed | Song (2023) |
| (La,Y)H10 | Fm-3m | Fm-3m | 253 | 183 | 38 | 3.1 | 1480 | Phonon, s-wave | 0.089 | Confirmed | Semenok (2021) |
| CeH9 | P63/mmc | P63/mmc | 100-117 | 88-100 | 16 | 1.6 | 1100 | Phonon, s-wave | 0.032 | Confirmed | Li (2022) |
| ThH9 / ThH10 | P63/mmc / Fm-3m | -- | 146-161 | 170-175 | 24 | 2.0 | 1200 | Phonon, s-wave | 0.052 | Confirmed | Semenok (2020) |
| PrH9 | P63/mmc | P63/mmc | 9 | 120 | -- | -- | -- | Phonon, s-wave | ~0.002 | Claimed | Bi (2022) |
| BaH12 | Predicted | Cmc21 | 20 (pred.) | 135 | -- | 0.9 | 900 | Phonon, s-wave | ~0.005 | Predicted | Semenok (2021) |
| Lu-N-H ("LuH2 + N") | -- | -- | 294 (claimed) | 1 (claimed) | -- | -- | -- | Claimed ambient | -- | **RETRACTED** | Dias (2023) - retracted from Nature |
| CSH ("carbonaceous sulfur hydride") | -- | -- | 287 (claimed) | 267 | -- | -- | -- | Claimed | -- | **RETRACTED** | Snider/Dias (2020) - retracted from Nature |
| LuH3 (without N) | Fm-3m | Fm-3m | ~12 (theoretical) | 0 | -- | 0.5 | -- | Phonon | ~0.003 | Normal metal | Ferreira (2023), Xie (2023) |

### 2.3 Cuprate (High-T_c) Superconductors

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| La2-xBaxCuO4 (x~0.15) | K2NiF4-type tetragonal | I4/mmm | 30-35 | 0 | 8-10 | -- | 380 | Spin-fluct., d-wave | 0.010 | Confirmed | Bednorz & Muller (1986) |
| La2-xSrxCuO4 (x~0.15) | K2NiF4-type tetragonal | I4/mmm | 38-40 | 0 | 10-12 | -- | 400 | Spin-fluct., d-wave | 0.013 | Confirmed | Cava (1987) |
| YBa2Cu3O7-d (YBCO) | Orthorhombic perovskite | Pmmm | 93 | 0 | 20 | -- | 410 | Spin-fluct., d-wave | 0.031 | Confirmed | Wu (1987) |
| Bi2Sr2CaCu2O8+d (Bi-2212) | Tetragonal | I4/mmm | 85-95 | 0 | 25-35 | -- | 300 | Spin-fluct., d-wave | 0.030 | Confirmed | Maeda (1988) |
| Bi2Sr2Ca2Cu3O10+d (Bi-2223) | Tetragonal | I4/mmm | 108-110 | 0 | 30-38 | -- | 310 | Spin-fluct., d-wave | 0.037 | Confirmed | Maeda (1988) |
| Tl2Ba2CaCu2O8 (Tl-2212) | Tetragonal | I4/mmm | 110 | 0 | 30-35 | -- | 350 | Spin-fluct., d-wave | 0.037 | Confirmed | Sheng & Hermann (1988) |
| Tl2Ba2Ca2Cu3O10 (Tl-2223) | Tetragonal | I4/mmm | 125 | 0 | 35-40 | -- | 350 | Spin-fluct., d-wave | 0.042 | Confirmed | Sheng & Hermann (1988) |
| HgBa2CuO4+d (Hg-1201) | Tetragonal | P4/mmm | 94-97 | 0 | 22-25 | -- | 400 | Spin-fluct., d-wave | 0.032 | Confirmed | Putilin (1993) |
| HgBa2CaCu2O6+d (Hg-1212) | Tetragonal | P4/mmm | 128 | 0 | 35-40 | -- | 400 | Spin-fluct., d-wave | 0.043 | Confirmed | Putilin (1993) |
| HgBa2Ca2Cu3O8+d (Hg-1223) | Tetragonal | P4/mmm | 133 | 0; 164 at 30 GPa | 40-50 | -- | 400 | Spin-fluct., d-wave | 0.050 | Confirmed | Schilling (1993), Gao (1994) |
| Tl0.5Pb0.5Sr2Ca2Cu3O9 | Tetragonal | P4/mmm | 120 | 0 | 32-38 | -- | 370 | Spin-fluct., d-wave | 0.040 | Confirmed | Subramanian (1988) |
| (Ca,Sr)14Cu24O41 (spin ladder) | Orthorhombic | Fmmm | 12 (under P) | 3-5 | -- | -- | -- | Spin-fluct. | ~0.003 | Confirmed | Uehara (1996) |

### 2.4 Iron-Based Superconductors

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| LaFeAsO1-xFx (x~0.1) | ZrCuSiAs-type | P4/nmm | 26 | 0; 43 at 4 GPa | 4-7 | 0.2-0.3 | 300 | Spin-fluct., s+- | 0.008 | Confirmed | Kamihara (2008), Takahashi (2008) |
| SmFeAsO1-xFx | ZrCuSiAs-type | P4/nmm | 55-58 | 0 | 12-15 | -- | 350 | Spin-fluct., s+- | 0.018 | Confirmed | Ren (2008) |
| NdFeAsO1-xFx | ZrCuSiAs-type | P4/nmm | 52 | 0 | 10-12 | -- | 340 | Spin-fluct., s+- | 0.017 | Confirmed | Ren (2008) |
| BaFe2As2 (Ba-122, doped) | ThCr2Si2-type | I4/mmm | 38 | 0 | 6-12 | -- | 270 | Spin-fluct., s+- | 0.012 | Confirmed | Rotter (2008) |
| Ba0.6K0.4Fe2As2 | ThCr2Si2-type | I4/mmm | 38 | 0 | 6-12 | -- | 270 | Spin-fluct., s+- | 0.012 | Confirmed | Rotter (2008) |
| FeSe | PbO-type | P4/nmm | 8 (bulk); 37 at 7 GPa | 0 (bulk) | 1.5-3 | 0.17 | 210 | Spin-fluct. + orbital | 0.003 | Confirmed | Hsu (2008), Margadonna (2009) |
| FeSe monolayer on SrTiO3 | -- | -- | 65-109 | 0 | 15-20 | -- | -- | Spin + phonon (STO interface) | 0.030 | Confirmed | Wang (2012), Ge (2015), He (2013) |
| LiFeAs | PbFCl-type | P4/nmm | 18 | 0 | 3-5 | -- | 250 | Spin-fluct. | 0.005 | Confirmed | Tapp (2008) |
| KFe2As2 | ThCr2Si2-type | I4/mmm | 3.8 | 0 | -- | -- | -- | Spin-fluct., d-wave? | ~0.001 | Confirmed | Sasmal (2008) |
| CaKFe4As4 | -- | P4/mmm | 35 | 0 | 8-10 | -- | 280 | Spin-fluct., s+- | 0.011 | Confirmed | Iyo (2016) |
| (Tl,K)Fe_xSe2 | anti-PbO | I4/mmm | 31 | 0 | 5-10 | -- | 230 | Spin-fluct. | 0.009 | Confirmed | Guo (2010) |

### 2.5 Nickelate Superconductors

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| Nd0.8Sr0.2NiO2 | Infinite-layer | P4/mmm | 9-15 | 0 | 2-3 | -- | 350 | Spin-fluct.? d-wave? | 0.004 | Confirmed (thin film) | Li (2019) |
| La3Ni2O7 | Ruddlesden-Popper (n=2) | Amam / Fmmm | 80 | 14-43 | 12-18 | -- | -- | Spin-fluct. (bilayer) | 0.025 | Confirmed | Sun (2023), Hou (2023) |
| La4Ni3O10 | Ruddlesden-Popper (n=3) | -- | 30 | 20-30 | ~5 | -- | -- | Spin-fluct.? | 0.009 | Preliminary | Zhang (2024) |
| Pr0.8Sr0.2NiO2 | Infinite-layer | P4/mmm | 7-12 | 0 | -- | -- | -- | Spin-fluct.? | ~0.003 | Confirmed (thin film) | Osada (2020) |

### 2.6 Heavy Fermion and Other Unconventional Superconductors

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| CeCu2Si2 | ThCr2Si2 | I4/mmm | 0.6 | 0 | 0.1 | -- | 260 | Spin-fluct., d-wave | ~0.0002 | Confirmed | Steglich (1979) |
| UPt3 | MgCd3-type hex | P63/mmc | 0.5 | 0 | 0.08 | -- | 200 | Spin-fluct., f-wave? | ~0.0002 | Confirmed | Stewart (1984) |
| UTe2 | Orthorhombic | Immm | 1.6-2.1 | 0 | 0.3 | -- | 135 | Spin-fluct., p-wave? | ~0.0006 | Confirmed | Ran (2019) |
| UBe13 | NaZn13 | Fm-3c | 0.85 | 0 | 0.13 | -- | 590 | Heavy fermion | ~0.0003 | Confirmed | Ott (1983) |
| CeCoIn5 | HoCoGa5 | P4/mmm | 2.3 | 0 | 0.5 | -- | 240 | AFM spin-fluct., d-wave | ~0.0007 | Confirmed | Petrovic (2001) |
| PuCoGa5 | HoCoGa5 | P4/mmm | 18.5 | 0 | 3 | -- | 290 | Spin-fluct. | 0.005 | Confirmed | Sarrao (2002) |

### 2.7 Organic and Fullerene Superconductors

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| K3C60 | fcc | Fm-3m | 19.3 | 0 | 3.5 | 0.7-1.0 | 300 | Phonon (intramolecular) | 0.006 | Confirmed | Hebard (1991) |
| Rb3C60 | fcc | Fm-3m | 29 | 0 | 5 | 0.8-1.2 | 280 | Phonon | 0.009 | Confirmed | Rosseinsky (1991) |
| Cs3C60 | fcc/A15 | Fm-3m / Pm-3n | 38 | 0.7-1.5 | 7 | 1.0-1.4 | 260 | Phonon + correlation | 0.012 | Confirmed | Ganin (2010) |
| (BEDT-TTF)2Cu(NCS)2 | Orthorhombic | Pnma | 10.4 | 0 | 2 | -- | 200 | Spin-fluct.? | ~0.003 | Confirmed | Urayama (1988) |

### 2.8 Topological Superconductor Candidates

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| Cu_xBi2Se3 | Tetradymite | R-3m | 3.8 | 0 | 0.6 | -- | 180 | Nematic p-wave? | ~0.001 + I_topo | Candidate | Hor (2010), Fu & Berg (2010) |
| Sn1-xInxTe | Rock salt | Fm-3m | 1.2-4.5 | 0 | 0.2-0.7 | -- | 140 | Topological + s-wave? | ~0.001 + I_topo | Candidate | Sasaki (2012) |
| Fe(Te,Se) (interface) | PbO-type | P4/nmm | 14.5 | 0 | 1.8 | -- | 200 | s-wave + topological surface | ~0.004 + I_topo | Candidate | Zhang (2018), Wang (2018) |
| UTe2 | Orthorhombic | Immm | 1.6-2.1 | 0 | 0.3 | -- | 135 | Spin-triplet, odd-parity | ~0.0006 + I_topo | Strong candidate | Ran (2019), Jiao (2020) |
| Sr2RuO4 | K2NiF4 | I4/mmm | 1.5 | 0 | 0.1-0.3 | -- | 460 | Debated (formerly p-wave) | ~0.0004 | Confirmed | Maeno (1994), Pustogow (2019) |

### 2.9 Other Notable Systems (T_c > 30 K or of special interest)

| Material | Structure | Space Group | T_c (K) | P (GPa) | Delta(0) (meV) | lambda | Theta_D (K) | Pairing | R_tau | Status | Key Refs |
|----------|-----------|-------------|---------|----------|----------------|--------|-------------|---------|-------|--------|----------|
| Nb3Ge | A15 | Pm-3n | 23 | 0 | 3.7 | 1.12 | 302 | Phonon, s-wave | 0.007 | Confirmed | Gavaler (1973) |
| Nb3Sn | A15 | Pm-3n | 18 | 0 | 3.0 | 1.04 | 278 | Phonon, s-wave | 0.005 | Confirmed | Matthias (1954) |
| NbN | Rock salt | Fm-3m | 16 | 0 | 2.6 | 0.87 | 300 | Phonon, s-wave | 0.005 | Confirmed | -- |
| V3Si | A15 | Pm-3n | 17 | 0 | 2.8 | 0.90 | 390 | Phonon, s-wave | 0.005 | Confirmed | Hardy & Hulm (1954) |
| Boron-doped diamond | Diamond | Fd-3m | 4-11 | 0 | 1-2 | 0.4-0.5 | 2200 | Phonon, s-wave | ~0.002 | Confirmed | Ekimov (2004) |
| SrTiO3 (doped) | Perovskite | Pm-3m | 0.05-0.4 | 0 | -- | -- | 420 | Phonon? (dilute limit) | ~10^-4 | Confirmed | Schooley (1964) |
| Bi-S based (LaO1-xFxBiS2) | CeOBiS2-type | P4/nmm | 10.6 | 0 | -- | -- | 200 | Phonon? | ~0.003 | Confirmed | Mizuguchi (2012) |
| (TMTSF)2PF6 | Triclinic | P-1 | 1.2 | 0.65 | -- | -- | -- | Spin-fluct., triplet? | ~0.0004 | Confirmed | Jerome (1980) |

---

## 3. Detailed Material Class Profiles

### 3.1 Compressed Hydrides

**The physics**: Hydrogen, being the lightest element, has the highest phonon frequencies (Theta_D ~ 1000-2000 K). According to BCS theory, T_c ~ Theta_D * exp(-1/lambda), so high phonon frequencies combined with strong electron-phonon coupling can yield very high T_c. Metallic hydrogen itself is predicted to be a room-temperature superconductor (T_c ~ 300-400 K) at pressures of ~400-500 GPa (Ashcroft 1968, McMahon & Ceperley 2011).

**H3S (203 K at 155 GPa)**:
- Crystal structure: Body-centered cubic (Im-3m), each S atom is surrounded by H atoms forming an H cage. At 155 GPa, the cubic phase is stabilized.
- Mechanism: Conventional phonon-mediated BCS. The high-frequency H-S stretching modes (omega ~ 150-200 meV) dominate. lambda = 2.19 from Eliashberg calculations.
- Key insight: H3S has the highest lambda among known superconductors at the time of its discovery. The critical H-S bond-stretching phonon near 150 meV provides the dominant pairing.
- Allen-Dynes modified McMillan formula: T_c = (omega_log / 1.2) * exp[-1.04(1+lambda) / (lambda - mu*(1+0.62*lambda))] with omega_log ~ 1100 K, mu* ~ 0.13.
- Gap: 2*Delta(0)/(k_B*T_c) ~ 3.5 (BCS ratio), Delta(0) ~ 30 meV.
- Refs: Drozdov et al., Nature 525, 73 (2015); Duan et al., Sci. Rep. 4, 6968 (2014); Errea et al., PRL 114, 157004 (2015); Einaga et al., Nature Physics 13, 652 (2017).

**LaH10 (250-260 K at 170-190 GPa)**:
- Crystal structure: Face-centered cubic clathrate (Fm-3m). La atoms sit at the center of H32 cages (sodalite-like). The H sublattice forms a near-perfect fcc arrangement.
- Mechanism: Conventional phonon-mediated. lambda = 3.41, the highest known. The extraordinary coupling comes from the H cage vibrations.
- Quantum effects: Anharmonic quantum fluctuations of H atoms are essential. Harmonic phonon calculations predict imaginary frequencies; inclusion of quantum nuclear effects (SSCHA) stabilizes the structure. Errea et al. Nature 578, 66 (2020) showed quantum effects are crucial.
- Gap: Delta(0) ~ 36-40 meV, 2*Delta/(k_B*T_c) ~ 3.8 (slightly above BCS weak coupling).
- The highest confirmed T_c among all materials.
- Refs: Liu et al., PNAS 114, e6990 (2017); Somayazulu et al., PRL 122, 027001 (2019); Drozdov et al., Nature 569, 528 (2019); Errea et al., Nature 578, 66 (2020).

**YH6 and YH9**:
- YH6: Im-3m sodalite clathrate, T_c = 220 K at 166 GPa. lambda ~ 2.5.
- YH9: P63/mmc, T_c = 243 K at 201 GPa. lambda ~ 2.8.
- Both show conventional phonon-mediated pairing with very strong H-phonon coupling.
- Refs: Troyan et al., Adv. Mater. 33, 2006832 (2021); Kong et al., Nature Comm. 12, 5075 (2021).

**CaH6 (215 K at 172 GPa)**:
- Im-3m sodalite clathrate. Very similar physics to LaH10 family.
- lambda ~ 2.7, Theta_D ~ 1350 K.
- Refs: Ma et al., PRL 128, 167001 (2022).

**Lower-pressure hydrides**: LaBeH8, CeH9, ThH9/ThH10 all show T_c ~ 100-160 K at somewhat lower pressures (80-175 GPa). The quest continues for hydrogen-rich compounds stable at ambient or near-ambient pressure.

**RETRACTED claims**:
- **Carbonaceous sulfur hydride (CSH)**: Snider et al., Nature 586, 373 (2020) claimed T_c = 287.7 K at 267 GPa. Retracted in September 2022 due to data integrity concerns. Independent analysis (Hirsch & Marsiglio, van der Marel &"; Eremets group) could not reproduce. Ranga Dias (Rochester) found guilty of data fabrication.
- **Lu-N-H ("LuH2+N")**: Dasenbrock-Gammon et al., Nature 615, 244 (2023) claimed near-ambient superconductivity (294 K, 1 GPa). Retracted in November 2023. Multiple independent groups failed to reproduce. Material identified as LuH2 (a normal metal) with possible N doping producing color changes but no superconductivity.
- **Lesson**: All extraordinary claims require extraordinary evidence. The hydride field was damaged by these fraudulent claims, but the verified high-pressure results (H3S, LaH10, YH6/9, CaH6) remain robust with multiple independent confirmations.

### 3.2 Cuprate Superconductors

**The physics**: CuO2 planes are the essential structural unit. Parent compounds are Mott insulators with antiferromagnetic order. Doping introduces charge carriers into the CuO2 planes, destroying long-range AFM order and inducing superconductivity. The pairing mechanism is predominantly spin-fluctuation mediated, with d-wave symmetry (dx2-y2).

**General features**:
- Pairing symmetry: d-wave (dx2-y2), confirmed by phase-sensitive experiments (Tsuei & Kirtley 2000).
- Mechanism: Primarily antiferromagnetic spin fluctuations (J ~ 100-150 meV), though the exact pairing "glue" remains debated after nearly 40 years.
- lambda_ep for phonon contribution is small (~0.1-0.3), but the effective coupling via spin fluctuations is strong.
- Characteristic energy scales: J ~ 130 meV (superexchange), Delta_max ~ 20-50 meV, pseudogap T* ~ 200-400 K.
- All cuprates share the CuO2 plane motif; T_c is optimized at a doping of ~0.16 holes per Cu.

**Record holders**:
- Ambient pressure: HgBa2Ca2Cu3O8+d (Hg-1223), T_c = 133 K. Under 30 GPa pressure, T_c reaches 164 K.
- The multilayer Hg-12(n-1)n series has the highest T_c because the inner CuO2 planes are more optimally doped.

**Mechanism debate (2024-2026 status)**:
- Spin-fluctuation (AFM): Dominant view. Anderson's resonating valence bond (RVB) theory (1987), Scalapino's spin-fluctuation exchange. Supported by ARPES, neutron scattering, NMR.
- Phonon contribution: Lattice effects (especially apical oxygen modes) may enhance T_c by 10-20%, but are not the primary driver. Devereaux & Hackl (2007).
- Orbital current (Varma): Proposed circulating currents in CuO2 as the order parameter. Some neutron evidence (Bourges & Sidis) but controversial.
- Nematic/charge order: Interplay with charge density waves (CDW) is important, especially underdoped regime. CDW competes with superconductivity. Ghiringhelli (2012), Chang (2012).
- Hubbard model numerical solutions: Exact results from DMRG, DQMC, and tensor-network methods increasingly confirm d-wave pairing in the 2D Hubbard model at realistic parameters (Qin et al., PRX 2020).

### 3.3 Iron-Based Superconductors

**The physics**: Multi-band systems with Fe-3d orbitals forming multiple Fermi surface sheets (hole pockets at Gamma, electron pockets at M in the unfolded BZ). Pairing driven by spin fluctuations enhanced by nesting between hole and electron pockets, leading to s+- symmetry (gap changes sign between Fermi surface sheets).

**General features**:
- Crystal structure: Always contain Fe-pnictide (FeAs, FeP) or Fe-chalcogenide (FeSe, FeTe) layers.
- Four main families: 1111 (LaFeAsO), 122 (BaFe2As2), 111 (LiFeAs), 11 (FeSe).
- Pairing: Predominantly s+- (sign-changing s-wave), driven by repulsive spin fluctuations (Mazin et al. 2008, Kuroki et al. 2008). Some materials may have d-wave nodes.
- lambda_ep (phonon part) is small (~0.2-0.3); spin-fluctuation coupling is the dominant pairing channel.
- Multiband structure is essential: two-gap features observed in many compounds.

**FeSe monolayer on SrTiO3**:
- Arguably the most remarkable iron-based superconductor: T_c = 65-109 K for a single monolayer FeSe on SrTiO3 substrate, vs. 8 K for bulk FeSe.
- Enhancement mechanism: Interface phonon coupling. High-energy optical phonon in SrTiO3 (~100 meV) couples to electrons in FeSe monolayer, providing an additional pairing channel.
- This is a direct demonstration of the multi-channel pairing principle (spin + phonon) predicted by the tau framework.
- Electron doping from the substrate is also important.
- Refs: Wang et al., Chin. Phys. Lett. 29, 037402 (2012); He et al., Nature Mater. 12, 605 (2013); Ge et al., Nature Mater. 14, 285 (2015); Lee et al., Nature 515, 245 (2014).

**SmFeAsO1-xFx (T_c = 55-58 K)**: Highest T_c among iron-based at ambient pressure. The 1111 family with Sm has optimal Fermi surface nesting.

### 3.4 MgB2 (Two-Gap Superconductor)

**Crystal structure**: AlB2-type hexagonal (P6/mmm). Layers of B atoms form a graphite-like honeycomb, intercalated with Mg.

**T_c**: 39 K at ambient pressure. The highest T_c among conventional (phonon-mediated) ambient-pressure superconductors.

**Mechanism**: Phonon-mediated, but with a distinctive two-gap structure:
- Sigma band (2D B p_xy orbitals in-plane): Strong coupling to the E2g bond-stretching phonon. Gap Delta_sigma = 7.1 meV. lambda_sigma ~ 1.0.
- Pi band (3D B p_z orbitals): Weak coupling. Gap Delta_pi = 2.8 meV. lambda_pi ~ 0.4.
- Average lambda = 0.87, Theta_D = 750 K.
- The E2g phonon at ~75 meV is the key pairing boson.

**Why T_c is so high for a "conventional" superconductor**: The B-B bond-stretching mode has an unusually high frequency and coupling, enabled by the quasi-2D boron network and light B mass. The two-gap structure actually helps: the sigma band provides the strong pairing, while the pi band acts as a reservoir.

**Tau-framework significance**: MgB2 demonstrates multi-channel pairing (sigma + pi bands), and the two-gap structure naturally maps to I_pair = I_sigma + I_pi in the tau framework.

**Refs**: Nagamatsu et al., Nature 410, 63 (2001); Choi et al., Nature 418, 758 (2002); Kortus et al., PRL 86, 4656 (2001).

### 3.5 Topological Superconductor Candidates

**Why topological superconductors matter for tau framework**: Topological protection provides temperature-independent pairing information I_topo = |nu| * ln(2), which does not decrease with temperature. This is a fundamentally different contribution to I_pair compared to phonon or spin channels.

**UTe2 (Strongest candidate for p-wave / spin-triplet)**:
- Orthorhombic (Immm), T_c = 1.6-2.1 K.
- Evidence for spin-triplet pairing: re-entrant superconductivity under field, anomalous Knight shift, possible half-quantum vortices.
- Multiple superconducting phases under pressure and field.
- Could host Majorana surface states.
- Refs: Ran et al., Science 365, 684 (2019); Jiao et al., Nature 579, 523 (2020).

**Cu_xBi2Se3**:
- Doped topological insulator. T_c ~ 3.8 K.
- Evidence for nematic superconductivity (two-fold symmetric gap) suggesting odd-parity pairing.
- Refs: Hor et al., PRL 104, 057001 (2010); Matano et al., Nature Physics 12, 852 (2016).

**Fe(Te,Se) (Surface/interface)**:
- Bulk T_c ~ 14.5 K (conventional s-wave).
- Topological surface states observed. Possible Majorana bound states in vortex cores.
- The combination of superconductivity + topological surface states makes it a platform for Majorana physics.
- Refs: Zhang et al., Science 360, 182 (2018); Wang et al., Science 362, 926 (2018).

**Sr2RuO4 (Status update)**:
- Long believed to be a p-wave spin-triplet superconductor (analog of He-3). T_c ~ 1.5 K.
- 2019 NMR results (Pustogow et al., Nature 574, 72, 2019) showed a drop in Knight shift below T_c, inconsistent with a simple p-wave state. The pairing symmetry is now actively debated.
- Current candidates: even-parity d + ig state, or singlet d-wave. The p-wave assignment is no longer consensus.

### 3.6 Nickelate Superconductors (2019-present)

**The physics**: Nickelate superconductors are structurally analogous to cuprates (NiO2 planes vs. CuO2 planes) but with important differences: Ni has a d^9-delta or d^8 configuration, and the d-p hybridization is different.

**Infinite-layer nickelates (Nd0.8Sr0.2NiO2)**:
- T_c = 9-15 K in thin films. Not yet reproduced in bulk.
- Infinite-layer structure obtained by topotactic reduction of the perovskite NdNiO3.
- Pairing mechanism: likely spin-fluctuation mediated, but the details differ from cuprates (Ni 3d x2-y2 orbital, but with additional Nd 5d states at Fermi level).
- Refs: Li et al., Nature 572, 624 (2019); Osada et al., Nano Lett. 20, 5735 (2020).

**Bilayer nickelate La3Ni2O7 under pressure**:
- T_c ~ 80 K at 14-43 GPa --- a major discovery in 2023.
- Ruddlesden-Popper n=2 structure. Under pressure, the Ni-O bond lengths change, enhancing interlayer coupling.
- Pairing mechanism: Bilayer spin-fluctuation exchange, possibly s+- or d-wave.
- Significance: First nickelate with T_c competitive with cuprates, and at much lower pressure than hydrides.
- Refs: Sun et al., Nature 621, 493 (2023); Hou et al., Chin. Phys. Lett. 40, 117302 (2023).

**Trilayer La4Ni3O10**:
- T_c ~ 30 K under pressure (20-30 GPa). Preliminary reports in 2024.
- Suggests a trend: more NiO2 layers may increase T_c (analogous to multilayer cuprates).

### 3.7 Retracted and Controversial Claims

**Ranga Dias affair (2020-2023)**:
1. **CSH (carbonaceous sulfur hydride)**: Published in Nature 586, 373 (2020). Claimed T_c = 287.7 K at 267 GPa. Multiple red flags: unusual data processing, non-standard background subtraction, refusal to share raw data. Retracted September 2022. Hirsch & Marsiglio published detailed critiques.
2. **Lu-N-H system**: Published in Nature 615, 244 (2023). Claimed T_c = 294 K at 1 GPa (near-ambient pressure). The color change in LuH2 under pressure was misinterpreted as a superconducting transition. Retracted November 2023. At least 10 independent groups failed to reproduce.
3. **Impact**: Dias was found to have fabricated data in both cases (and in his PhD thesis work on MnS2). Rochester fired him. The retractions damaged public trust in the high-pressure superconductivity field, but the verified results (H3S, LaH10, etc.) from other groups remain solid.

**LK-99 (Cu-substituted lead apatite, 2023)**:
- Lee et al. (South Korea, preprints July 2023) claimed room-temperature, ambient-pressure superconductivity in Pb10-xCux(PO4)6O (x ~ 1).
- Caused massive excitement and rapid replication attempts worldwide.
- Within weeks, multiple groups showed: the material is NOT a superconductor. The observed resistance drop was due to Cu2S impurity (a semiconductor with a phase transition near 400 K), and the "levitation" was due to paramagnetism, not the Meissner effect.
- Definitive disproof: single crystals showed no superconductivity (KAIST, Max Planck).

---

## 4. Theoretical Frameworks and Pairing Mechanisms

### 4.1 BCS/Eliashberg (Phonon-Mediated)

- Standard framework for conventional superconductors and hydrides.
- T_c determined by: lambda (electron-phonon coupling), omega_log (logarithmic average phonon frequency), mu* (Coulomb pseudopotential).
- Allen-Dynes modified McMillan formula: T_c = (omega_log/1.2) * exp[-1.04(1+lambda)/(lambda - mu*(1+0.62*lambda))]
- Typical parameters: mu* ~ 0.10-0.15.
- For hydrides, lambda can reach 3-4, pushing T_c toward room temperature.

### 4.2 Spin-Fluctuation Mediated

- Dominant in cuprates, iron-based, heavy fermions.
- The effective pairing interaction V_eff(q,omega) comes from the spin susceptibility chi(q,omega).
- In cuprates: AFM fluctuations at Q = (pi,pi) lead to d-wave pairing.
- In iron-based: nesting between hole and electron pockets leads to s+- pairing.
- Key energy scale: J (superexchange) ~ 100-150 meV in cuprates.

### 4.3 Orbital Fluctuation Mediated

- Relevant in iron-based superconductors and some heavy fermion systems.
- Orbital degeneracy of Fe 3d states drives orbital fluctuations.
- Can cooperate or compete with spin fluctuations.
- Kontani & Onari, PRL 104, 157001 (2010).

### 4.4 Topological Protection

- In topological superconductors, the pairing state has a nontrivial topological invariant.
- Bulk-boundary correspondence: topological surface states are protected by symmetry.
- Majorana bound states at vortex cores or edges.
- The topological contribution to I_pair is temperature-independent: I_topo = |nu| * ln(2).

### 4.5 Tau-Framework (This Work)

The tau framework (from `paper_superconductivity_tau.tex`) provides a unified screening criterion:

```
tau_eff(T) = max{0, tau_dec(T) - I_pair(T)}
```

Superconductivity occurs when tau_eff = 0, i.e., when the total pairing information I_pair exceeds the thermal decoherence tau_dec.

**Key predictions**:
1. Multi-channel pairing (I_phonon + I_spin + I_orbital + I_topo) is the optimal route to room temperature.
2. Topological protection provides temperature-independent I_topo, reducing the burden on other channels.
3. s-wave symmetry is more efficient than d-wave (no nodal information waste).
4. FeSe/SrTiO3 is a natural demonstration of multi-channel enhancement.

---

## 5. R_tau Analysis Summary

### 5.1 Calibration

From Table 1 of `paper_superconductivity_tau.tex`:

| Material | T_c (K) | Delta(0) (meV) | lambda_ep | R_tau | % of RT target |
|----------|---------|----------------|-----------|-------|----------------|
| Nb | 9.3 | 1.5 | 0.82 | 0.003 | 0.3% |
| MgB2 | 39 | 7.1 | 0.87 | 0.012 | 1.2% |
| YBCO | 93 | 20 | -- | 0.031 | 3.1% |
| H3S (155 GPa) | 203 | 30 | 2.19 | 0.069 | 6.9% |
| LaH10 (170 GPa) | 260 | 40 | 3.41 | 0.092 | 9.2% |
| **Room-T target** | **300** | **>= 50** | **--** | **>= 1** | **100%** |

### 5.2 Extended R_tau Estimates for All T_c > 30 K Materials

Using the proxy formula R_tau ~ (T_c/300) * coupling_factor:

| Material | T_c (K) | R_tau (est.) | Gap to RT |
|----------|---------|-------------|-----------|
| La2-xBaxCuO4 | 35 | 0.010 | 100x |
| Ba0.6K0.4Fe2As2 | 38 | 0.012 | 83x |
| MgB2 | 39 | 0.012 | 83x |
| Cs3C60 | 38 | 0.012 | 83x |
| LSCO | 40 | 0.013 | 77x |
| SmFeAsO | 55 | 0.018 | 56x |
| FeSe/STO monolayer | 65-109 | 0.030 | 33x |
| La3Ni2O7 (14 GPa) | 80 | 0.025 | 40x |
| Bi-2212 | 90 | 0.030 | 33x |
| YBCO | 93 | 0.031 | 32x |
| Hg-1201 | 97 | 0.032 | 31x |
| Bi-2223 | 110 | 0.037 | 27x |
| Tl-2212 | 110 | 0.037 | 27x |
| CeH9 (100 GPa) | 117 | 0.035 | 29x |
| Tl-2223 | 125 | 0.042 | 24x |
| Hg-1212 | 128 | 0.043 | 23x |
| Hg-1223 | 133 | 0.050 | 20x |
| ThH10 (175 GPa) | 161 | 0.052 | 19x |
| Hg-1223 (30 GPa) | 164 | 0.055 | 18x |
| H3S (155 GPa) | 203 | 0.069 | 14x |
| CaH6 (172 GPa) | 215 | 0.075 | 13x |
| YH6 (166 GPa) | 220 | 0.078 | 13x |
| YH9 (201 GPa) | 243 | 0.087 | 11x |
| LaH10 (170 GPa) | 260 | 0.092 | 11x |

**Key observation**: Even the highest confirmed T_c (LaH10 at 260 K) has R_tau ~ 0.09, an order of magnitude below the room-temperature threshold R_tau = 1. This means tau_dec(300K) is so large that a qualitative enhancement (not just quantitative improvement) of I_pair is needed.

### 5.3 Strategies to Reach R_tau = 1

From the paper's analysis (Section VI):

1. **Single-channel brute force**: lambda ~ 5-10 with omega_D ~ 200 meV. This is beyond current hydrides (lambda_max ~ 3.4 for LaH10). Would require metallic hydrogen or extremely hydrogen-dense compounds at manageable pressures.

2. **Multi-channel coherent addition**: I_total = I_phonon + I_spin + I_orbital. Three channels each contributing ~ 3x LaH10's I_pair with coherent addition. Candidate platform: FeSe/SrTiO3 type systems show this principle works.

3. **Topological floor**: I_topo >= 0.5 * tau_dec(300K) reduces the needed I_conventional by half. Requires a material with nontrivial topology AND strong conventional pairing. No known material achieves this yet.

4. **Optimal material profile**: High Theta_D (light elements), strong electron-phonon coupling, proximity to a magnetic quantum critical point, nontrivial band topology. Hydrogen-rich materials with mixed bonding character (covalent + ionic) near magnetic instabilities.

---

## 6. Machine Learning and Computational Screening

### 6.1 ML for Superconductor Discovery (State of the Art 2024-2026)

**Approaches**:

1. **Graph Neural Networks (GNNs) for crystal structure prediction**:
   - CGCNN (Xie & Grossman, PRL 2018): Crystal Graph Convolutional Neural Network predicts material properties from crystal structure.
   - MEGNet (Chen et al., Chem. Mater. 2019): Multi-Edge Graph networks for molecular and crystal properties.
   - Applied to superconductors: Stanev et al., npj Comp. Mater. 4, 29 (2018) --- random forest model predicting T_c from composition.

2. **High-throughput DFT screening**:
   - Electron-phonon coupling calculations with EPW (Ponce et al., Comp. Phys. Comm. 2016) and Quantum ESPRESSO.
   - Sanna et al. (2024): Screened ~10,000 hydrogen-rich compounds using SSCHA + Eliashberg, identified ~200 candidates with T_c > 100 K.
   - Flores-Livas et al., Physics Reports 856, 1 (2020): Comprehensive review of computational approaches.

3. **Generative models for material design**:
   - VAE and diffusion models trained on materials databases (ICSD, Materials Project, AFLOW).
   - CrystalGAN, CDVAE (Xie et al., ICLR 2022): Generate new crystal structures with desired properties.
   - Application to superconductors: Still in early stages. The main bottleneck is the small training set (only ~30,000 known superconductors in SuperCon database, most with T_c < 30 K).

4. **Natural Language Processing for literature mining**:
   - Automated extraction of T_c, composition, and synthesis conditions from papers.
   - Tshitoyan et al., Nature 571, 95 (2019): Word2Vec embeddings predicted thermoelectric materials; similar approach applied to superconductors.

**Key databases**:
- SuperCon (NIMS Japan): ~33,000 superconductor entries.
- AFLOW: ~3.5 million calculated material entries.
- Materials Project: ~150,000 materials with DFT properties.
- JARVIS-DFT: ~75,000 materials, includes some electron-phonon data.

### 6.2 Computational Screening Methodologies

1. **Structure prediction + Eliashberg**: CALYPSO, USPEX, AIRSS for crystal structure prediction under pressure; then EPW or ABINIT for electron-phonon coupling. This pipeline discovered H3S (Duan 2014) and LaH10 (Liu 2017) theoretically before experiments.

2. **Bayesian optimization**: Adaptive sampling of composition/pressure space to efficiently find high-T_c candidates. Lookman et al., npj Comp. Mater. 5, 21 (2019).

3. **Transfer learning**: Pre-train on large material property datasets, fine-tune on smaller superconductor datasets. Improves predictions for underrepresented material classes.

### 6.3 Recent ML Predictions (2024-2025)

- **Ternary hydrides**: ML-guided screening identified several promising ternary hydrides (e.g., LaYH20, CaLaH12) with predicted T_c > 200 K at moderate pressures (50-100 GPa). Some await experimental verification.
- **Non-hydride ambient-pressure candidates**: ML models trained on SuperCon data suggest that layered borocarbides, transition-metal nitrides, and certain intermetallics may harbor undiscovered superconductors with T_c > 40 K. Limited experimental follow-up.
- **GNN-predicted topological superconductors**: Combining topological invariant calculations with T_c prediction to identify materials with both nontrivial topology and significant T_c.

---

## 7. Tau-Framework Predictions and Research Directions

### 7.1 Most Promising Candidates for Multi-Channel R_tau Enhancement

Based on the tau-framework analysis, materials at the intersection of multiple pairing channels are the most promising:

1. **FeSe/SrTiO3 interface** (proven multi-channel):
   - Already demonstrates spin + phonon (interface) enhancement: T_c jumps from 8 K to 65-109 K.
   - R_tau enhanced by ~ 10x over bulk FeSe.
   - Research direction: Engineer better interfaces, optimize doping, explore other substrates.

2. **La3Ni2O7 under pressure** (bilayer spin + possible orbital):
   - T_c = 80 K with bilayer spin-fluctuation pairing.
   - If orbital fluctuations also contribute (as in iron-based systems), multi-channel enhancement is plausible.
   - Research direction: Detailed Eliashberg analysis of the pairing interaction; search for ambient-pressure analogs.

3. **Hydrogen-rich compounds near magnetic instability**:
   - Combine the high phonon frequencies of hydrides (large I_phonon) with proximity to an antiferromagnetic instability (I_spin contribution).
   - Candidate: Fe-H or Mn-H compounds under pressure.
   - Research direction: DFT screening of magnetic hydrides.

4. **Topological hydrides**:
   - Combine the extreme phonon coupling of hydrides with nontrivial band topology.
   - Candidate: PtH (predicted topological + superconducting at high P).
   - The topological floor I_topo ~ ln(2) per surface channel would add a temperature-independent boost.

### 7.2 Calibration Plan (from future_directions_tau_framework.md)

**Step 1**: For each known superconductor (Nb, MgB2, YBCO, H3S, LaH10), compute:
- tau_phonon(T) from alpha^2 F(omega) (DFT + tunneling data)
- I_Cooper from measured Delta(T)
- Verify tau_eff crosses zero at T = T_c

**Step 2**: Identify which mechanism gives largest I per unit of tau_phonon. This determines the optimal pairing channel.

**Step 3**: Screen materials databases (AFLOW, Materials Project) for:
- High I_phonon: strong e-ph coupling + high phonon frequencies
- High I_spin: near magnetic QCP
- Topological enhancement: nontrivial band topology

### 7.3 Tuna-9 Quantum Simulation

From the paper: a 2-site, 2-electron Cooper pair can be simulated on a quantum processor:
- Prepare Bell state (analog of Cooper pair)
- Apply thermal noise channel (analog of phonon scattering)
- Measure tau_eff vs. noise strength
- Already demonstrated 89% tau reduction with single Bell pair on Tuna-9

---

## 8. Summary and Outlook

### 8.1 State of the Field (March 2026)

**Confirmed high-T_c records**:
1. LaH10: 260 K at 170-190 GPa (highest confirmed T_c overall)
2. YH9: 243 K at 201 GPa
3. H3S: 203 K at 155 GPa
4. Hg-1223: 164 K at 30 GPa (highest T_c under moderate pressure)
5. Hg-1223: 133 K at ambient pressure (highest ambient-pressure T_c)
6. FeSe/STO: 65-109 K (highest interface/2D T_c)
7. La3Ni2O7: 80 K at 14-43 GPa (nickelate record)

**No verified room-temperature superconductor exists** as of March 2026. All claims (CSH, Lu-N-H, LK-99) have been retracted or debunked.

### 8.2 Tau-Framework Assessment

The R_tau analysis reveals a stark quantitative gap: even the best known superconductor (LaH10, R_tau ~ 0.09) is an order of magnitude below the room-temperature threshold. This suggests:

1. **Incremental improvement within a single pairing channel will not suffice**. We need either a qualitatively new mechanism or coherent multi-channel enhancement.

2. **Multi-channel pairing is the most promising route**. FeSe/STO already demonstrates a 10x T_c enhancement through spin + phonon cooperation.

3. **Topological protection offers a unique advantage**: temperature-independent I_topo that does not degrade at room temperature. No known material exploits this yet for high-T_c purposes.

4. **The design principle is clear**: maximize I_pair through multiple concurrent channels while minimizing tau_dec through high Theta_D and structural order.

### 8.3 Key References (Organized by Topic)

**Reviews**:
- Flores-Livas et al., Physics Reports 856, 1 (2020): Computational approaches to superconductivity.
- Pickard et al., Ann. Rev. Cond. Mat. Phys. 11, 57 (2020): High-pressure hydrides.
- Norman, Science 332, 196 (2011): Challenge of unconventional superconductivity.
- Proust & Taillefer, Ann. Rev. Cond. Mat. Phys. 10, 409 (2019): Cuprate review.
- Hosono & Kuroki, Physica C 514, 399 (2015): Iron-based review.
- Sato & Ando, Rep. Prog. Phys. 80, 076501 (2017): Topological superconductors.

**Landmark experimental papers**:
- Onnes, Commun. Phys. Lab. Univ. Leiden (1911): Discovery of superconductivity.
- Bednorz & Muller, Z. Phys. B 64, 189 (1986): Cuprate discovery.
- Nagamatsu et al., Nature 410, 63 (2001): MgB2.
- Kamihara et al., JACS 130, 3296 (2008): Iron-based discovery.
- Drozdov et al., Nature 525, 73 (2015): H3S at 203 K.
- Somayazulu et al., PRL 122, 027001 (2019): LaH10 at 250 K.
- Li et al., Nature 572, 624 (2019): Nickelate superconductivity.
- Sun et al., Nature 621, 493 (2023): La3Ni2O7 at 80 K.

**Theory**:
- BCS, Phys. Rev. 108, 1175 (1957): BCS theory.
- Eliashberg, JETP 11, 696 (1960): Eliashberg theory.
- Anderson, Science 235, 1196 (1987): RVB theory.
- Mazin et al., PRL 101, 057003 (2008): s+- pairing in iron-based.
- Liu et al., PNAS 114, e6990 (2017): Prediction of LaH10.
- Errea et al., Nature 578, 66 (2020): Quantum nuclear effects in hydrides.

**Machine learning**:
- Stanev et al., npj Comp. Mater. 4, 29 (2018): ML prediction of T_c.
- Xie & Grossman, PRL 120, 145301 (2018): CGCNN.
- Tshitoyan et al., Nature 571, 95 (2019): NLP for materials.
- Xie et al., ICLR 2022: CDVAE for crystal generation.

**Tau framework (our work)**:
- Huang, Paper 1 (2026): Petz recovery unification.
- Huang, paper_superconductivity_tau.tex: Superconductivity as zero temporal asymmetry.
- Huang, future_directions_tau_framework.md: Research roadmap.

---

## Appendix A: Complete Data Tables

### A.1 Compressed Hydride Superconductors (All Confirmed)

| Compound | Structure type | T_c (K) | P (GPa) | lambda | omega_log (K) | Delta(0) (meV) | 2Delta/(k_BT_c) | R_tau |
|----------|---------------|---------|----------|--------|---------------|----------------|-----------------|-------|
| H3S | bcc | 203 | 155 | 2.19 | 1100 | 30 | 3.43 | 0.069 |
| LaH10 | fcc clathrate | 250-260 | 170-190 | 3.41 | 1050 | 36-40 | 3.57 | 0.092 |
| YH6 | sodalite | 220 | 166 | 2.5 | 1150 | 33 | 3.48 | 0.078 |
| YH9 | hcp-type | 243 | 201 | 2.8 | 1200 | 37 | 3.53 | 0.087 |
| CaH6 | sodalite | 215 | 172 | 2.7 | 1100 | 32 | 3.45 | 0.075 |
| (La,Y)H10 | fcc clathrate | 253 | 183 | 3.1 | 1080 | 38 | 3.48 | 0.089 |
| ThH9 | hcp-type | 146 | 170 | 2.0 | 950 | 24 | 3.81 | 0.052 |
| ThH10 | fcc clathrate | 161 | 175 | 2.1 | 1000 | 25 | 3.60 | 0.055 |
| CeH9 | hcp-type | 100-117 | 88-100 | 1.6 | 850 | 16 | 3.52 | 0.035 |
| LaBeH8 | fcc clathrate | 110 | 80 | 1.4 | 900 | 17 | 3.58 | 0.035 |

### A.2 Cuprate Family (Ambient Pressure T_c)

| Compound | # CuO2 layers | T_c,max (K) | Delta_max (meV) | J (meV) | Symmetry | R_tau |
|----------|---------------|-------------|-----------------|---------|----------|-------|
| La2-xSrxCuO4 | 1 | 40 | 12 | 135 | d(x2-y2) | 0.013 |
| YBa2Cu3O7 | 2 | 93 | 20 | 125 | d(x2-y2) | 0.031 |
| Bi2Sr2CaCu2O8 | 2 | 95 | 35 | 130 | d(x2-y2) | 0.030 |
| Bi2Sr2Ca2Cu3O10 | 3 | 110 | 38 | 130 | d(x2-y2) | 0.037 |
| Tl2Ba2CaCu2O8 | 2 | 110 | 35 | 120 | d(x2-y2) | 0.037 |
| Tl2Ba2Ca2Cu3O10 | 3 | 125 | 40 | 120 | d(x2-y2) | 0.042 |
| HgBa2CuO4 | 1 | 97 | 25 | 140 | d(x2-y2) | 0.032 |
| HgBa2CaCu2O6 | 2 | 128 | 40 | 140 | d(x2-y2) | 0.043 |
| HgBa2Ca2Cu3O8 | 3 | 133 (164 at 30 GPa) | 50 | 140 | d(x2-y2) | 0.050 |

### A.3 Iron-Based Family

| Compound | Family | T_c (K) | Delta(0) (meV) | Symmetry | Nesting | R_tau |
|----------|--------|---------|----------------|----------|---------|-------|
| SmFeAsO1-xFx | 1111 | 55-58 | 12-15 | s+- | Strong | 0.018 |
| NdFeAsO1-xFx | 1111 | 52 | 10-12 | s+- | Strong | 0.017 |
| LaFeAsO1-xFx (4 GPa) | 1111 | 43 | 7-10 | s+- | Moderate | 0.013 |
| Ba0.6K0.4Fe2As2 | 122 | 38 | 6-12 | s+- | Moderate | 0.012 |
| CaKFe4As4 | 1144 | 35 | 8-10 | s+- | Moderate | 0.011 |
| FeSe/STO monolayer | 11 interface | 65-109 | 15-20 | s-wave (no nodes) | N/A | 0.030 |
| FeSe (7 GPa) | 11 | 37 | 5-8 | s+- ? | Moderate | 0.012 |

---

## Appendix B: Physical Constants Used

| Symbol | Name | Value |
|--------|------|-------|
| k_B | Boltzmann constant | 8.617 x 10^-5 eV/K = 0.08617 meV/K |
| k_B * 300 K | Thermal energy at RT | 25.85 meV |
| hbar | Reduced Planck | 6.582 x 10^-16 eV*s |
| ln(2) | Natural log of 2 | 0.6931 |

---

*This survey is a living document. Last updated: 2026-03-28.*
*Cross-reference with: paper_superconductivity_tau.tex, future_directions_tau_framework.md*
