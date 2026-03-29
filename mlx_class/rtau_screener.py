#!/usr/bin/env python3
"""
R_tau Superconductor Screening Tool
====================================

Screens for potential high-Tc superconductor candidates using the R_tau criterion
from the Petz recovery / retrodictability framework.

Key idea:
  R_tau = I_Cooper(max) / tau_phonon(T)
  Proxy: R_tau ~ (T_c / T) * (Theta_D / T)

Multi-channel extension:
  R_tau_multi = alpha * R_tau_single
  where alpha = number of independent pairing channels

Usage:
  python rtau_screener.py                        # built-in dataset
  python rtau_screener.py --api-key YOUR_KEY     # Materials Project API
  python rtau_screener.py --custom               # interactive custom input
  python rtau_screener.py --output results.csv   # save CSV

Requires: numpy (mandatory), mp-api (optional, for Mode 1)
"""

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy is required. Install with: pip install numpy")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Material:
    """A candidate superconductor material."""
    material_id: str
    formula: str
    space_group: str = "unknown"
    theta_D: float = 0.0            # Debye temperature (K)
    lambda_ep: float = 0.0          # electron-phonon coupling
    lambda_source: str = "estimate"  # how lambda was obtained
    Tc_experimental: float = 0.0     # experimental T_c if known (K)
    Tc_estimated: float = 0.0        # McMillan estimate (K)
    band_gap: float = 0.0            # eV
    e_above_hull: float = 0.0        # eV/atom
    density: float = 0.0             # g/cm^3
    pressure_GPa: float = 0.0        # synthesis pressure
    channels: List[str] = field(default_factory=lambda: ["phonon"])
    notes: str = ""
    # Computed
    R_tau: float = 0.0
    R_tau_multi: float = 0.0
    channels_needed_300K: int = 0


# ---------------------------------------------------------------------------
# Physics: McMillan / Allen-Dynes T_c
# ---------------------------------------------------------------------------

def mcmillan_Tc(theta_D: float, lambda_ep: float, mu_star: float = 0.13) -> float:
    """
    McMillan formula for T_c.

    T_c = (theta_D / 1.45) * exp[-1.04*(1+lambda) / (lambda - mu*(1+0.62*lambda))]

    Returns T_c in Kelvin.  Returns 0 if denominator is non-positive
    (no superconductivity).
    """
    if lambda_ep <= 0 or theta_D <= 0:
        return 0.0
    denom = lambda_ep - mu_star * (1.0 + 0.62 * lambda_ep)
    if denom <= 0:
        return 0.0
    exponent = -1.04 * (1.0 + lambda_ep) / denom
    if exponent < -50:
        return 0.0
    return (theta_D / 1.45) * math.exp(exponent)


def allen_dynes_Tc(theta_D: float, lambda_ep: float,
                   mu_star: float = 0.13, omega_log: Optional[float] = None) -> float:
    """
    Allen-Dynes modification of the McMillan formula.

    Uses omega_log if provided; otherwise falls back to theta_D / 1.2 as proxy
    for omega_log.

    T_c = (omega_log / 1.2) * exp[-1.04*(1+lambda) / (lambda - mu*(1+0.62*lambda))]
    with strong-coupling correction factors f1, f2.
    """
    if lambda_ep <= 0 or theta_D <= 0:
        return 0.0

    if omega_log is None:
        omega_log = theta_D / 1.2   # rough proxy

    denom = lambda_ep - mu_star * (1.0 + 0.62 * lambda_ep)
    if denom <= 0:
        return 0.0

    # Strong-coupling corrections (Allen & Dynes 1975)
    Lambda_1 = 2.46 * (1.0 + 3.8 * mu_star)
    Lambda_2 = 1.82 * (1.0 + 6.3 * mu_star) * (omega_log / theta_D)

    f1 = (1.0 + (lambda_ep / Lambda_1) ** (3.0 / 2.0)) ** (1.0 / 3.0)
    f2_arg = 1.0 + (lambda_ep ** 2 / (Lambda_2 ** 2 + lambda_ep ** 2))
    # Prevent overflow
    f2 = 1.0 + (lambda_ep ** 2) / (lambda_ep ** 2 + Lambda_2 ** 2) * \
         ((1.0 - lambda_ep * omega_log / theta_D) if lambda_ep * omega_log / theta_D < 1 else 0.0)
    f2 = max(f2, 0.1)

    exponent = -1.04 * (1.0 + lambda_ep) / denom
    if exponent < -50:
        return 0.0

    return f1 * f2 * (omega_log / 1.2) * math.exp(exponent)


# ---------------------------------------------------------------------------
# R_tau computation
# ---------------------------------------------------------------------------

def compute_R_tau(Tc: float, theta_D: float, T: float = 300.0) -> float:
    """
    R_tau proxy:
      R_tau = (T_c / T) * (Theta_D / T)

    This is a dimensionless measure of the Cooper pairing strength
    relative to thermal decoherence at temperature T.

    Physical meaning:
      - I_Cooper(max) ~ T_c * Theta_D  (pairing interaction * phonon frequency)
      - tau_phonon(T) ~ T^2             (thermal scattering rate)
      - R_tau = I_Cooper / tau_phonon ~ T_c * Theta_D / T^2
    """
    if T <= 0:
        return float('inf')
    return (Tc * theta_D) / (T * T)


def compute_R_tau_at_temperatures(Tc: float, theta_D: float,
                                  temps: Optional[List[float]] = None) -> dict:
    """Compute R_tau at several temperatures."""
    if temps is None:
        temps = [273.15, 293.15, 310.15, 363.15]  # 0C, 20C, 37C, 90C
    return {T: compute_R_tau(Tc, theta_D, T) for T in temps}


def lambda_needed_for_Rtau1(theta_D: float, T: float = 300.0,
                            mu_star: float = 0.13) -> float:
    """
    What electron-phonon coupling lambda is needed so that R_tau = 1 at
    temperature T?

    R_tau = 1  =>  T_c = T^2 / theta_D
    Then invert McMillan:
      T_c = (theta_D/1.45) exp[-1.04(1+lambda)/(lambda - mu*(1+0.62*lambda))]
    """
    Tc_target = T * T / theta_D
    if Tc_target <= 0 or theta_D <= 0:
        return float('inf')

    ratio = Tc_target / (theta_D / 1.45)
    if ratio >= 1.0:
        # Need very strong coupling; McMillan breaks down.
        # Return approximate: lambda ~ 3+ for extreme cases
        return 3.0 + ratio

    if ratio <= 0:
        return float('inf')

    ln_ratio = math.log(ratio)
    if ln_ratio >= 0:
        return 3.0

    # Invert: ln_ratio = -1.04*(1+lambda)/(lambda - mu*(1+0.62*lambda))
    # Let x = lambda
    # ln_ratio * (x - mu*(1+0.62*x)) = -1.04*(1+x)
    # ln_ratio * x - ln_ratio*mu - ln_ratio*0.62*mu*x = -1.04 - 1.04*x
    # x*(ln_ratio - 0.62*mu*ln_ratio + 1.04) = -1.04 + mu*ln_ratio
    # x = (-1.04 + mu*ln_ratio) / (ln_ratio - 0.62*mu*ln_ratio + 1.04)

    A = ln_ratio * (1.0 - 0.62 * mu_star) + 1.04
    B = mu_star * ln_ratio - 1.04

    if abs(A) < 1e-12:
        return float('inf')

    lam = B / A
    if lam < 0:
        return float('inf')
    return lam


def channels_needed(R_tau_single: float, target_Tc: float = 300.0,
                    T: float = 300.0, theta_D: float = 300.0) -> int:
    """
    How many independent pairing channels alpha are needed so that
    R_tau_multi = alpha * R_tau_single >= R_tau_target?

    R_tau_target = target_Tc * theta_D / T^2
    """
    if R_tau_single <= 0:
        return 999
    R_tau_target = (target_Tc * theta_D) / (T * T)
    return max(1, math.ceil(R_tau_target / R_tau_single))


# ---------------------------------------------------------------------------
# Empirical lambda estimator
# ---------------------------------------------------------------------------

def estimate_lambda(formula: str, theta_D: float, notes: str = "") -> Tuple[float, str]:
    """
    Estimate electron-phonon coupling lambda from composition heuristics.

    Returns (lambda, source_description).
    """
    f = formula.upper()

    # Known materials with measured lambda
    known = {
        "H3S":    (2.2,  "exp: Drozdov+2015"),
        "LAH10":  (2.6,  "calc: Peng+2017"),
        "YH6":    (2.0,  "calc: Troyan+2021"),
        "YH9":    (2.5,  "calc: Kong+2021"),
        "CAH6":   (2.2,  "calc: Wang+2012"),
        "THH10":  (2.4,  "calc: Semenok+2020"),
        "PRH9":   (1.8,  "calc: Zhou+2020"),
        "CEH9":   (1.7,  "calc: Salke+2019"),
        "BAH12":  (2.1,  "calc: theory"),
        "MGH6":   (1.6,  "calc: Feng+2015"),
        "MGB2":   (0.87, "exp: Choi+2002"),
        "NB":     (1.04, "exp: measured"),
        "NB3SN":  (1.8,  "exp: measured"),
        "NB3GE":  (1.7,  "exp: measured"),
        "PBMO6S8": (1.3, "exp: Chevrel"),
        "PB":     (1.55, "exp: measured"),
        "V3SI":   (1.3,  "exp: measured"),
        "K3C60":  (0.9,  "exp: fulleride"),
        "CS3C60": (1.2,  "exp: fulleride"),
    }

    # Normalize formula for lookup
    f_clean = f.replace(" ", "").replace("_", "")
    for key, (lam, src) in known.items():
        if f_clean == key or f_clean.replace("(", "").replace(")", "") == key:
            return lam, src

    # Heuristic classification
    has_H = "H" in f
    h_count = 0
    # Simple H counting
    import re
    h_matches = re.findall(r'H(\d*)', formula)
    for m in h_matches:
        h_count += int(m) if m else 1

    has_light = any(el in f for el in ["B", "C", "N", "O"])

    # High-pressure hydrides (many H atoms)
    if has_H and h_count >= 6:
        lam = 1.8 + 0.1 * min(h_count - 6, 6)  # 1.8 to 2.4
        return lam, f"hydride heuristic (H_count={h_count})"

    # Medium hydrides
    if has_H and h_count >= 3:
        lam = 1.2 + 0.15 * (h_count - 3)
        return lam, f"hydride heuristic (H_count={h_count})"

    # Light-H compounds
    if has_H and has_light:
        lam = 0.8 + 0.1 * h_count
        return lam, "light-element hydride heuristic"

    # Metal hydrides
    if has_H:
        lam = 0.6 + 0.15 * h_count
        return lam, "metal hydride heuristic"

    # Borides / carbides / nitrides
    if "B" in f:
        return 0.7, "boride heuristic"
    if "C" in f and "60" in formula:  # fullerides
        return 1.0, "fulleride heuristic"
    if "C" in f:
        return 0.6, "carbide heuristic"
    if "N" in f:
        return 0.5, "nitride heuristic"

    # Cuprates (not really BCS, but for comparison)
    if "CU" in f and "O" in f:
        return 1.5, "cuprate (non-BCS, effective lambda)"

    # Nickelates
    if "NI" in f and "O" in f:
        return 1.2, "nickelate (non-BCS, effective lambda)"

    # Heavy fermion
    if any(el in f for el in ["CE", "U", "YB", "PR"]):
        return 0.4, "heavy-fermion heuristic"

    # Default metallic
    return 0.5, "generic metal heuristic"


# ---------------------------------------------------------------------------
# Pairing channel identification
# ---------------------------------------------------------------------------

def identify_channels(formula: str, notes: str = "") -> List[str]:
    """
    Identify potential pairing channels beyond phonon-mediated.

    Returns list of channel types: phonon, spin, topological, charge, orbital
    """
    f = formula.upper()
    ch = ["phonon"]  # always present in metals

    # Magnetic / spin-fluctuation channel
    magnetic_elements = {"FE", "CO", "NI", "MN", "CR", "CE", "PR",
                         "ND", "SM", "EU", "GD", "TB", "DY", "HO",
                         "ER", "TM", "YB", "U", "NP", "PU"}
    for el in magnetic_elements:
        if el in f:
            ch.append("spin")
            break

    # Topological channel (heavy elements with strong SOC)
    topo_elements = {"BI", "SB", "TE", "SE", "SN", "PB", "HG", "TL",
                     "IR", "PT", "AU", "W", "TA", "RE", "OS"}
    for el in topo_elements:
        if el in f:
            ch.append("topological")
            break

    # Charge-density-wave competing order
    cdw_elements = {"NB", "TA", "TI", "MO", "W"}
    has_chalcogen = any(el in f for el in ["S", "SE", "TE"])
    for el in cdw_elements:
        if el in f and has_chalcogen:
            ch.append("charge")
            break

    # Orbital fluctuations (multi-orbital systems)
    multiorbital = {"FE", "RU", "OS", "IR", "NI"}
    for el in multiorbital:
        if el in f:
            if "orbital" not in ch:
                ch.append("orbital")
            break

    return ch


# ---------------------------------------------------------------------------
# Built-in dataset (Mode 2)
# ---------------------------------------------------------------------------

def get_builtin_dataset() -> List[Material]:
    """
    Returns a curated dataset of ~40 known and predicted superconductors.
    """
    materials = []

    def add(mid, formula, sg, theta_D, lam, lam_src, Tc_exp,
            pressure=0.0, notes="", channels=None):
        m = Material(
            material_id=mid,
            formula=formula,
            space_group=sg,
            theta_D=theta_D,
            lambda_ep=lam,
            lambda_source=lam_src,
            Tc_experimental=Tc_exp,
            pressure_GPa=pressure,
            notes=notes,
        )
        if channels:
            m.channels = channels
        else:
            m.channels = identify_channels(formula)
        materials.append(m)

    # === High-pressure hydrides ===
    add("hyd-001", "H3S",    "Im-3m",   1330, 2.20, "exp",   203,   155,
        "Drozdov+2015, first high-Tc hydride")
    add("hyd-002", "LaH10",  "Fm-3m",   1200, 2.60, "calc",  260,   170,
        "Somayazulu+2019, near room-T")
    add("hyd-003", "YH6",    "Im-3m",   1020, 2.00, "calc",  224,   166,
        "Troyan+2021")
    add("hyd-004", "YH9",    "P63/mmc", 1100, 2.50, "calc",  243,   201,
        "Kong+2021")
    add("hyd-005", "CaH6",   "Im-3m",   1100, 2.20, "calc",  215,   172,
        "Ma+2022")
    add("hyd-006", "ThH10",  "Fm-3m",   1050, 2.40, "calc",  161,   175,
        "Semenok+2020")
    add("hyd-007", "PrH9",   "F-43m",    900, 1.80, "calc",    9,   120,
        "Zhou+2020, requires verification")
    add("hyd-008", "CeH9",   "P63/mmc",  850, 1.70, "calc",  100,   100,
        "Salke+2019")
    add("hyd-009", "BaH12",  "Cmc21",   1050, 2.10, "calc",   20,   150,
        "predicted, not yet observed at stated Tc")
    add("hyd-010", "MgH6",   "Im-3m",   1150, 1.60, "calc",  271,   300,
        "Feng+2015, predicted")
    add("hyd-011", "SnH12",  "I4/mmm",  1000, 1.90, "calc",   70,   250,
        "predicted")
    add("hyd-012", "AcH16",  "P-1",     1100, 2.80, "calc",   80,   150,
        "predicted by Semenok+2021")

    # === Lower-pressure or ambient superconductors ===
    add("amb-001", "MgB2",   "P6/mmm",   750, 0.87, "exp",    39,     0,
        "Nagamatsu+2001, ambient pressure")
    add("amb-002", "Nb",     "Im-3m",    275, 1.04, "exp",   9.3,     0,
        "elemental Nb")
    add("amb-003", "Nb3Sn",  "Pm-3n",    290, 1.80, "exp",    18,     0,
        "A15 phase")
    add("amb-004", "Nb3Ge",  "Pm-3n",    300, 1.70, "exp",    23,     0,
        "A15 phase")
    add("amb-005", "Pb",     "Fm-3m",    105, 1.55, "exp",   7.2,     0,
        "elemental Pb")
    add("amb-006", "V3Si",   "Pm-3n",    340, 1.30, "exp",    17,     0,
        "A15 phase")
    add("amb-007", "K3C60",  "Fm-3m",    200, 0.90, "exp",    19,     0,
        "fulleride, molecular SC")
    add("amb-008", "Cs3C60", "A15",      150, 1.20, "exp",    38,    0.7,
        "fulleride under pressure")

    # === Cuprates (unconventional, d-wave) ===
    add("cup-001", "YBa2Cu3O7",   "Pmmm",    410, 1.50, "eff",  92,    0,
        "YBCO, d-wave, non-BCS", ["phonon", "spin"])
    add("cup-002", "Bi2Sr2Ca2Cu3O10", "I4/mmm", 350, 1.60, "eff", 110, 0,
        "BSCCO-2223", ["phonon", "spin"])
    add("cup-003", "HgBa2Ca2Cu3O8", "P4/mmm",  380, 1.70, "eff", 133,  0,
        "Highest-Tc cuprate at ambient", ["phonon", "spin"])
    add("cup-004", "La2CuO4",  "I4/mmm",   360, 1.00, "eff",   38,     0,
        "Parent cuprate, doped", ["phonon", "spin"])

    # === Nickelates ===
    add("nic-001", "La3Ni2O7",  "Amam",   450, 1.20, "eff",  80,    14,
        "Sun+2023, bilayer nickelate", ["phonon", "spin", "orbital"])
    add("nic-002", "Nd6Ni5O12", "I4/mmm", 400, 1.00, "eff",  15,     0,
        "infinite-layer nickelate", ["phonon", "spin", "orbital"])

    # === Iron-based ===
    add("fe-001",  "LaFeAsOF",    "P4/nmm", 320, 1.10, "eff",  26,    0,
        "LaOFeAs-type, 1111 family", ["phonon", "spin", "orbital"])
    add("fe-002",  "FeSe",        "P4/nmm", 230, 0.80, "eff",   8,    0,
        "FeSe monolayer: ~65K on STO", ["phonon", "spin", "orbital"])

    # === tau-framework predicted candidates ===
    add("tau-001", "Mg2IrH6",  "Fm-3m",  750, 1.40, "tau-pred", 0,    50,
        "Predicted: Ir SOC + H phonon + Mg light", ["phonon", "topological"])
    add("tau-002", "SrAuH3",   "Pm-3m",  600, 1.10, "tau-pred", 0,    30,
        "Predicted: Au SOC + H phonon", ["phonon", "topological"])
    add("tau-003", "KB3C3",    "P6/mmm", 800, 1.00, "tau-pred", 0,     0,
        "Predicted: boron-carbon cage, multi-band", ["phonon", "charge"])
    add("tau-004", "La3Ni2O7H", "Amam",  500, 1.40, "tau-pred", 0,    14,
        "Proposed: H-doped nickelate", ["phonon", "spin", "orbital"])
    add("tau-005", "ScH12",    "Im-3m", 1200, 2.30, "tau-pred", 0,   200,
        "Predicted: light rare-earth hydride")
    add("tau-006", "LiB3H8",  "P-1",    900, 1.50, "tau-pred", 0,    50,
        "Predicted: lightweight borohydride")
    add("tau-007", "CaBH5",   "Cmcm",   950, 1.70, "tau-pred", 0,   100,
        "Predicted: Ca borohydride")
    add("tau-008", "BeH2",    "Ibam",   1200, 0.80, "calc",    0,    60,
        "Predicted: lightest hydride, high theta_D")
    add("tau-009", "LiH6",    "R-3m",   1100, 1.80, "calc",    0,   150,
        "Predicted: Li polyhydride")
    add("tau-010", "NaH9",    "Immm",   1050, 1.60, "calc",    0,   200,
        "Predicted: Na polyhydride")

    return materials


# ---------------------------------------------------------------------------
# Materials Project API (Mode 1)
# ---------------------------------------------------------------------------

def fetch_from_materials_project(api_key: str) -> List[Material]:
    """
    Query Materials Project for hydrogen-containing metallic/near-metallic
    stable compounds with Debye temperature data.
    """
    try:
        from mp_api.client import MPRester
    except ImportError:
        print("WARNING: mp-api not installed. Install with: pip install mp-api")
        print("Falling back to built-in dataset.\n")
        return get_builtin_dataset()

    materials = []

    print("Connecting to Materials Project API...")
    try:
        with MPRester(api_key) as mpr:
            # Search for H-containing, metallic or near-metallic, stable compounds
            docs = mpr.materials.summary.search(
                elements=["H"],
                band_gap=(0, 0.5),
                energy_above_hull=(0, 0.1),
                fields=[
                    "material_id", "formula_pretty", "symmetry",
                    "band_gap", "energy_above_hull",
                    "density", "volume",
                ],
                num_chunks=5,
            )

            print(f"Found {len(docs)} candidate materials from MP.")

            for doc in docs:
                mid = str(doc.material_id)
                formula = doc.formula_pretty
                sg = doc.symmetry.symbol if doc.symmetry else "unknown"
                bg = doc.band_gap if doc.band_gap is not None else 0.0
                eah = doc.energy_above_hull if doc.energy_above_hull is not None else 0.0

                # Estimate Debye temperature from density heuristic if not available
                # theta_D ~ (density)^{1/3} * C, rough scaling
                density = doc.density if doc.density else 3.0
                theta_D_est = 200.0 * (density / 3.0) ** (1.0 / 3.0)

                lam, lam_src = estimate_lambda(formula, theta_D_est)
                ch = identify_channels(formula)

                m = Material(
                    material_id=mid,
                    formula=formula,
                    space_group=sg,
                    theta_D=theta_D_est,
                    lambda_ep=lam,
                    lambda_source=lam_src,
                    band_gap=bg,
                    e_above_hull=eah,
                    density=density,
                    channels=ch,
                )
                materials.append(m)

    except Exception as e:
        print(f"API error: {e}")
        print("Falling back to built-in dataset.\n")
        return get_builtin_dataset()

    if not materials:
        print("No materials found via API. Using built-in dataset.\n")
        return get_builtin_dataset()

    return materials


# ---------------------------------------------------------------------------
# Custom input (Mode 3)
# ---------------------------------------------------------------------------

def custom_input_mode() -> List[Material]:
    """Interactive mode: user provides material parameters."""
    materials = []
    print("\n=== CUSTOM MATERIAL INPUT ===")
    print("Enter material parameters (type 'done' to finish):\n")

    idx = 0
    while True:
        idx += 1
        formula = input(f"Material {idx} formula (or 'done'): ").strip()
        if formula.lower() == 'done':
            break

        try:
            theta_D = float(input("  Debye temperature Theta_D (K): ").strip())
        except ValueError:
            print("  Invalid. Skipping.")
            continue

        lam_input = input("  Electron-phonon coupling lambda (or 'auto'): ").strip()
        if lam_input.lower() == 'auto':
            lam, lam_src = estimate_lambda(formula, theta_D)
            print(f"  -> Auto-estimated lambda = {lam:.2f} ({lam_src})")
        else:
            try:
                lam = float(lam_input)
                lam_src = "user-input"
            except ValueError:
                lam, lam_src = estimate_lambda(formula, theta_D)
                print(f"  -> Invalid input, auto-estimated lambda = {lam:.2f}")

        nef_input = input("  N(E_F) in states/eV/atom (or 'skip'): ").strip()
        notes = ""
        if nef_input.lower() != 'skip':
            try:
                nef = float(nef_input)
                notes = f"N(E_F)={nef}"
            except ValueError:
                pass

        channels_input = input("  Pairing channels (comma-separated, or 'auto'): ").strip()
        if channels_input.lower() == 'auto':
            ch = identify_channels(formula)
        else:
            ch = [c.strip() for c in channels_input.split(",")]

        pressure = 0.0
        p_input = input("  Pressure (GPa, or 0 for ambient): ").strip()
        try:
            pressure = float(p_input)
        except ValueError:
            pass

        m = Material(
            material_id=f"custom-{idx:03d}",
            formula=formula,
            theta_D=theta_D,
            lambda_ep=lam,
            lambda_source=lam_src,
            pressure_GPa=pressure,
            channels=ch,
            notes=notes,
        )
        materials.append(m)
        print()

    return materials


# ---------------------------------------------------------------------------
# Analysis pipeline
# ---------------------------------------------------------------------------

def analyze_materials(materials: List[Material], T_ref: float = 300.0,
                      mu_star: float = 0.13) -> List[Material]:
    """
    For each material:
      1. Estimate T_c via McMillan (unless experimental T_c is known)
      2. Compute R_tau at reference temperature
      3. Compute multi-channel R_tau
      4. Compute channels needed for R_tau = 1 at 300 K
    """
    for m in materials:
        # Estimate T_c
        Tc_mcm = mcmillan_Tc(m.theta_D, m.lambda_ep, mu_star)
        Tc_ad = allen_dynes_Tc(m.theta_D, m.lambda_ep, mu_star)
        m.Tc_estimated = max(Tc_mcm, Tc_ad)

        # Use experimental if available and higher (for unconventional SC)
        Tc_use = max(m.Tc_experimental, m.Tc_estimated)
        if Tc_use == 0:
            Tc_use = m.Tc_estimated

        # Compute R_tau
        m.R_tau = compute_R_tau(Tc_use, m.theta_D, T_ref)

        # Multi-channel R_tau
        alpha = len(m.channels)
        m.R_tau_multi = alpha * m.R_tau

        # Channels needed for room-T SC
        if m.R_tau > 0:
            m.channels_needed_300K = channels_needed(
                m.R_tau, target_Tc=300.0, T=T_ref, theta_D=m.theta_D
            )
        else:
            m.channels_needed_300K = 999

    # Sort by multi-channel R_tau, descending
    materials.sort(key=lambda m: m.R_tau_multi, reverse=True)
    return materials


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def print_report(materials: List[Material], T_ref: float = 300.0,
                 top_n: int = 20, mu_star: float = 0.13):
    """Print the screening report."""

    temps_C = [0, 20, 37, 90]
    temps_K = [t + 273.15 for t in temps_C]

    sep = "=" * 100
    print(f"\n{sep}")
    print("  R_tau SUPERCONDUCTOR SCREENING REPORT")
    print(f"  Reference temperature: {T_ref:.1f} K")
    print(f"  Coulomb pseudopotential mu* = {mu_star}")
    print(f"  Total materials screened: {len(materials)}")
    print(f"{sep}\n")

    # ---- Main table ----
    print(f"Top {min(top_n, len(materials))} candidates (sorted by R_tau_multi):\n")

    header = (f"{'Rank':>4}  {'Material':<16} {'SG':<10} {'theta_D':>7} "
              f"{'lambda':>6} {'lam_src':<12} {'Tc_est':>7} {'Tc_exp':>7} "
              f"{'R_tau':>8} {'R_multi':>8} {'#Ch':>3} {'Ch_need':>7} "
              f"{'P(GPa)':>7}")
    print(header)
    print("-" * len(header))

    for i, m in enumerate(materials[:top_n]):
        Tc_use = max(m.Tc_experimental, m.Tc_estimated)
        print(f"{i+1:>4}  {m.formula:<16} {m.space_group:<10} {m.theta_D:>7.0f} "
              f"{m.lambda_ep:>6.2f} {m.lambda_source:<12} {m.Tc_estimated:>7.1f} "
              f"{m.Tc_experimental:>7.1f} {m.R_tau:>8.4f} {m.R_tau_multi:>8.4f} "
              f"{len(m.channels):>3} {m.channels_needed_300K:>7} "
              f"{m.pressure_GPa:>7.0f}")

    # ---- Temperature dependence ----
    print(f"\n\n{'=' * 80}")
    print("  R_tau at different temperatures (top 10)")
    print(f"{'=' * 80}\n")

    header2 = (f"{'Material':<16} " +
               " ".join(f"{'R@' + str(t) + 'C':>10}" for t in temps_C) +
               "  " +
               " ".join(f"{'lam@' + str(t) + 'C':>10}" for t in temps_C))
    print(header2)
    print("-" * len(header2))

    for m in materials[:10]:
        Tc_use = max(m.Tc_experimental, m.Tc_estimated)
        r_vals = [compute_R_tau(Tc_use, m.theta_D, T) for T in temps_K]
        lam_vals = [lambda_needed_for_Rtau1(m.theta_D, T, mu_star) for T in temps_K]

        r_str = " ".join(f"{r:>10.4f}" for r in r_vals)
        l_str = " ".join(f"{l:>10.2f}" if l < 100 else f"{'inf':>10}" for l in lam_vals)
        print(f"{m.formula:<16} {r_str}  {l_str}")

    # ---- Multi-channel analysis ----
    print(f"\n\n{'=' * 80}")
    print("  MULTI-CHANNEL ANALYSIS")
    print(f"{'=' * 80}\n")

    print(f"{'Material':<16} {'Channels':<35} {'alpha':>5} {'R_single':>10} "
          f"{'R_multi':>10} {'alpha_300K':>10}")
    print("-" * 100)

    for m in materials[:top_n]:
        ch_str = "+".join(m.channels)
        alpha = len(m.channels)
        alpha_300 = m.channels_needed_300K
        print(f"{m.formula:<16} {ch_str:<35} {alpha:>5} {m.R_tau:>10.4f} "
              f"{m.R_tau_multi:>10.4f} {alpha_300:>10}")

    # ---- Room-T feasibility ----
    print(f"\n\n{'=' * 80}")
    print("  ROOM-TEMPERATURE FEASIBILITY ANALYSIS (T = 300 K)")
    print(f"{'=' * 80}\n")

    print("For R_tau >= 1 at 300 K, we need Tc * theta_D >= 90000 K^2\n")

    # Find materials closest to R_tau = 1
    print("Materials closest to R_tau = 1 (single channel):\n")
    for m in materials[:10]:
        Tc_use = max(m.Tc_experimental, m.Tc_estimated)
        deficit = 1.0 - m.R_tau
        Tc_needed = 300.0 * 300.0 / m.theta_D if m.theta_D > 0 else float('inf')
        lambda_for_Tc = lambda_needed_for_Rtau1(m.theta_D, 300.0, mu_star)

        print(f"  {m.formula:<16}  R_tau={m.R_tau:.4f}  "
              f"deficit={deficit:+.4f}  "
              f"Tc_needed={Tc_needed:.0f}K  "
              f"lambda_needed={lambda_for_Tc:.2f}")

    # ---- Best multi-channel route ----
    print(f"\n\nBest route to room-temperature SC via multi-channel enhancement:\n")

    # Sort by channels_needed
    by_channels = sorted(materials, key=lambda m: m.channels_needed_300K)
    for m in by_channels[:10]:
        Tc_use = max(m.Tc_experimental, m.Tc_estimated)
        print(f"  {m.formula:<16}  alpha_current={len(m.channels):>2}  "
              f"alpha_needed={m.channels_needed_300K:>3}  "
              f"R_tau_single={m.R_tau:.4f}  "
              f"Tc={'%.0f' % Tc_use:>5}K  "
              f"theta_D={m.theta_D:.0f}K  "
              f"P={m.pressure_GPa:.0f} GPa")

    # ---- Key insights ----
    print(f"\n\n{'=' * 80}")
    print("  KEY INSIGHTS FROM R_tau FRAMEWORK")
    print(f"{'=' * 80}\n")

    # Find best single-channel
    best_single = max(materials, key=lambda m: m.R_tau)
    best_multi = max(materials, key=lambda m: m.R_tau_multi)
    best_ambient = max(
        [m for m in materials if m.pressure_GPa < 1],
        key=lambda m: m.R_tau_multi,
        default=materials[0]
    )

    print(f"  Highest R_tau (single channel):  {best_single.formula} "
          f"(R_tau = {best_single.R_tau:.4f})")
    print(f"  Highest R_tau (multi-channel):   {best_multi.formula} "
          f"(R_tau_multi = {best_multi.R_tau_multi:.4f})")
    print(f"  Best ambient-pressure candidate: {best_ambient.formula} "
          f"(R_tau_multi = {best_ambient.R_tau_multi:.4f})")

    # R_tau = 1 threshold
    above_1 = [m for m in materials if m.R_tau >= 1.0]
    above_1_multi = [m for m in materials if m.R_tau_multi >= 1.0]
    print(f"\n  Materials with R_tau >= 1 (single): {len(above_1)}")
    for m in above_1:
        print(f"    - {m.formula} (R_tau = {m.R_tau:.4f}, P = {m.pressure_GPa:.0f} GPa)")
    print(f"  Materials with R_tau_multi >= 1:     {len(above_1_multi)}")
    for m in above_1_multi:
        if m not in above_1:
            print(f"    - {m.formula} (R_tau_multi = {m.R_tau_multi:.4f})")

    print(f"\n  The R_tau framework suggests that room-temperature SC requires:")
    print(f"    - theta_D > 800 K (light elements, strong bonds)")
    print(f"    - lambda > 2.0 (strong electron-phonon coupling)")
    print(f"    - Multiple pairing channels (alpha >= 2-3)")
    print(f"    - Or: theta_D > 1200 K with lambda > 2.5 for single-channel\n")

    return


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def save_csv(materials: List[Material], filename: str):
    """Save results to CSV file."""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "rank", "material_id", "formula", "space_group",
            "theta_D_K", "lambda", "lambda_source",
            "Tc_estimated_K", "Tc_experimental_K",
            "R_tau", "R_tau_multi", "n_channels", "channels",
            "channels_needed_300K", "pressure_GPa", "notes"
        ])
        for i, m in enumerate(materials):
            writer.writerow([
                i + 1, m.material_id, m.formula, m.space_group,
                f"{m.theta_D:.1f}", f"{m.lambda_ep:.3f}", m.lambda_source,
                f"{m.Tc_estimated:.2f}", f"{m.Tc_experimental:.1f}",
                f"{m.R_tau:.6f}", f"{m.R_tau_multi:.6f}",
                len(m.channels), "+".join(m.channels),
                m.channels_needed_300K, f"{m.pressure_GPa:.1f}", m.notes
            ])
    print(f"\nResults saved to: {filename}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="R_tau Superconductor Screening Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rtau_screener.py                        # built-in dataset
  python rtau_screener.py --api-key mp-XXXXX     # Materials Project API
  python rtau_screener.py --custom               # interactive input
  python rtau_screener.py --output results.csv   # save CSV
  python rtau_screener.py --temperature 77       # screen at 77 K (LN2)
  python rtau_screener.py --mu-star 0.10         # different mu*
        """
    )
    parser.add_argument("--api-key", type=str, default=None,
                        help="Materials Project API key")
    parser.add_argument("--custom", action="store_true",
                        help="Interactive custom material input mode")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output CSV filename")
    parser.add_argument("--temperature", "-T", type=float, default=300.0,
                        help="Reference temperature in K (default: 300)")
    parser.add_argument("--mu-star", type=float, default=0.13,
                        help="Coulomb pseudopotential mu* (default: 0.13)")
    parser.add_argument("--top", type=int, default=20,
                        help="Number of top candidates to display (default: 20)")
    parser.add_argument("--all", action="store_true",
                        help="Show all materials, not just top N")

    args = parser.parse_args()

    # Select data source
    if args.custom:
        materials = custom_input_mode()
        if not materials:
            print("No materials entered. Exiting.")
            return
    elif args.api_key:
        materials = fetch_from_materials_project(args.api_key)
    else:
        print("Using built-in dataset of known and predicted superconductors.")
        print("(Use --api-key YOUR_KEY to query Materials Project)\n")
        materials = get_builtin_dataset()

    # Run analysis
    materials = analyze_materials(materials, T_ref=args.temperature,
                                  mu_star=args.mu_star)

    # Print report
    top_n = len(materials) if args.all else args.top
    print_report(materials, T_ref=args.temperature, top_n=top_n,
                 mu_star=args.mu_star)

    # Save CSV
    if args.output:
        save_csv(materials, args.output)
    else:
        # Default CSV output
        default_csv = args.output or "rtau_screening_results.csv"
        save_csv(materials, default_csv)


if __name__ == "__main__":
    main()
