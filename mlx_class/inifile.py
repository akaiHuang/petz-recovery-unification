"""
inifile.py -- CLASS-compatible .ini parameter file parser for mlx_class.

Parses CLASS-style .ini files and maps CLASS parameter names to mlx_class
parameter names, then creates and runs a UnifiedSolver.

CLASS .ini format:
    h = 0.6736
    omega_b = 0.02237
    omega_cdm = 0.12
    output = tCl,pCl,lCl,mPk
    l_max_scalars = 2500

Lines starting with '#' or ';' are comments.
Empty lines are ignored.
Inline comments (after '#') are stripped.

Usage:
    from mlx_class.inifile import load_ini, run_from_ini

    params = load_ini('my_cosmology.ini')
    result = run_from_ini('my_cosmology.ini')

Author: Sheng-Kai Huang, 2026
"""
import os
import sys
from typing import Dict, Any, Optional, Tuple

# ---------------------------------------------------------------------------
# CLASS -> mlx_class parameter name mapping
# ---------------------------------------------------------------------------
# Keys: CLASS .ini name (case-sensitive, as CLASS uses)
# Values: mlx_class UnifiedSolver kwarg name
_CLASS_TO_MLX = {
    # Cosmological parameters
    'h':                'h',
    'omega_b':          'omega_b',
    'omega_cdm':        'omega_cdm',
    'T_cmb':            'T_CMB',
    'A_s':              'A_s',
    'n_s':              'n_s',
    'tau_reio':         'tau_reio',
    'N_ur':             'N_ur',
    'N_eff':            'N_ur',           # alias
    'N_ncdm':           'N_ncdm',
    'm_ncdm':           'sum_mnu',        # CLASS m_ncdm is per-species; we sum
    # Dark energy (CPL)
    'w0_fld':           'w0',
    'wa_fld':           'wa',
    'w0':               'w0',             # alternative name
    'wa':               'wa',             # alternative name
    # Curvature
    'Omega_k':          'Omega_k',
    # Tensors
    'r':                'tensor_to_scalar_ratio',
    # Precision / grid
    'l_max_scalars':    'ell_max',
    'l_max':            'ell_max',        # alias
    # Solver mode
    'modes':            'modes',
    'output':           'output',
}

# Parameters that should be parsed as float
_FLOAT_PARAMS = {
    'h', 'omega_b', 'omega_cdm', 'T_CMB', 'A_s', 'n_s', 'tau_reio',
    'N_ur', 'sum_mnu', 'w0', 'wa', 'Omega_k', 'tensor_to_scalar_ratio',
}

# Parameters that should be parsed as int
_INT_PARAMS = {
    'ell_max', 'N_k', 'N_ncdm',
}

# Parameters that remain as string
_STRING_PARAMS = {
    'output', 'modes',
}


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def load_ini(filepath: str) -> Dict[str, Any]:
    """
    Parse a CLASS-style .ini file and return mlx_class parameters.

    Parameters
    ----------
    filepath : str
        Path to the .ini file.

    Returns
    -------
    dict
        Dictionary of mlx_class parameter names -> values.
        Includes a special key '_output_flags' with parsed output requests:
        {'tCl': bool, 'pCl': bool, 'lCl': bool, 'mPk': bool, 'tCl_lensed': bool}

    Raises
    ------
    FileNotFoundError
        If the .ini file does not exist.
    ValueError
        If a line cannot be parsed.
    """
    filepath = os.path.expanduser(filepath)
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"INI file not found: {filepath}")

    raw_params = _parse_raw(filepath)
    mlx_params = _map_to_mlx(raw_params)
    return mlx_params


def _parse_raw(filepath: str) -> Dict[str, str]:
    """
    Read a CLASS .ini file and return raw key-value string pairs.

    Handles:
    - Comment lines (# or ;)
    - Inline comments (# after value)
    - Whitespace around = sign
    - Empty lines
    - Values with scientific notation (2.1e-9)
    - Comma-separated lists (output = tCl,pCl,lCl)
    """
    params = {}
    with open(filepath, 'r') as f:
        for line_no, line in enumerate(f, 1):
            # Strip trailing whitespace/newline
            line = line.strip()

            # Skip empty lines and comment lines
            if not line or line.startswith('#') or line.startswith(';'):
                continue

            # Remove inline comments (but be careful with # in values)
            # CLASS convention: # always starts a comment
            comment_idx = line.find('#')
            if comment_idx >= 0:
                line = line[:comment_idx].strip()

            if not line:
                continue

            # Split on first '='
            if '=' not in line:
                # CLASS also supports lines without '=' for some directives,
                # but for standard parameters they always have '='.
                # Skip gracefully.
                continue

            key, _, value = line.partition('=')
            key = key.strip()
            value = value.strip()

            if not key:
                continue

            params[key] = value

    return params


def _map_to_mlx(raw_params: Dict[str, str]) -> Dict[str, Any]:
    """
    Map CLASS raw parameter names/values to mlx_class typed parameters.
    """
    mlx_params = {}
    unmapped = {}

    # Parse output flags first (needed for solver configuration)
    output_flags = {
        'tCl': False,
        'pCl': False,
        'lCl': False,
        'mPk': False,
    }

    for class_key, raw_value in raw_params.items():
        # Look up the mlx_class name
        mlx_key = _CLASS_TO_MLX.get(class_key)

        if mlx_key is None:
            # Check case-insensitive fallback for common variations
            mlx_key = _CLASS_TO_MLX.get(class_key.lower())

        if mlx_key is None:
            # Store unmapped parameters for diagnostics (not an error --
            # CLASS has many parameters we don't support yet)
            unmapped[class_key] = raw_value
            continue

        # Special handling for 'output'
        if mlx_key == 'output':
            tokens = [t.strip() for t in raw_value.split(',')]
            for token in tokens:
                t_lower = token.lower()
                if t_lower in ('tcl', 'tcl'):
                    output_flags['tCl'] = True
                elif t_lower in ('pcl',):
                    output_flags['pCl'] = True
                elif t_lower in ('lcl',):
                    output_flags['lCl'] = True
                elif t_lower in ('mpk',):
                    output_flags['mPk'] = True
            continue

        # Special handling for 'modes' (e.g., modes = s,t)
        if mlx_key == 'modes':
            mlx_params['_modes'] = [m.strip() for m in raw_value.split(',')]
            continue

        # Special handling for m_ncdm (may be comma-separated list of masses)
        if class_key == 'm_ncdm':
            masses = [float(m.strip()) for m in raw_value.split(',') if m.strip()]
            mlx_params['sum_mnu'] = sum(masses)
            continue

        # Type conversion
        if mlx_key in _FLOAT_PARAMS:
            mlx_params[mlx_key] = float(raw_value)
        elif mlx_key in _INT_PARAMS:
            mlx_params[mlx_key] = int(float(raw_value))
        else:
            mlx_params[mlx_key] = raw_value

    # Store output flags
    mlx_params['_output_flags'] = output_flags

    # Store unmapped keys for diagnostics
    if unmapped:
        mlx_params['_unmapped'] = unmapped

    return mlx_params


# ---------------------------------------------------------------------------
# Solver creation and execution
# ---------------------------------------------------------------------------

def create_solver(params: Dict[str, Any]):
    """
    Create a UnifiedSolver from parsed .ini parameters.

    Parameters
    ----------
    params : dict
        Output of load_ini().

    Returns
    -------
    UnifiedSolver
        Configured solver ready to call .compute_all().
    """
    from mlx_class.solver_unified import UnifiedSolver

    # Extract the subset of params that UnifiedSolver accepts
    solver_kwargs = {}

    # Direct numerical parameters
    _direct_map = {
        'h':         'h',
        'omega_b':   'omega_b',
        'omega_cdm': 'omega_cdm',
        'A_s':       'A_s',
        'n_s':       'n_s',
        'tau_reio':  'tau_reio',
        'N_ur':      'N_ur',
        'sum_mnu':   'sum_mnu',
        'w0':        'w0',
        'wa':        'wa',
        'ell_max':   'l_max',
        'N_k':       'N_k',
    }

    for ini_key, solver_key in _direct_map.items():
        if ini_key in params:
            solver_kwargs[solver_key] = params[ini_key]

    # Output flags -> solver configuration
    output_flags = params.get('_output_flags', {})

    # lCl -> compute_lensing
    if output_flags.get('lCl', False):
        solver_kwargs['compute_lensing'] = True
    elif 'lCl' in output_flags:
        solver_kwargs['compute_lensing'] = output_flags['lCl']

    # mPk -> compute_pk
    if output_flags.get('mPk', False):
        solver_kwargs['compute_pk'] = True

    # Khronon mode (not a standard CLASS parameter, but we support it)
    if params.get('khronon', False):
        solver_kwargs['khronon'] = True

    # Fast mode
    if params.get('fast_mode', False):
        solver_kwargs['fast_mode'] = True

    # Print configuration summary
    print("=" * 72)
    print("  INI FILE CONFIGURATION")
    print("=" * 72)
    for k, v in sorted(solver_kwargs.items()):
        print(f"  {k:20s} = {v}")

    # Report unmapped CLASS parameters
    unmapped = params.get('_unmapped', {})
    if unmapped:
        print(f"\n  [INFO] {len(unmapped)} CLASS parameter(s) not mapped "
              f"(ignored):")
        for k, v in sorted(unmapped.items()):
            print(f"    {k} = {v}")
    print("=" * 72)

    solver = UnifiedSolver(**solver_kwargs)
    return solver


def run_from_ini(filepath: str, quiet: bool = False):
    """
    Load a CLASS .ini file and run the full pipeline.

    Parameters
    ----------
    filepath : str
        Path to the .ini file.
    quiet : bool
        If True, suppress most output.

    Returns
    -------
    UnifiedResult
        The computation result with Dl_TT, Dl_TE, Dl_EE, etc.

    Examples
    --------
    >>> from mlx_class.inifile import run_from_ini
    >>> result = run_from_ini('planck_bestfit.ini')
    >>> print(result.Dl_TT[:5])
    """
    params = load_ini(filepath)

    if not quiet:
        print(f"\n[inifile] Loaded: {os.path.abspath(filepath)}")
        output_flags = params.get('_output_flags', {})
        active = [k for k, v in output_flags.items() if v]
        if active:
            print(f"[inifile] Output: {', '.join(active)}")

    solver = create_solver(params)
    result = solver.compute_all()
    return result


def params_to_solver_kwargs(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert parsed .ini params to UnifiedSolver keyword arguments.

    This is useful when you want to inspect or modify parameters
    before creating the solver.

    Parameters
    ----------
    params : dict
        Output of load_ini().

    Returns
    -------
    dict
        Keyword arguments suitable for UnifiedSolver(**kwargs).
    """
    _direct_map = {
        'h':         'h',
        'omega_b':   'omega_b',
        'omega_cdm': 'omega_cdm',
        'A_s':       'A_s',
        'n_s':       'n_s',
        'tau_reio':  'tau_reio',
        'N_ur':      'N_ur',
        'sum_mnu':   'sum_mnu',
        'w0':        'w0',
        'wa':        'wa',
        'ell_max':   'l_max',
        'N_k':       'N_k',
    }

    kwargs = {}
    for ini_key, solver_key in _direct_map.items():
        if ini_key in params:
            kwargs[solver_key] = params[ini_key]

    output_flags = params.get('_output_flags', {})
    if output_flags.get('lCl', False):
        kwargs['compute_lensing'] = True
    if output_flags.get('mPk', False):
        kwargs['compute_pk'] = True

    return kwargs


# ---------------------------------------------------------------------------
# Comparison / verification utility
# ---------------------------------------------------------------------------

def verify_ini_vs_api(filepath: str) -> Tuple[bool, str]:
    """
    Load a .ini file and compare results with the equivalent Python API call.

    This runs the solver twice -- once from .ini and once from direct API --
    and checks that the TT spectrum matches to within numerical precision.

    Parameters
    ----------
    filepath : str
        Path to the .ini file.

    Returns
    -------
    (passed, report) : (bool, str)
        Whether the comparison passed, and a diagnostic report.
    """
    import numpy as np
    from mlx_class.solver_unified import UnifiedSolver

    # Run from .ini
    print("[verify] Running from .ini file...")
    params = load_ini(filepath)
    solver_ini = create_solver(params)
    result_ini = solver_ini.compute_all()

    # Run from Python API with equivalent parameters
    print("[verify] Running from Python API...")
    kwargs = params_to_solver_kwargs(params)
    solver_api = UnifiedSolver(**kwargs)
    result_api = solver_api.compute_all()

    # Compare TT spectra
    lines = []
    lines.append(f"INI file: {os.path.abspath(filepath)}")
    lines.append(f"INI  ell range: [{result_ini.ell[0]}, {result_ini.ell[-1]}], "
                 f"N_ell = {len(result_ini.ell)}")
    lines.append(f"API  ell range: [{result_api.ell[0]}, {result_api.ell[-1]}], "
                 f"N_ell = {len(result_api.ell)}")

    # Find common ell range
    ell_common = np.intersect1d(result_ini.ell, result_api.ell)
    if len(ell_common) == 0:
        return False, "No common ell values found."

    # Extract TT at common ells
    ini_idx = np.searchsorted(result_ini.ell, ell_common)
    api_idx = np.searchsorted(result_api.ell, ell_common)

    Dl_ini = result_ini.Dl_TT[ini_idx]
    Dl_api = result_api.Dl_TT[api_idx]

    # Compute relative difference
    mask = np.abs(Dl_api) > 1.0  # avoid division by near-zero
    if np.any(mask):
        rel_diff = np.abs(Dl_ini[mask] - Dl_api[mask]) / np.abs(Dl_api[mask])
        max_rel = np.max(rel_diff)
        mean_rel = np.mean(rel_diff)
    else:
        max_rel = 0.0
        mean_rel = 0.0

    lines.append(f"Common ells: {len(ell_common)}")
    lines.append(f"Max  relative difference (Dl_TT): {max_rel:.2e}")
    lines.append(f"Mean relative difference (Dl_TT): {mean_rel:.2e}")

    # Pass if max difference < 1e-10 (should be identical -- same code path)
    passed = max_rel < 1e-10
    lines.append(f"PASSED: {passed}")

    report = '\n'.join(lines)
    return passed, report


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def cli_main():
    """
    Command-line interface for running mlx_class from a .ini file.

    Usage:
        python -m mlx_class my_cosmology.ini [--verify] [--no-plot] [--quiet]
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog='mlx_class',
        description='mlx_class: GPU-accelerated CMB Boltzmann solver. '
                    'Run from a CLASS-compatible .ini file.',
        epilog='Example: python -m mlx_class planck_bestfit.ini',
    )
    parser.add_argument(
        'inifile',
        nargs='?',
        default=None,
        help='Path to a CLASS-compatible .ini parameter file.',
    )
    parser.add_argument(
        '--verify',
        action='store_true',
        help='Verify that .ini results match Python API results.',
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress most output.',
    )
    parser.add_argument(
        '--dump-params',
        action='store_true',
        help='Parse the .ini file and print parameters without running.',
    )

    args = parser.parse_args()

    if args.inifile is None:
        parser.print_help()
        sys.exit(0)

    filepath = args.inifile

    # Resolve relative paths: check current directory, then mlx_class directory
    if not os.path.isabs(filepath) and not os.path.isfile(filepath):
        # Try looking in the mlx_class package directory
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        pkg_path = os.path.join(pkg_dir, filepath)
        if os.path.isfile(pkg_path):
            filepath = pkg_path

    if args.dump_params:
        params = load_ini(filepath)
        print(f"\nParsed parameters from: {os.path.abspath(filepath)}")
        print("-" * 50)
        for k, v in sorted(params.items()):
            if k.startswith('_'):
                continue
            print(f"  {k:25s} = {v}")
        print("-" * 50)
        flags = params.get('_output_flags', {})
        print(f"  Output flags: {flags}")
        unmapped = params.get('_unmapped', {})
        if unmapped:
            print(f"  Unmapped: {unmapped}")
        sys.exit(0)

    if args.verify:
        passed, report = verify_ini_vs_api(filepath)
        print("\n" + "=" * 60)
        print("  VERIFICATION REPORT")
        print("=" * 60)
        print(report)
        print("=" * 60)
        sys.exit(0 if passed else 1)

    # Normal run
    result = run_from_ini(filepath, quiet=args.quiet)

    # Print summary
    if not args.quiet and result.peak_ells is not None and len(result.peak_ells) > 0:
        print(f"\nPeak positions: {result.peak_ells[:5]}")
        if result.timing:
            print(f"Total time: {result.timing.get('total', 0):.2f}s")


if __name__ == '__main__':
    cli_main()
