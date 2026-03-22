"""
__main__.py -- Entry point for `python -m mlx_class`.

Supports three modes:
  1. Unified API:  python -m mlx_class [--backend sync|fast] [--N_k 180] [--quiet]
  2. INI file mode: python -m mlx_class my_cosmology.ini [--verify] [--quiet]
  3. Legacy mode:   python -m mlx_class --ode [--implicit] [--khronon] ...

If the first positional argument is a .ini file, delegates to inifile.cli_main().
If --backend is specified, uses the unified CMBSolver API.
Otherwise, falls back to the legacy main.main() entry point.

Author: Sheng-Kai Huang, 2026
"""
import sys
import os


def _is_ini_invocation():
    """Check if the user is passing a .ini file as the first argument."""
    for arg in sys.argv[1:]:
        if arg.startswith('-'):
            # Skip flags -- but check for --dump-params and --verify
            # which are ini-mode flags
            if arg in ('--dump-params', '--verify'):
                return True
            continue
        # First positional argument: check if it looks like an .ini file
        if arg.endswith('.ini'):
            return True
        # Also check if it's a file that exists (even without .ini extension)
        if os.path.isfile(arg):
            return True
        # Check in the package directory
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        if os.path.isfile(os.path.join(pkg_dir, arg)):
            return True
        # Not a file -- probably a legacy flag value
        return False
    return False


def _is_unified_invocation():
    """Check if the user passed --backend, indicating unified API mode."""
    return '--backend' in sys.argv[1:]


def _run_unified():
    """Run via the unified CMBSolver API."""
    import argparse
    from .solver_api import CMBSolver

    parser = argparse.ArgumentParser(description='mlx_class CMB solver')
    parser.add_argument('--backend', default='sync', choices=['sync', 'fast'])
    parser.add_argument('--N_k', type=int, default=None)
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()

    solver = CMBSolver(backend=args.backend, N_k=args.N_k,
                       verbose=not args.quiet)
    result = solver.run()

    print(f"\nPeaks: {result.peaks['ell'][:5]}")
    print(f"Backend: {result.backend}")


if _is_ini_invocation():
    from mlx_class.inifile import cli_main
    cli_main()
elif _is_unified_invocation():
    _run_unified()
else:
    from mlx_class.main import main
    main()
