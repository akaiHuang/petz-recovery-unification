"""
__main__.py -- Entry point for `python -m mlx_class`.

Supports two modes:
  1. INI file mode:  python -m mlx_class my_cosmology.ini [--verify] [--quiet]
  2. Legacy mode:    python -m mlx_class --ode [--implicit] [--khronon] ...

If the first positional argument is a .ini file, delegates to inifile.cli_main().
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


if _is_ini_invocation():
    from mlx_class.inifile import cli_main
    cli_main()
else:
    from mlx_class.main import main
    main()
