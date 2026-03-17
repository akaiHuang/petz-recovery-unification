"""tau-chrono command-line interface."""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    """Entry point for the tau-chrono CLI."""
    parser = argparse.ArgumentParser(
        prog="tau-chrono",
        description="Bayesian Noise Tracker for Quantum Circuits",
    )
    parser.add_argument(
        "--version", action="store_true", help="Show version and exit"
    )
    sub = parser.add_subparsers(dest="command")

    # --- benchmark ---
    bench_parser = sub.add_parser(
        "benchmark", help="Run built-in benchmarks"
    )
    bench_parser.add_argument(
        "--suite",
        choices=["1q", "2q", "all"],
        default="1q",
        help="Benchmark suite to run (default: 1q)",
    )

    # --- characterize ---
    char_parser = sub.add_parser(
        "characterize", help="Characterize gate noise on a backend"
    )
    char_parser.add_argument(
        "--backend", type=str, default="simulator",
        help="Backend name (default: simulator)",
    )
    char_parser.add_argument(
        "--gates", type=str, default="h,x,sx,rz,s",
        help="Comma-separated gate list",
    )
    char_parser.add_argument(
        "--shots", type=int, default=4096,
        help="Shots per circuit (default: 4096)",
    )

    # --- analyze ---
    analyze_parser = sub.add_parser(
        "analyze", help="Analyze a gate sequence"
    )
    analyze_parser.add_argument(
        "gates", type=str,
        help="Comma-separated gate/noise sequence, e.g. 'ad:0.05,dep:0.08,ad:0.12'",
    )

    args = parser.parse_args(argv)

    if args.version:
        from tau_chrono import __version__
        print(f"tau-chrono {__version__}")
        return 0

    if args.command == "benchmark":
        return _run_benchmark(args.suite)
    elif args.command == "characterize":
        return _run_characterize(args.backend, args.gates, args.shots)
    elif args.command == "analyze":
        return _run_analyze(args.gates)
    else:
        parser.print_help()
        return 0


def _run_benchmark(suite: str) -> int:
    """Run built-in benchmarks."""
    if suite in ("1q", "all"):
        print("Running single-qubit benchmark...")
        from tau_chrono.benchmarks import single_qubit
        single_qubit.main()
    if suite in ("2q", "all"):
        print("\nRunning two-qubit benchmark...")
        from tau_chrono.benchmarks import two_qubit
        two_qubit.main()
    return 0


def _run_characterize(backend_name: str, gates_str: str, shots: int) -> int:
    """Characterize gates on a backend."""
    print(f"Characterizing gates [{gates_str}] on {backend_name} ({shots} shots)")
    print("(Hardware characterization requires qiskit. Use: pip install tau-chrono[qiskit])")
    # TODO: implement when backend adapters are ready
    return 0


def _run_analyze(gates_str: str) -> int:
    """Analyze a gate/noise sequence."""
    import numpy as np
    from tau_chrono import (
        amplitude_damping, depolarizing, dephasing,
        bayesian_compose,
    )

    channel_map = {
        "ad": amplitude_damping,
        "dep": depolarizing,
        "deph": dephasing,
    }

    channels = []
    names = []
    for spec in gates_str.split(","):
        spec = spec.strip()
        parts = spec.split(":")
        if len(parts) != 2:
            print(f"Error: invalid gate spec '{spec}'. Use format 'type:param', e.g. 'ad:0.05'")
            return 1
        ch_type, param = parts[0], float(parts[1])
        if ch_type not in channel_map:
            print(f"Error: unknown channel type '{ch_type}'. Available: {list(channel_map.keys())}")
            return 1
        channels.append(channel_map[ch_type](param))
        names.append(f"{ch_type}({param})")

    rho = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
    sigma = np.diag([0.8, 0.2]).astype(complex)

    result = bayesian_compose(channels, sigma, rho, channel_names=names)

    print(f"\n{'Gate':<15} {'tau_naive':>10} {'tau_eff':>10} {'Classification'}")
    print("-" * 50)
    for gr in result.gate_results:
        print(f"{gr.channel_name:<15} {gr.tau_naive:>10.6f} {gr.tau_eff:>10.6f} {gr.classification}")

    print(f"\n{'Total tau (Bayesian):':<30} {result.tau_bayesian_total:.6f}")
    print(f"{'Total tau (multiplicative):':<30} {result.tau_multiplicative_total:.6f}")
    print(f"{'Improvement:':<30} {result.improvement_percent:.1f}%")
    print(f"{'Composition inequality:':<30} {'PASS' if result.composition_holds else 'FAIL'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
