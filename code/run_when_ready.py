#!/usr/bin/env python3
"""
Auto-run benchmark when Tuna-9 becomes available.

Polls backend status every 30 seconds. When it goes idle,
automatically runs the full benchmark suite.

Usage:
    python code/run_when_ready.py
    python code/run_when_ready.py --backend "Tuna-9" --poll 30
"""

import argparse
import subprocess
import sys
import time
import os


def check_status(backend_name):
    """Check backend status via QI API."""
    try:
        from qiskit_quantuminspire.qi_provider import QIProvider
        provider = QIProvider()
        for b in provider.backends():
            if b.name == backend_name:
                return b.status.value
        return "not_found"
    except Exception as e:
        return f"error: {e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", default="Tuna-9")
    parser.add_argument("--poll", type=int, default=30, help="Poll interval (seconds)")
    parser.add_argument("--shots", type=int, default=4096)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    python = sys.executable
    script_dir = os.path.dirname(os.path.abspath(__file__))
    benchmark = os.path.join(script_dir, "hardware_benchmark.py")

    print(f"Waiting for {args.backend} to become idle...")
    print(f"Polling every {args.poll} seconds.")
    print(f"Will run: {python} {benchmark} --backend '{args.backend}' "
          f"--shots {args.shots} --timeout {args.timeout}")
    print("-" * 60)

    attempt = 0
    while True:
        attempt += 1
        status = check_status(args.backend)
        ts = time.strftime("%H:%M:%S")
        print(f"[{ts}] Attempt {attempt}: {args.backend} status = {status}")

        if status == "idle":
            print(f"\n{'='*60}")
            print(f"{args.backend} is IDLE! Starting benchmark...")
            print(f"{'='*60}\n")

            cmd = [
                python, benchmark,
                "--backend", args.backend,
                "--shots", str(args.shots),
                "--timeout", str(args.timeout),
            ]
            result = subprocess.run(cmd, cwd=os.path.dirname(script_dir))

            if result.returncode == 0:
                print(f"\n{'='*60}")
                print("BENCHMARK COMPLETED SUCCESSFULLY!")
                print(f"{'='*60}")
            else:
                print(f"\nBenchmark exited with code {result.returncode}")
                print("Check results/ directory for partial output.")

            return result.returncode

        time.sleep(args.poll)


if __name__ == "__main__":
    main()
