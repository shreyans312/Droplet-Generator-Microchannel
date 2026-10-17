#!/usr/bin/env python3
"""Wait for a CFD case, then measure droplets and update the raw comparison.

This never modifies the physics or labels an unconverged run validated.
"""

from __future__ import annotations

import argparse
import gzip
import shlex
import subprocess
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MEASURE = ROOT / "cases/flexible-rigid/measure.py"


def completed(log: Path) -> bool:
    if not log.exists():
        return False
    if log.suffix == ".gz":
        # Called after the solver pipeline has closed the gzip file.
        ending = b""
        with gzip.open(log, "rb") as handle:
            while block := handle.read(1024 * 1024):
                ending = (ending + block)[-4096:]
    else:
        with log.open("rb") as handle:
            handle.seek(max(0, log.stat().st_size - 4096))
            ending = handle.read()
    if b"FOAM FATAL" in ending or b"Floating point exception" in ending:
        raise RuntimeError(f"solver failed; inspect {log}")
    return b"\nEnd\n" in ending


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--parallel", action="store_true")
    parser.add_argument("--deadline-hours", type=float, default=24)
    parser.add_argument("--skip-comparison", action="store_true", help="do not overwrite baseline comparison for a tagged diagnostic case")
    args = parser.parse_args()
    case = args.case.resolve()
    if not (case / "case.json").exists():
        raise FileNotFoundError(case / "case.json")
    log = case / args.log
    deadline = time.monotonic() + 3600 * args.deadline_hours
    while not completed(log):
        if time.monotonic() > deadline:
            raise TimeoutError(f"case did not complete before deadline: {log}")
        time.sleep(30)
    print(f"solver complete: {case}", flush=True)
    if args.parallel:
        command = (
            "source /opt/openfoam13/etc/bashrc && "
            f"reconstructPar -case {shlex.quote(str(case))} -newTimes -noFunctionObjects"
        )
        subprocess.run(["/bin/bash", "-lc", command], check=True)
    subprocess.run(
        [
            "/usr/bin/python3", str(MEASURE), "--case", str(case),
            "--x-min-um", "450", "--x-max-um", "990",
            "--y-min-um", "-90", "--y-max-um", "90",
            "--measurement-x-um", "750",
        ],
        check=True,
    )
    with (case / "alpha_diagnostic.txt").open("w") as handle:
        subprocess.run(
            ["/usr/bin/python3", str(HERE / "diagnose_alpha.py"), str(case), "--time", "1e9"],
            stdout=handle,
            check=True,
        )
    subprocess.run(["/usr/bin/python3", str(HERE / "audit.py"), str(case)], check=True)
    if not args.skip_comparison:
        mesh = case.parent.name.removeprefix("runs-")
        subprocess.run([sys.executable, str(HERE / "compare.py"), "--mesh", mesh], check=True)


if __name__ == "__main__":
    main()
