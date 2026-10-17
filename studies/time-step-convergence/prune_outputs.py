#!/usr/bin/env python3
"""Remove unneeded generated fields while keeping alpha history and restart state."""

from __future__ import annotations

import argparse
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNS = (HERE / "runs").resolve()
REMOVABLE_FIELDS = ("U", "p", "p_rgh", "rho", "phi", "alphaPhi.water")


def time_directories(case: Path) -> list[Path]:
    found: list[tuple[float, Path]] = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0:
            found.append((value, path))
    return [path for _, path in sorted(found)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--include-latest",
        action="store_true",
        help="also prune the newest snapshot after a completed run",
    )
    args = parser.parse_args()
    if RUNS.parent != HERE.resolve() or RUNS.name != "runs":
        raise RuntimeError(f"refusing unexpected runs path: {RUNS}")
    removed_files = 0
    removed_bytes = 0
    for variant in ("coarse_40us", "baseline_20us", "fine_10us", "ultrafine_5us"):
        for case_name in ("rigid", "flexible"):
            case = RUNS / variant / case_name
            if not case.exists():
                continue
            snapshots = time_directories(case)
            targets = snapshots if args.include_latest else snapshots[:-1]
            for snapshot in targets:
                for field in REMOVABLE_FIELDS:
                    path = snapshot / field
                    if path.is_file():
                        removed_bytes += path.stat().st_size
                        path.unlink()
                        removed_files += 1
    print(f"removed {removed_files} auxiliary field files ({removed_bytes / 1024**2:.1f} MiB)")


if __name__ == "__main__":
    main()
