#!/usr/bin/env python3
"""Run generated cases and prepare raw top-view measurements."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


HERE = Path(__file__).resolve().parent
REMOVABLE_FIELDS = ("U", "p", "p_rgh", "rho", "phi", "alphaPhi.water")


def compress_log(path: Path) -> Path:
    """Compress a completed solver log without holding it in memory."""
    destination = path.with_name(path.name + ".gz")
    temporary = destination.with_name(destination.name + ".tmp")
    with path.open("rb") as source, gzip.open(temporary, "wb", compresslevel=6) as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    temporary.replace(destination)
    path.unlink()
    return destination


def prune_auxiliary_snapshots(case: Path) -> tuple[int, int]:
    """Keep alpha history for measurement and one complete restart snapshot."""
    snapshots: list[tuple[float, Path]] = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            time_value = float(path.name)
        except ValueError:
            continue
        if time_value > 0:
            snapshots.append((time_value, path))
    removed_files = 0
    removed_bytes = 0
    for _, snapshot in sorted(snapshots)[:-1]:
        for field in REMOVABLE_FIELDS:
            path = snapshot / field
            if path.is_file():
                removed_bytes += path.stat().st_size
                path.unlink()
                removed_files += 1
    return removed_files, removed_bytes


def selected_ids(group: str) -> list[str]:
    with (HERE / "matrix.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if group == "all":
        return [row["case_id"] for row in rows]
    if group in {"fig6_r05", "fig6_r10", "fig6_r15"}:
        ratio = {"fig6_r05": "0.5", "fig6_r10": "1.0", "fig6_r15": "1.5"}[group]
        return [row["case_id"] for row in rows if row["flow_ratio"] == ratio]
    requested = set(group.split(","))
    return [row["case_id"] for row in rows if row["case_id"] in requested]


def run(case_id: str, mesh_level: str, prepare_only: bool = False) -> tuple[str, bool, str]:
    case = HERE / f"runs-{mesh_level}" / case_id
    if not case.exists():
        return case_id, False, "case has not been generated"
    completion_marker = case / "RUN_COMPLETE"
    if completion_marker.exists() and not prepare_only:
        return case_id, True, "already complete"
    preparation_commands = (
        (["blockMesh"], "log.blockMesh"),
        (["checkMesh"], "log.checkMesh"),
        (["setFields"], "log.setFields"),
    )
    commands = preparation_commands if prepare_only else preparation_commands + (
        (["foamRun", "-solver", "incompressibleVoF"], "log.foamRun"),
    )
    for command, log_name in commands:
        log_path = case / log_name
        with log_path.open("w") as log:
            result = subprocess.run(command, cwd=case, env=os.environ.copy(), stdout=log, stderr=subprocess.STDOUT)
        reported_log = log_name
        if log_name == "log.foamRun":
            reported_log = compress_log(log_path).name
        if result.returncode:
            return case_id, False, f"{' '.join(command)} failed; see {reported_log}"
    if "Mesh OK" not in (case / "log.checkMesh").read_text(errors="ignore"):
        return case_id, False, "checkMesh did not report Mesh OK"
    if prepare_only:
        return case_id, True, "mesh and initial fields prepared"
    removed_files, removed_bytes = prune_auxiliary_snapshots(case)
    completion_marker.write_text(
        f"complete\npruned_files={removed_files}\npruned_bytes={removed_bytes}\n"
    )
    return case_id, True, f"complete; pruned {removed_bytes / 1024**2:.1f} MiB"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", default="fig6_r05")
    parser.add_argument("--mesh", choices=("screening", "medium", "fine"), default="screening")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [
            pool.submit(run, case_id, args.mesh, args.prepare_only)
            for case_id in selected_ids(args.group)
        ]
        for future in as_completed(futures):
            case_id, ok, message = future.result()
            print(f"{case_id}: {'OK' if ok else 'FAILED'} - {message}", flush=True)


if __name__ == "__main__":
    main()
