#!/usr/bin/env python3
"""Conservative cleanup of verified duplicate CFD outputs, not unique results.

Default is a read-only plan. --apply verifies fields again before deleting any
historical rank directories. Global reconstructed fields and the newest rank
checkpoint remain available. Large logs are losslessly compressed and verified.
"""
import argparse
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np
from audit import scalar_values

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEMP = Path("/tmp/droplet-validation.qL9S2V")
CASES = [HERE / "runs-focused/unstretched_robust",
         HERE / "runs-focused/unstretched_robust_auto1_sharpened_compression"]
CASES += sorted((HERE / "static-drop").glob("*"))


def byte_size(path):
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file() and not p.is_symlink())


def field_values(path, cell_count):
    data = path.read_bytes()
    if b"class volScalarField" in data or re.search(rb"class\s+volScalarField", data):
        values = scalar_values(path)
        return np.repeat(values, cell_count) if len(values) == 1 else values
    match = re.search(rb"internalField\s+uniform\s+\(([^)]+)\)", data)
    if match:
        return np.tile(np.fromstring(match[1].decode(), sep=" "), (cell_count, 1))
    match = re.search(rb"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(", data)
    if not match:
        raise ValueError(f"unsupported vector field: {path}")
    count, start = int(match[1]), match.end()
    if count != cell_count:
        raise ValueError(f"cell count mismatch: {path}")
    if re.search(rb"format\s+binary", data[:start]):
        end = start + count * 3 * 8
        if data[end:end+1] != b")":
            raise ValueError(f"unexpected binary vector architecture: {path}")
        return np.frombuffer(data, dtype="=f8", count=count*3, offset=start).reshape((-1, 3))
    raise ValueError(f"nonuniform ASCII vector comparison not implemented: {path}")


def verify(case, time_name, ranks):
    root_time = case / time_name
    parts = [rank / time_name for rank in ranks]
    if not root_time.is_dir() or not all(p.is_dir() for p in parts):
        return False
    if not all((p / "alpha.water").is_file() for p in [root_time, *parts]):
        return False
    root_alpha = scalar_values(root_time / "alpha.water")
    count = len(root_alpha)
    part_counts = [len(scalar_values(p / "alpha.water")) for p in parts]
    if any(n <= 1 for n in part_counts) or sum(part_counts) != count:
        return False
    # All saved rank field names must have a nonempty reconstructed counterpart.
    names = {f.name for p in parts for f in p.iterdir() if f.is_file()}
    if any(not (root_time / name).is_file() or (root_time / name).stat().st_size == 0 for name in names):
        return False
    for name in ("alpha.water", "U", "p_rgh", "rho", "p"):
        if name not in names:
            continue
        global_values = field_values(root_time / name, count).ravel()
        parallel_values = np.concatenate([field_values(p / name, n).ravel() for p, n in zip(parts, part_counts)])
        if global_values.size != parallel_values.size or not np.array_equal(np.sort(global_values), np.sort(parallel_values)):
            return False
    return True


def duplicate_directories():
    selected = []
    for case in CASES:
        if not case.is_dir() or case.is_symlink():
            continue
        ranks = sorted([p for p in case.iterdir() if p.is_dir() and re.fullmatch(r"processor\d+", p.name)])
        if len(ranks) != 8:
            continue
        if any(p.is_symlink() for p in ranks):
            raise RuntimeError(f"symlinked processor directory: {case}")
        times = [p.name for p in ranks[0].iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name)]
        newest = max(times, key=float)
        for name in sorted(times, key=float):
            if float(name) == 0 or name == newest:
                continue
            if verify(case, name, ranks):
                selected.extend(rank / name for rank in ranks)
        print(f"Verified duplicate history: {case.name}; newest checkpoint {newest} retained", flush=True)
    return selected


def digest(stream):
    sha = hashlib.sha256()
    while chunk := stream.read(1024*1024):
        sha.update(chunk)
    return sha.hexdigest()


def main(apply=False):
    before = shutil.disk_usage(HERE).free
    directories = duplicate_directories()
    logs = sorted(p for p in HERE.rglob("log.*") if p.is_file() and not p.is_symlink() and p.suffix != ".gz" and p.stat().st_size >= 5*1024*1024 and not p.with_name(p.name+".gz").exists())
    caches = [HERE / "__pycache__", HERE / "curvature/Make/linux64GccDPInt32Opt", HERE / "curvatureSolver/Make/linux64GccDPInt32Opt"]
    caches = [p for p in caches if p.is_dir()]
    report = {"applied": apply, "duplicate_directories": [str(p) for p in directories],
              "duplicate_bytes": sum(byte_size(p) for p in directories),
              "logs_to_compress": [str(p) for p in logs],
              "uncompressed_log_bytes": sum(p.stat().st_size for p in logs),
              "regenerable_caches": [str(p) for p in caches],
              "temporary_directory": str(TEMP), "free_bytes_before": before}
    print(json.dumps({k:v for k,v in report.items() if k not in ("duplicate_directories", "logs_to_compress")}, indent=2), flush=True)
    if not apply:
        (HERE / "cleanup_plan.json").write_text(json.dumps(report, indent=2)+"\n")
        return
    status = HERE / "automatic_followup_status.json"
    if status.exists() and json.loads(status.read_text())["state"] in ("running_channel_diagnostic", "analysing_completed_candidate", "preparing_channel_diagnostic"):
        raise RuntimeError("active follow-up workflow; cleanup refused")
    # Only resolved rank time directories inside the explicitly selected cases.
    for path in directories:
        case = path.parent.parent
        if not path.is_dir():
            continue
        # The plan was already deeply verified; check original global fields
        # still exist immediately before deleting each rank copy.
        if not (case / path.name / "alpha.water").is_file():
            raise RuntimeError(f"reconstructed result disappeared: {path}")
        if path.is_symlink() or case not in CASES:
            raise RuntimeError(f"invalid cleanup target: {path}")
        shutil.rmtree(path)
    for path in logs:
        compressed = path.with_name(path.name + ".gz")
        with path.open("rb") as source, gzip.open(compressed, "xb", compresslevel=6) as target:
            shutil.copyfileobj(source, target, 1024*1024)
        with path.open("rb") as source, gzip.open(compressed, "rb") as target:
            if digest(source) != digest(target):
                raise RuntimeError(f"log compression verification failed: {path}")
        path.unlink()
        print(f"Losslessly compressed: {path}", flush=True)
    for path in caches:
        shutil.rmtree(path)
    if TEMP.exists():
        if TEMP.is_symlink():
            raise RuntimeError("symlinked temporary directory; cleanup refused")
        allowed = {"roshan2025.pdf", "roshan-page6.png"}
        if {p.name for p in TEMP.iterdir()} - allowed:
            raise RuntimeError("unexpected temporary files; leaving folder untouched")
        source = TEMP / "roshan2025.pdf"
        destination = ROOT / "references/roshan2025-stretchable-droplet-generation.pdf"
        if source.exists():
            if not destination.exists():
                shutil.copy2(source, destination)
            with source.open("rb") as original, destination.open("rb") as preserved:
                if digest(original) != digest(preserved):
                    raise RuntimeError("paper preservation verification failed")
        for path in TEMP.iterdir():
            path.unlink()
        TEMP.rmdir()
        report["paper_preserved_at"] = str(destination)
    report.update(free_bytes_after=shutil.disk_usage(HERE).free)
    report["reclaimed_bytes"] = report["free_bytes_after"] - before
    (HERE / "cleanup_report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ("reclaimed_bytes", "free_bytes_after")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    main(parser.parse_args().apply)
