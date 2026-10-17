#!/usr/bin/env python3
"""Run MPI CFD with compressed raw logs and a durable live status summary.

Numerical failures or low disk space request a checkpointed stop, not a new
calibration. Intended to run inside a persistent user service.
"""
import argparse
import gzip
import json
import re
import shutil
import subprocess
import time
from pathlib import Path


def run(case: Path) -> int:
    control = (case / "system/controlDict").read_text()
    target = float(re.search(r"\bendTime\s+([\d.eE+-]+)", control)[1])
    state = {"case": str(case), "state": "running", "target_time_s": target,
             "validation_pass": False, "warnings": []}
    start = time.monotonic()
    updated = start - 11
    bad_bounds = 0
    stop_requested = False
    proc = subprocess.Popen(["mpirun", "--oversubscribe", "-np", "8", "foamRun",
                             "-case", str(case), "-parallel"],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    status = case / "live_status.json"
    with gzip.open(case / "log.foamRun.monitored.gz", "wb", compresslevel=1) as raw, (case / "log.monitor").open("w", buffering=1) as live:
        for line in iter(proc.stdout.readline, b""):
            raw.write(line)
            decoded = line.decode(errors="replace").strip()
            if decoded.startswith("Time = "):
                state["time_s"] = float(decoded.split()[2].rstrip("s"))
            if decoded.startswith("deltaT = "):
                state["delta_t_s"] = float(decoded.split()[2])
            if "Min(alpha.water) =" in decoded:
                values = re.search(r"Min\(alpha.water\) = (\S+)\s+Max\(alpha.water\) = (\S+)", decoded)
                if values:
                    lo, hi = map(float, values.groups())
                    state.update(alpha_min=lo, alpha_max=hi)
                    bad_bounds = bad_bounds + 1 if lo < -.01 or hi > 1.01 else 0
            if "FOAM FATAL" in decoded or "Floating point exception" in decoded:
                state["warnings"].append(decoded)
                live.write(decoded + "\n")
            if decoded == "End":
                state["solver_end_marker"] = True
            now = time.monotonic()
            if now - updated >= 10:
                state["elapsed_wall_seconds"] = round(now - start, 1)
                state["disk_free_gb"] = round(shutil.disk_usage(case).free / 1e9, 3)
                if (state["disk_free_gb"] < .8 or bad_bounds >= 100) and not stop_requested:
                    reason = "low disk space" if state["disk_free_gb"] < .8 else "persistent unbounded phase fraction"
                    state["warnings"].append(reason)
                    subprocess.run(["foamDictionary", str(case / "system/controlDict"), "-entry", "stopAt", "-set", "writeNow"], check=True, stdout=subprocess.DEVNULL)
                    state["state"] = "checkpoint_stop_requested"
                    stop_requested = True
                status.write_text(json.dumps(state, indent=2) + "\n")
                live.write(json.dumps(state) + "\n")
                raw.flush()
                updated = now
    code = proc.wait()
    state["exit_code"] = code
    reached_target = abs(state.get("time_s", float("inf")) - target) <= 1e-8
    state["state"] = "completed" if code == 0 and state.get("solver_end_marker") and reached_target and not stop_requested else "stopped_or_failed"
    state["elapsed_wall_seconds"] = round(time.monotonic() - start, 1)
    status.write_text(json.dumps(state, indent=2) + "\n")
    return 0 if state["state"] == "completed" else code or 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    raise SystemExit(run(args.case.resolve()))
