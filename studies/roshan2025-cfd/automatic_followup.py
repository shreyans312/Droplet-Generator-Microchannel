#!/usr/bin/env python3
"""Bounded, uncalibrated numerical-method follow-up after verification tests.

Runs at most two independent candidates from the same saved rigid checkpoint.
Passes here are screening gates, never a claim of publication-ready validation.
Physics and geometry are not tuned to experiment. Old cases are never changed.
"""
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STATUS = HERE / "automatic_followup_status.json"
METHODS = {"fine_smoothed": ("smoothed", "compression"),
           "fine_sharpened": ("sharpened", "compression"),
           "fine_mplic_smoothed": ("smoothed", "mplic")}


def save(state):
    STATUS.write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state), flush=True)


def candidates(results):
    final = {r["case"]: r for r in results if r.get("time_s") == .003}
    stock = final["fine_capped"]
    improved = []
    for name in METHODS:
        row = final[name]
        pressure = abs(row["pressure_error_relative_volume_radius_percent"])
        speed = row["maximum_speed_m_s"]
        initial = next(r for r in results if r.get("case") == name and r.get("time_s") == 0)
        mass_error = abs(row["phase_volume_m3"] / initial["phase_volume_m3"] - 1)
        if pressure < abs(stock["pressure_error_relative_volume_radius_percent"]) and speed < stock["maximum_speed_m_s"] and mass_error < 1e-5:
            # Physical verification metrics only: no experimental size fitting.
            score = pressure / 5 + speed / .0001
            improved.append({"benchmark": name, "surface_force": METHODS[name][0],
                             "advection": METHODS[name][1], "score": score,
                             "pressure_error_percent": pressure, "maximum_speed_m_s": speed,
                             "static_verification_pass": pressure <= 5 and speed <= .0001})
    return sorted(improved, key=lambda r: r["score"])


def screen(audit):
    stats = audit["post_washout_statistics"]["axial_diameter_um"]
    mass = audit["solver_step_conservation"]
    blocks = audit["last_time_blocks"]
    means = [r["diameter_mean_um"] for r in blocks if r["count"] >= 2]
    freqs = [r["frequency_hz"] for r in blocks if r["count"] >= 2]
    diameter = stats["diameter_mean_um"]
    frequency = stats["frequency_hz"]
    checks = {"completed": audit["completed_output"], "enough_droplets": stats["count"] >= 6,
              "mass_conservation": mass is not None and abs(mass["residual_percent_of_injected_phase"]) <= .5,
              "low_size_cv": stats["diameter_cv_percent"] is not None and stats["diameter_cv_percent"] <= 5,
              "block_size_stationarity": len(means) == 3 and (max(means)-min(means))/sum(means)*3 <= .05,
              "block_frequency_stationarity": len(freqs) == 3 and (max(freqs)-min(freqs))/sum(freqs)*3 <= .05,
              "provisional_axial_size_agreement": diameter is not None and abs(diameter/88-1) <= .1,
              "frequency_agreement": frequency is not None and abs(frequency/14-1) <= .1}
    return {"screening_pass": all(checks.values()), "checks": checks,
            "diameter_um": diameter, "frequency_hz": frequency,
            "validation_pass": False, "reason_not_validation": "Independent spatial/temporal convergence and experimental diameter/geometry correspondence remain required."}


def main():
    state = {"state": "waiting_for_verification_tests", "validation_pass": False, "runs": []}
    save(state)
    deadline = time.monotonic() + 3600
    def benchmark_finished(name):
        case = HERE / "static-drop" / name
        log = case / "log.reconstruct"
        return (case / "0.003/U").exists() and log.exists() and "\nEnd\n" in log.read_text()[-4096:]
    while not all(benchmark_finished(name) for name in METHODS):
        if time.monotonic() > deadline:
            state.update(state="verification_tests_incomplete", reason="No channel candidate started; inspect benchmark service logs.")
            save(state)
            return
        time.sleep(10)
    subprocess.run(["/usr/bin/python3", str(HERE / "static_drop.py"), "analyse"], check=True, stdout=subprocess.DEVNULL)
    choices = candidates(json.loads((HERE / "static-drop/results.json").read_text()))[:2]
    state["candidates"] = choices
    if not choices:
        state.update(state="no_improved_candidate", reason="All tested numerical methods failed the improvement gate; do not fit a diameter factor.")
        save(state)
        return
    from generate import generate
    original = HERE / "runs-focused/unstretched_robust"
    for index, choice in enumerate(choices, 1):
        if shutil.disk_usage(HERE).free < 2.5e9:
            state.update(state="stopped_for_disk_space", reason="Preserving original results; need at least 2.5 GB free for another candidate.")
            save(state)
            return
        state.update(state="preparing_channel_diagnostic", active_candidate=choice)
        save(state)
        case = generate("unstretched", "focused", 1.8, tag=f"auto{index}_{choice['surface_force']}_{choice['advection']}",
                        alpha_scheme=choice["advection"], start_from_latest=True,
                        capillary_safety=1.0, surface_force=choice["surface_force"])
        shutil.copytree(original / "constant/polyMesh", case / "constant/polyMesh")
        shutil.copytree(original / "1.2", case / "1.2")
        metadata = json.loads((case / "case.json").read_text())
        metadata.update(restart_checkpoint_s=1.2, restart_source=str(original), analysis_start_time_s=1.4,
                        verification_benchmark=choice, purpose="Uncalibrated method-isolation diagnostic; not a validated production result.")
        (case / "case.json").write_text(json.dumps(metadata, indent=2) + "\n")
        for utility in ("checkMesh", "decomposePar"):
            command = [utility, "-case", str(case)] + (["-latestTime"] if utility == "decomposePar" else [])
            with (case / f"log.{utility}").open("w") as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        state.update(state="running_channel_diagnostic", active_case=str(case))
        save(state)
        code = subprocess.run([sys.executable, str(HERE / "monitor_solver.py"), str(case)]).returncode
        if code:
            state.update(state="solver_stopped_or_failed", reason=f"Candidate returned {code}; original case preserved.")
            save(state)
            return
        state["state"] = "analysing_completed_candidate"
        save(state)
        with (case / "log.finalize").open("w") as log:
            subprocess.run([sys.executable, str(HERE / "finalize.py"), "--case", str(case), "--log", "log.foamRun.monitored.gz", "--parallel", "--skip-comparison"], check=True, stdout=log, stderr=subprocess.STDOUT)
        verdict = screen(json.loads((case / "audit.json").read_text()))
        (case / "screening.json").write_text(json.dumps(verdict, indent=2) + "\n")
        state["runs"].append({"case": str(case), **verdict})
        save(state)
        if verdict["screening_pass"] and choice["static_verification_pass"]:
            state.update(state="screening_passed_convergence_still_required")
            save(state)
            return
    state.update(state="tested_candidates_need_further_correction", reason="Bounded diagnostic iteration completed; inspect raw audits before choosing further geometry or solver changes.")
    save(state)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        state = json.loads(STATUS.read_text()) if STATUS.exists() else {}
        state.update(state="workflow_failed", error=f"{type(error).__name__}: {error}")
        save(state)
        raise
