#!/usr/bin/env python3
"""Generate live SVG/CSV progress reports for the 18-point Figure 6 sweep."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUN_ROOT = HERE / "runs-medium"
TIME_PATTERN = re.compile(rb"(?m)^Time = ([0-9.eE+-]+)s$")
COLOURS = {"ready": "#15803d", "running": "#2563eb", "queued": "#9ca3af", "failed": "#dc2626"}


def latest_time(log: Path) -> float:
    if not log.exists():
        return 0.0
    with log.open("rb") as handle:
        size = handle.seek(0, 2)
        handle.seek(max(0, size - 1024 * 1024))
        matches = TIME_PATTERN.findall(handle.read())
    return float(matches[-1]) if matches else 0.0


def rows() -> list[dict[str, object]]:
    with (HERE / "matrix.csv").open(newline="") as handle:
        matrix = list(csv.DictReader(handle))
    result: list[dict[str, object]] = []
    for item in matrix:
        case = RUN_ROOT / item["case_id"]
        metadata = json.loads((case / "case.json").read_text())
        target = float(metadata["end_time_s"])
        if (case / "RUN_COMPLETE").exists():
            state = "ready"
            current = target
        elif (case / "log.foamRun").exists():
            state = "running"
            current = latest_time(case / "log.foamRun")
        elif (case / "log.foamRun.gz").exists():
            state = "failed"
            current = 0.0
        else:
            state = "queued"
            current = 0.0
        result.append({
            "case_id": item["case_id"],
            "flow_ratio": float(item["flow_ratio"]),
            "kp": float(item["kp"]),
            "state": state,
            "simulated_time_s": current,
            "target_time_s": target,
            "progress_percent": min(100.0, 100 * current / target),
        })
    return result


def write_csv(data: list[dict[str, object]]) -> None:
    with (HERE / "sweep_progress.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def write_point_progress(data: list[dict[str, object]]) -> None:
    width, height = 1400, 880
    left, right, top, bottom = 235, 210, 110, 70
    plot_width = width - left - right
    row_height = (height - top - bottom) / len(data)
    ready = sum(row["state"] == "ready" for row in data)
    running = sum(row["state"] == "running" for row in data)
    queued = sum(row["state"] == "queued" for row in data)
    failed = sum(row["state"] == "failed" for row in data)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111827}.grid{stroke:#d1d5db;stroke-width:1}.minor{stroke:#e5e7eb;stroke-width:1}</style>',
        f'<text x="{width/2:.1f}" y="34" font-size="24" font-weight="600" text-anchor="middle">Figure 6 CFD sweep progress</text>',
        f'<text x="{width/2:.1f}" y="62" font-size="15" text-anchor="middle">Ready {ready}/18  |  Running {running}  |  Queued {queued}  |  Failed {failed}</text>',
        f'<text x="{width/2:.1f}" y="86" font-size="13" text-anchor="middle">A point is ready only after its solver reaches the configured final time.</text>',
    ]
    for tick in range(0, 101, 10):
        px = left + plot_width * tick / 100
        klass = "grid" if tick % 20 == 0 else "minor"
        parts += [
            f'<line class="{klass}" x1="{px:.1f}" y1="{top}" x2="{px:.1f}" y2="{height-bottom}"/>',
            f'<text x="{px:.1f}" y="{height-bottom+26}" font-size="13" text-anchor="middle">{tick}%</text>',
        ]
    for index, row in enumerate(data):
        cy = top + (index + 0.5) * row_height
        progress = float(row["progress_percent"])
        state = str(row["state"])
        colour = COLOURS[state]
        label = f'r={float(row["flow_ratio"]):.1f}, Kp={float(row["kp"]):+.5f}'
        parts += [
            f'<text x="{left-14}" y="{cy+5:.1f}" font-size="13" text-anchor="end">{label}</text>',
            f'<rect x="{left}" y="{cy-9:.1f}" width="{plot_width}" height="18" rx="3" fill="#f3f4f6"/>',
            f'<rect x="{left}" y="{cy-9:.1f}" width="{plot_width*progress/100:.1f}" height="18" rx="3" fill="{colour}"/>',
            f'<text x="{left+plot_width+10}" y="{cy+5:.1f}" font-size="12">{state}: {float(row["simulated_time_s"]):.3f}/{float(row["target_time_s"]):.1f} s</text>',
        ]
    parts.append('</svg>')
    (HERE / "sweep_point_progress.svg").write_text("\n".join(parts) + "\n")


def write_counts(data: list[dict[str, object]]) -> None:
    states = ("ready", "running", "queued", "failed")
    counts = {state: sum(row["state"] == state for row in data) for state in states}
    width, height = 760, 520
    left, right, top, bottom = 95, 50, 90, 85
    plot_height = height - top - bottom
    bar_width = 105
    gap = 42
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111827}.grid{stroke:#d1d5db;stroke-width:1}.axis{stroke:#111827;stroke-width:2}</style>',
        '<text x="380" y="34" font-size="23" font-weight="600" text-anchor="middle">CFD data-point readiness</text>',
        f'<text x="380" y="62" font-size="15" text-anchor="middle">{counts["ready"]} of 18 points currently ready for plotting</text>',
    ]
    for tick in range(0, 19, 3):
        py = height - bottom - plot_height * tick / 18
        parts += [
            f'<line class="grid" x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}"/>',
            f'<text x="{left-12}" y="{py+5:.1f}" font-size="13" text-anchor="end">{tick}</text>',
        ]
    parts += [
        f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>',
    ]
    for index, state in enumerate(states):
        x = left + 45 + index * (bar_width + gap)
        count = counts[state]
        y = height - bottom - plot_height * count / 18
        parts += [
            f'<rect x="{x}" y="{y:.1f}" width="{bar_width}" height="{height-bottom-y:.1f}" rx="5" fill="{COLOURS[state]}"/>',
            f'<text x="{x+bar_width/2:.1f}" y="{y-12:.1f}" font-size="22" font-weight="600" text-anchor="middle">{count}</text>',
            f'<text x="{x+bar_width/2:.1f}" y="{height-bottom+29}" font-size="14" text-anchor="middle">{state.capitalize()}</text>',
        ]
    parts.append('</svg>')
    (HERE / "sweep_readiness_counts.svg").write_text("\n".join(parts) + "\n")


def main() -> None:
    data = rows()
    write_csv(data)
    write_point_progress(data)
    write_counts(data)
    counts = {state: sum(row["state"] == state for row in data) for state in COLOURS}
    print(" ".join(f"{state}={count}" for state, count in counts.items()))


if __name__ == "__main__":
    main()
