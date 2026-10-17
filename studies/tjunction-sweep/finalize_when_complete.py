#!/usr/bin/env python3
"""Finalize the 18-point medium sweep after every solver case completes."""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_ROOT = HERE / "runs-medium"
EXPECTED_CASES = 18
POLL_SECONDS = 60


def completed_count() -> int:
    return sum(1 for _ in RUN_ROOT.glob("*/RUN_COMPLETE"))


def main() -> None:
    while completed_count() < EXPECTED_CASES:
        time.sleep(POLL_SECONDS)
    subprocess.run(
        [sys.executable, str(HERE / "analyze.py"), "--mesh", "medium"],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        [sys.executable, str(ROOT / "studies" / "raj2016-comparison" / "generate.py")],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        [sys.executable, str(HERE / "progress.py")],
        cwd=ROOT,
        check=True,
    )
    (HERE / "SWEEP_FINALIZED").write_text(
        "18 solver cases complete; measurements and publication plots regenerated.\n"
    )


if __name__ == "__main__":
    main()
