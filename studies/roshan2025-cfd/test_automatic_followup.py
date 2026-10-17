import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from automatic_followup import METHODS, candidates, screen
from finalize import completed
from monitor_solver import run
import generate as model_generator


class FollowupTests(unittest.TestCase):
    def test_generated_cases_preserve_materials_and_checkpoint_start(self):
        with tempfile.TemporaryDirectory(prefix="cfd-generator-test-") as directory:
            root = Path(directory) / "root"
            study = root / "studies/roshan2025-cfd"
            common = root / "cases/flexible-rigid/common"
            study.mkdir(parents=True)
            common.mkdir(parents=True)
            source = model_generator.HERE.parents[1] / "cases/flexible-rigid/common"
            for filename in ("fvSchemes", "fvSolution"):
                (common / filename).write_text((source / filename).read_text())
            with patch.object(model_generator, "HERE", study):
                case = model_generator.generate("unstretched", "focused", 1.8, tag="test",
                                                surface_force="sharpened", start_from_latest=True,
                                                capillary_safety=1)
            control = (case / "system/controlDict").read_text()
            self.assertIn("startFrom latestTime;", control)
            self.assertIn("solver incompressibleVoFCurvature;", control)
            self.assertIn("boundaryFlux", control)
            self.assertIn("rho 1630;", (case / "constant/physicalProperties.water").read_text())
            self.assertEqual(json.loads((case / "case.json").read_text())["cell_count"], 42192)

    def test_candidates_do_not_use_experimental_diameter(self):
        rows = [{"case": "fine_capped", "time_s": .003,
                 "pressure_error_relative_volume_radius_percent": -16,
                 "maximum_speed_m_s": .002}]
        for index, name in enumerate(METHODS):
            rows.extend([{"case": name, "time_s": 0, "phase_volume_m3": 1e-12},
                         {"case": name, "time_s": .003, "phase_volume_m3": 1e-12,
                          "pressure_error_relative_volume_radius_percent": -4-index,
                          "maximum_speed_m_s": .00005*(index+1)}])
        choice = candidates(rows)[0]
        self.assertEqual(choice["benchmark"], "fine_smoothed")
        self.assertTrue(choice["static_verification_pass"])
        rows[-1]["phase_volume_m3"] = 1.1e-12
        self.assertNotIn("fine_mplic_smoothed", [r["benchmark"] for r in candidates(rows)])

    def test_frequency_failure_cannot_be_hidden_by_size_agreement(self):
        stats = {"count": 20, "diameter_mean_um": 88, "frequency_hz": 28, "diameter_cv_percent": 2}
        report = {"post_washout_statistics": {"axial_diameter_um": stats},
                  "solver_step_conservation": {"residual_percent_of_injected_phase": .01},
                  "last_time_blocks": [stats]*3, "completed_output": True}
        verdict = screen(report)
        self.assertFalse(verdict["screening_pass"])
        self.assertFalse(verdict["checks"]["frequency_agreement"])
        self.assertFalse(verdict["validation_pass"])

    def test_monitor_and_compressed_completion(self):
        with tempfile.TemporaryDirectory(prefix="cfd-monitor-test-") as directory:
            case = Path(directory)
            (case / "system").mkdir()
            (case / "system/controlDict").write_text("endTime 1.8;\n")
            output = b"deltaT = 1e-6\nTime = 1.8s\nPhase-1 volume fraction = 0.1 Min(alpha.water) = 0 Max(alpha.water) = 1\nEnd\n"
            process = SimpleNamespace(stdout=io.BytesIO(output), wait=lambda: 0)
            with patch("monitor_solver.subprocess.Popen", return_value=process):
                self.assertEqual(run(case), 0)
            self.assertEqual(json.loads((case / "live_status.json").read_text())["time_s"], 1.8)
            self.assertTrue(completed(case / "log.foamRun.monitored.gz"))
            with gzip.open(case / "log.foamRun.monitored.gz", "rb") as log:
                self.assertEqual(log.read(), output)


if __name__ == "__main__":
    unittest.main()
