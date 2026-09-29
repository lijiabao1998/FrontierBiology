#!/usr/bin/env python3
"""Bounded regression checks for evidence correction, not new research."""
import copy
import importlib.util
import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


verifier = load("independent_leakage_verifier")
manifest = load("verify_manifest")


class EvidenceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads((ROOT / "results/r1/bio001_r1_results.json").read_text(encoding="utf-8"))

    def test_negative_result_matches_independent_reconstruction(self):
        result = verifier.verify(self.result)
        self.assertEqual(result["verdict"], "INDEPENDENT_VERIFICATION_MATCH")
        self.assertTrue(all(value is True for value in result["checks"].values()))
        self.assertEqual(result["independent_counts"]["donor_test_size"], 1105)
        self.assertTrue(result["independent_counts"]["donor_no_row_overlap"])

    def test_corrupted_metrics_and_pass_flags_are_rejected(self):
        # Each mutation represents a prior review failure, not a mirrored unit check.
        mutations = [
            ("verdict", "LEAKAGE_DEMO_PASS"),
            ("synthetic_demo/donor_leak_gap_vs_donorclean", 0.905),
            ("synthetic_demo/frozen_train_size", 5230),
            ("synthetic_demo/group_leak_gap", 0.8776),
            ("pre_registered_checks/E2_group_leak_gap_ge_0.2", True),
            ("synthetic_demo/donor_design/n_test", 4509),
            ("synthetic_demo/leaky_pearson", 0.99),
        ]
        for path, value in mutations:
            with self.subTest(path=path):
                altered = copy.deepcopy(self.result)
                fields = path.split("/")
                target = altered
                for field in fields[:-1]:
                    target = target[field]
                target[fields[-1]] = value
                self.assertEqual(verifier.verify(altered)["verdict"], "INDEPENDENT_VERIFICATION_MISMATCH")

    def test_verifier_cli_rejects_missing_or_false_evidence(self):
        # The preserved v4 verifier output is deliberately not a valid r1 result.
        result = subprocess.run([sys.executable, "-B", str(HERE / "independent_leakage_verifier.py"),
            "--check-only", "--result", str(ROOT / "results/r1/history/independent_leakage_verification.v4.json")],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn("INDEPENDENT_VERIFICATION_MISMATCH", result.stdout)

    def test_manifest_detects_changed_bytes(self):
        statuses = manifest.check_entries(["0" * 64 + " *result.json"], lambda name: b"changed\n")
        self.assertEqual(statuses, [("MISMATCH", "result.json")])
        self.assertEqual(manifest.check_entries([], lambda name: b""), [("EMPTY_MANIFEST", manifest.MANIFEST)])

    def test_r2_is_exploratory_and_has_no_pass_threshold(self):
        result = json.loads((ROOT / "results/r2/bio001_r2_ablation_results.json").read_text(encoding="utf-8"))
        self.assertEqual(result["verdict"], "EXPLORATORY_UNVERIFIED")
        self.assertIs(result["independently_verified"], False)
        self.assertIsNone(result["success_threshold"])
        self.assertNotIn("generated", result)

    def test_existing_producer_replays_preserve_exact_bytes(self):
        pairs = [(HERE / "bio001_split_evaluator.py", ROOT / "results/r1/bio001_r1_results.json"),
                 (HERE.parent / "bio001_r2_ablation/eval_model_ablation.py", ROOT / "results/r2/bio001_r2_ablation_results.json"),
                 (HERE / "independent_leakage_verifier.py", ROOT / "results/r1/independent_leakage_verification.json")]
        for script, output in pairs:
            with self.subTest(script=script.name):
                before = output.read_bytes()
                result = subprocess.run([sys.executable, "-B", str(script)], capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", errors="replace"))
                self.assertEqual(output.read_bytes(), before)
                self.assertNotIn(b"\r\n", before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
