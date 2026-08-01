import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import types
import unittest


try:
    import femm  # noqa: F401
except ModuleNotFoundError:
    sys.modules.setdefault("femm", types.ModuleType("femm"))

import femm_spm_template as femm_model
import generate_readme_assets as readme_assets
from motor_config import (
    canonical_config_sha256,
    config_to_dict,
    load_motor_config,
    physical_model_fingerprint,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ResolvedConfigProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.config = load_motor_config(REPO_ROOT / "motor_config.json")
        self.config_hash = canonical_config_sha256(self.config)
        self.physical_hash = physical_model_fingerprint(
            self.config.machine,
            self.config.materials,
            self.config.conventions,
        )

    def _payload(self) -> dict:
        return {
            "shared_config_sha256": self.config_hash,
            "physical_model_sha256": self.physical_hash,
            "shared_config": config_to_dict(self.config),
        }

    def _write_payload(self, payload: dict) -> Path:
        destination = self.directory / "resolved_femm_config.json"
        destination.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )
        return destination

    def test_load_resolved_config_recomputes_both_hashes(self):
        path = self._write_payload(self._payload())

        actual_config_hash, actual_physical_hash = (
            readme_assets._load_resolved_config(path, "test result")
        )

        self.assertEqual(actual_config_hash, self.config_hash)
        self.assertEqual(actual_physical_hash, self.physical_hash)

    def test_load_resolved_config_rejects_tampered_shared_config(self):
        payload = copy.deepcopy(self._payload())
        payload["shared_config"]["sweep"]["requested_torque_nm"] += 0.25
        path = self._write_payload(payload)

        with self.assertRaisesRegex(ValueError, "configuration hash mismatch"):
            readme_assets._load_resolved_config(path, "test result")

    def test_load_resolved_config_rejects_tampered_declared_config_hash(self):
        payload = self._payload()
        payload["shared_config_sha256"] = "0" * 64
        path = self._write_payload(payload)

        with self.assertRaisesRegex(ValueError, "configuration hash mismatch"):
            readme_assets._load_resolved_config(path, "test result")

    def test_load_resolved_config_rejects_tampered_physical_hash(self):
        payload = self._payload()
        payload["physical_model_sha256"] = "0" * 64
        path = self._write_payload(payload)

        with self.assertRaisesRegex(ValueError, "physical-model hash mismatch"):
            readme_assets._load_resolved_config(path, "test result")

    def test_require_matching_resolved_config_accepts_matching_physical_model(self):
        path = self._write_payload(self._payload())

        actual = readme_assets._require_matching_resolved_config(
            path,
            expected_physical_hash=self.physical_hash,
            description="test result",
        )

        self.assertEqual(actual, path)

    def test_require_matching_resolved_config_rejects_other_physical_model(self):
        path = self._write_payload(self._payload())

        with self.assertRaisesRegex(ValueError, "different FEMM physical model"):
            readme_assets._require_matching_resolved_config(
                path,
                expected_physical_hash="f" * 64,
                description="test result",
            )


class ResultFileProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        config = load_motor_config(REPO_ROOT / "motor_config.json")
        self.physical_hash = physical_model_fingerprint(
            config.machine,
            config.materials,
            config.conventions,
        )
        self.result = self.directory / "result.csv"
        self.result.write_text("angle,torque\n0,1.25\n", encoding="utf-8")
        self.manifest = self.directory / "femm_run_manifest.json"
        self.manifest.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "runs": [
                        {
                            "physical_model_sha256": self.physical_hash,
                            "analyses": ["validate"],
                            "arguments": {
                                "out": str(self.directory),
                                "mesh_level": "medium",
                            },
                            "files": [
                                {
                                    "path": self.result.name,
                                    "sha256": _sha256(self.result),
                                    "size_bytes": self.result.stat().st_size,
                                }
                            ],
                        }
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    def test_require_provenanced_file_accepts_exact_file_and_model_hash(self):
        actual = readme_assets._require_provenanced_file(
            self.result, self.physical_hash
        )

        self.assertEqual(actual, self.manifest)

    def test_require_provenanced_file_rejects_other_physical_model(self):
        with self.assertRaisesRegex(ValueError, "No successful FEMM run"):
            readme_assets._require_provenanced_file(self.result, "f" * 64)

    def test_require_provenanced_file_rejects_file_modified_after_run(self):
        self.result.write_text("angle,torque\n0,9.99\n", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "No successful FEMM run"):
            readme_assets._require_provenanced_file(
                self.result, self.physical_hash
            )

    def test_require_provenanced_file_enforces_analysis_arguments(self):
        actual = readme_assets._require_provenanced_file(
            self.result,
            self.physical_hash,
            required_analyses={"validate"},
            expected_arguments={"mesh_level": "medium"},
        )
        self.assertEqual(actual, self.manifest)

        with self.assertRaisesRegex(ValueError, "required analysis contract"):
            readme_assets._require_provenanced_file(
                self.result,
                self.physical_hash,
                required_analyses={"validate"},
                expected_arguments={"mesh_level": "fine"},
            )

    def test_require_provenanced_file_rejects_moved_output_directory(self):
        payload = json.loads(self.manifest.read_text(encoding="utf-8"))
        payload["runs"][0]["arguments"]["out"] = str(
            self.directory / "original-output"
        )
        self.manifest.write_text(json.dumps(payload), encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "required analysis contract"):
            readme_assets._require_provenanced_file(
                self.result, self.physical_hash
            )


class FemmRunManifestWriterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.config = load_motor_config(REPO_ROOT / "motor_config.json")
        self.args = argparse.Namespace(
            analysis=["validate", "convergence"],
            config=Path("motor_config.json"),
            show_femm=False,
        )

    def test_writer_records_only_created_or_replaced_files_and_appends_runs(self):
        unchanged = self.directory / "unchanged.csv"
        replaced = self.directory / "replaced.csv"
        unchanged.write_text("unchanged\n", encoding="utf-8")
        replaced.write_text("before\n", encoding="utf-8")
        before = femm_model._snapshot_output_files(self.directory)

        replaced.write_text("after with different size\n", encoding="utf-8")
        created = self.directory / "nested" / "created.csv"
        created.parent.mkdir()
        created.write_text("new result\n", encoding="utf-8")
        manifest = femm_model.write_femm_run_manifest(
            self.directory,
            before=before,
            started_at_utc="2026-08-01T00:00:00+00:00",
            analyses={"validate", "convergence"},
            args=self.args,
            config=self.config,
        )

        payload = json.loads(manifest.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(len(payload["runs"]), 1)
        first_run = payload["runs"][0]
        self.assertEqual(first_run["analyses"], ["convergence", "validate"])
        self.assertEqual(
            first_run["shared_config_sha256"],
            canonical_config_sha256(self.config),
        )
        self.assertEqual(
            first_run["physical_model_sha256"],
            physical_model_fingerprint(
                self.config.machine,
                self.config.materials,
                self.config.conventions,
            ),
        )
        self.assertEqual(first_run["arguments"]["config"], "motor_config.json")
        files = {item["path"]: item for item in first_run["files"]}
        self.assertEqual(set(files), {"nested/created.csv", "replaced.csv"})
        self.assertEqual(files["replaced.csv"]["sha256"], _sha256(replaced))
        self.assertEqual(files["nested/created.csv"]["sha256"], _sha256(created))
        self.assertNotIn("femm_run_manifest.json", files)
        self.assertNotIn("unchanged.csv", files)

        second_before = femm_model._snapshot_output_files(self.directory)
        second_result = self.directory / "second.csv"
        second_result.write_text("second run\n", encoding="utf-8")
        femm_model.write_femm_run_manifest(
            self.directory,
            before=second_before,
            started_at_utc="2026-08-01T00:01:00+00:00",
            analyses={"field"},
            args=self.args,
            config=self.config,
        )

        updated = json.loads(manifest.read_text(encoding="utf-8"))
        self.assertEqual(len(updated["runs"]), 2)
        self.assertEqual(updated["runs"][1]["analyses"], ["field"])
        self.assertEqual(
            [item["path"] for item in updated["runs"][1]["files"]],
            ["second.csv"],
        )

    def test_writer_rejects_a_run_without_changed_outputs(self):
        unchanged = self.directory / "unchanged.csv"
        unchanged.write_text("unchanged\n", encoding="utf-8")
        before = femm_model._snapshot_output_files(self.directory)

        with self.assertRaisesRegex(RuntimeError, "did not create or replace"):
            femm_model.write_femm_run_manifest(
                self.directory,
                before=before,
                started_at_utc="2026-08-01T00:00:00+00:00",
                analyses={"validate"},
                args=self.args,
                config=self.config,
            )


if __name__ == "__main__":
    unittest.main()
