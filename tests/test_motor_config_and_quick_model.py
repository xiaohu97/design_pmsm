import copy
from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest

import design_pmsm as quick
import motor_config as shared


def _install_optional_femm_import_stubs() -> None:
    """Allow pure adapter tests without plotting, progress, or COM packages."""
    try:
        import matplotlib.pyplot  # noqa: F401
    except ModuleNotFoundError:
        matplotlib = types.ModuleType("matplotlib")
        matplotlib.__path__ = []
        pyplot = types.ModuleType("matplotlib.pyplot")
        sys.modules.setdefault("matplotlib", matplotlib)
        sys.modules.setdefault("matplotlib.pyplot", pyplot)

    try:
        import tqdm  # noqa: F401
    except ModuleNotFoundError:
        tqdm_module = types.ModuleType("tqdm")
        tqdm_module.tqdm = lambda iterable=None, *args, **kwargs: iterable
        sys.modules.setdefault("tqdm", tqdm_module)

    try:
        import femm  # noqa: F401
    except ModuleNotFoundError:
        sys.modules.setdefault("femm", types.ModuleType("femm"))


class MotorConfigJsonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.path = shared.default_config_path()
        self.config = shared.load_motor_config(self.path)
        self.raw = shared.config_to_dict(self.config)

    def _load_text(self, text: str) -> shared.MotorConfig:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "motor.json"
            path.write_text(text, encoding="utf-8")
            return shared.load_motor_config(path)

    def test_default_json_round_trips_without_implicit_defaults(self) -> None:
        loaded_raw = json.loads(self.path.read_text(encoding="utf-8-sig"))
        self.assertEqual(shared.config_to_dict(self.config), loaded_raw)
        self.assertEqual(len(shared.canonical_config_sha256(self.config)), 64)

    def test_strict_json_rejects_nonfinite_numbers_and_invalid_syntax(self) -> None:
        valid_text = self.path.read_text(encoding="utf-8-sig")
        for token in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(token=token):
                text = valid_text.replace('"dc_bus_v": 48.0', f'"dc_bus_v": {token}', 1)
                with self.assertRaisesRegex(shared.ConfigError, "Non-finite JSON number"):
                    self._load_text(text)

        invalid = valid_text.replace('"schema_version": 1,', '"schema_version": 1,,', 1)
        with self.assertRaisesRegex(shared.ConfigError, "Invalid JSON"):
            self._load_text(invalid)

    def test_strict_json_rejects_duplicate_object_names(self) -> None:
        valid_text = self.path.read_text(encoding="utf-8-sig").lstrip()
        duplicate_name = '{\n  "name": "shadowed",' + valid_text[1:]
        with self.assertRaises(shared.ConfigError):
            self._load_text(duplicate_name)

    def test_unknown_and_missing_top_level_fields_are_rejected(self) -> None:
        unknown = copy.deepcopy(self.raw)
        unknown["unexpected"] = 1
        with self.assertRaisesRegex(shared.ConfigError, "Unknown top-level field"):
            shared.motor_config_from_dict(unknown)

        missing = copy.deepcopy(self.raw)
        del missing["drive"]
        with self.assertRaisesRegex(shared.ConfigError, "Missing top-level field"):
            shared.motor_config_from_dict(missing)

    def test_unknown_and_missing_section_fields_are_rejected(self) -> None:
        unknown = copy.deepcopy(self.raw)
        unknown["machine"]["unexpected"] = 1
        with self.assertRaisesRegex(shared.ConfigError, "Unknown machine field"):
            shared.motor_config_from_dict(unknown)

        missing = copy.deepcopy(self.raw)
        del missing["electromagnetic"]["lq_h"]
        with self.assertRaisesRegex(shared.ConfigError, "Missing electromagnetic field"):
            shared.motor_config_from_dict(missing)

    def test_json_booleans_are_not_accepted_as_integer_fields(self) -> None:
        invalid = copy.deepcopy(self.raw)
        invalid["machine"]["pole_pairs"] = True
        with self.assertRaisesRegex(shared.ConfigError, "must be an integer"):
            shared.motor_config_from_dict(invalid)


class SharedPhysicalConfigurationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = shared.load_motor_config()

    def test_geometry_fingerprint_matches_validated_electromagnetic_source(self) -> None:
        fingerprint = shared.physical_model_fingerprint(
            self.config.machine,
            self.config.materials,
            self.config.conventions,
        )
        self.assertEqual(
            fingerprint,
            self.config.electromagnetic.source_physical_model_sha256,
        )
        self.assertEqual(
            fingerprint,
            "cae1e23d408cbddc550a4c5cca0629572ad61aa277933f937a2d2c386ade03e0",
        )

    def test_geometry_change_invalidates_electromagnetic_parameters(self) -> None:
        changed_machine = replace(
            self.config.machine,
            airgap_m=self.config.machine.airgap_m * 1.01,
        )
        changed_fingerprint = shared.physical_model_fingerprint(
            changed_machine,
            self.config.materials,
            self.config.conventions,
        )
        self.assertNotEqual(
            changed_fingerprint,
            self.config.electromagnetic.source_physical_model_sha256,
        )
        with self.assertRaisesRegex(shared.ConfigError, "geometry changed"):
            replace(self.config, machine=changed_machine).validate()

    def test_validated_machine_is_18_slot_20_pole_and_75_turns_per_phase(self) -> None:
        machine = self.config.machine
        self.assertEqual(machine.slots, 18)
        self.assertEqual(2 * machine.pole_pairs, 20)
        self.assertEqual(machine.turns_per_slot, 25)
        self.assertEqual(machine.series_coils_per_phase, 3)
        self.assertEqual(machine.turns_per_phase, 75)
        self.assertAlmostEqual(machine.fundamental_winding_factor, 0.9452136366)


class QuickModelPhysicsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = shared.load_motor_config()

    def assert_inside_constraints(self, point: quick.OperatingPoint) -> None:
        self.assertLessEqual(
            point.phase_current_peak_a,
            self.config.drive.current_peak_limit_a * (1.0 + 1e-9),
        )
        self.assertLessEqual(
            point.phase_voltage_peak_v,
            quick.phase_voltage_limit(self.config) * (1.0 + 1e-9),
        )

    def test_five_amp_torque_and_20c_copper_loss_use_peak_dq_convention(self) -> None:
        torque = quick.electromagnetic_torque_nm(self.config, 0.0, 5.0)
        copper_loss_20c = (
            1.5
            * self.config.electromagnetic.phase_resistance_20c_ohm
            * 5.0**2
        )
        self.assertAlmostEqual(float(torque), 2.4011183685, places=9)
        self.assertAlmostEqual(copper_loss_20c, 8.25, places=12)

        point = quick.evaluate_operating_point(self.config, 0.0, float(torque))
        self.assertAlmostEqual(point.id_peak_a, 0.0, places=12)
        self.assertAlmostEqual(point.iq_peak_a, 5.0, places=9)
        self.assertAlmostEqual(
            point.copper_loss_w,
            1.5 * point.phase_resistance_ohm * 5.0**2,
            places=9,
        )

    def test_low_speed_minimum_current_solution_has_negligible_d_axis_current(self) -> None:
        point = quick.evaluate_operating_point(self.config, 100.0, 5.0)
        self.assertTrue(point.feasible)
        id_grid_step = self.config.drive.current_peak_limit_a / (
            self.config.sweep.current_search_points - 1
        )
        self.assertLessEqual(abs(point.id_peak_a), id_grid_step)
        self.assertGreater(point.iq_peak_a, 0.0)
        self.assert_inside_constraints(point)

    def test_1000_rpm_uses_field_weakening_inside_current_and_voltage_limits(self) -> None:
        point = quick.evaluate_operating_point(self.config, 1000.0, 3.0)
        self.assertTrue(point.feasible)
        self.assertLess(point.id_peak_a, -0.1)
        self.assertGreater(point.iq_peak_a, 0.0)
        self.assert_inside_constraints(point)

    def test_infeasible_torque_is_clipped_to_a_constrained_point(self) -> None:
        point = quick.evaluate_operating_point(self.config, 0.0, 100.0)
        self.assertFalse(point.feasible)
        self.assertLess(point.shaft_torque_nm, point.requested_shaft_torque_nm)
        self.assertTrue(point.limit_reason.startswith("clipped:"))
        self.assertIn("current", point.limit_reason)
        self.assert_inside_constraints(point)

    def test_every_feasible_envelope_point_obeys_both_constraints(self) -> None:
        sweep = replace(
            self.config.sweep,
            speed_points=13,
            current_search_points=1001,
        )
        config = replace(self.config, sweep=sweep)
        result = quick.run_quick_analysis(config)
        self.assertEqual(len(result.envelope_points), 13)
        self.assertTrue(any(point.feasible for point in result.envelope_points))
        for point in result.envelope_points:
            with self.subTest(speed_rpm=point.speed_rpm):
                self.assertLessEqual(
                    point.phase_current_peak_a,
                    config.drive.current_peak_limit_a * (1.0 + 1e-9),
                )
                if point.feasible:
                    self.assertLessEqual(
                        point.phase_voltage_peak_v,
                        quick.phase_voltage_limit(config) * (1.0 + 1e-9),
                    )
                else:
                    self.assertTrue(point.limit_reason.startswith("clipped:"))

    def test_dq_and_system_power_balances_close(self) -> None:
        for speed_rpm, torque_nm in ((100.0, 5.0), (1000.0, 3.0)):
            with self.subTest(speed_rpm=speed_rpm):
                point = quick.evaluate_operating_point(
                    self.config, speed_rpm, torque_nm
                )
                scale = max(1.0, abs(point.dc_input_power_w))
                self.assertLessEqual(abs(point.dq_power_balance_error_w), 1e-10 * scale)
                self.assertLessEqual(
                    abs(point.system_power_balance_error_w), 1e-10 * scale
                )

        calibrated = replace(
            self.config,
            drive=replace(self.config.drive, inverter_efficiency=0.95),
            losses=replace(
                self.config.losses,
                calibrated=True,
                core_loss_ref_w=8.0,
                magnet_loss_ref_w=2.0,
                mechanical_loss_ref_w=1.0,
            ),
        )
        point = quick.evaluate_operating_point(calibrated, 500.0, 3.0)
        self.assertGreater(point.core_loss_w, 0.0)
        self.assertGreater(point.magnet_loss_w, 0.0)
        self.assertGreater(point.mechanical_loss_w, 0.0)
        self.assertGreater(point.inverter_loss_w, 0.0)
        scale = max(1.0, abs(point.dc_input_power_w))
        self.assertLessEqual(abs(point.dq_power_balance_error_w), 1e-10 * scale)
        self.assertLessEqual(abs(point.system_power_balance_error_w), 1e-10 * scale)

    def test_currents_above_femm_validation_point_are_marked_extrapolated(self) -> None:
        torque_per_amp = float(quick.electromagnetic_torque_nm(self.config, 0.0, 1.0))
        inside = quick.evaluate_operating_point(
            self.config, 0.0, 4.0 * torque_per_amp
        )
        outside = quick.evaluate_operating_point(
            self.config, 0.0, 6.0 * torque_per_amp
        )
        self.assertAlmostEqual(inside.phase_current_peak_a, 4.0, places=8)
        self.assertAlmostEqual(outside.phase_current_peak_a, 6.0, places=8)
        self.assertFalse(inside.extrapolated)
        self.assertTrue(outside.extrapolated)


class FemmConfigurationMappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # Importing the adapter is read-only and does not open FEMM/COM.
        _install_optional_femm_import_stubs()
        import femm_spm_template as femm_model

        cls.femm_model = femm_model
        cls.config = shared.load_motor_config()

    def test_machine_mapping_converts_si_to_millimetres_without_drift(self) -> None:
        machine = self.femm_model.machine_from_motor_config(self.config)
        source = self.config.machine
        self.assertEqual(machine.pole_pairs, 10)
        self.assertEqual(machine.slots, 18)
        self.assertEqual(machine.turns_per_slot, 25)
        self.assertAlmostEqual(machine.stack_length_mm, source.stack_length_m * 1000.0)
        self.assertAlmostEqual(machine.r_shaft_mm, source.shaft_radius_m * 1000.0)
        self.assertAlmostEqual(machine.r_rotor_mm, source.rotor_radius_m * 1000.0)
        self.assertAlmostEqual(machine.mag_thickness_mm, source.magnet_thickness_m * 1000.0)
        self.assertAlmostEqual(machine.airgap_mm, source.airgap_m * 1000.0)
        self.assertAlmostEqual(machine.r_stator_inner_mm, 32.0)
        self.assertAlmostEqual(machine.r_stator_outer_mm, 60.0)
        self.assertAlmostEqual(machine.r_air_outer_mm, 70.0)
        self.assertAlmostEqual(machine.r_slot_outer_mm, 40.0)
        self.assertEqual(machine.steel_material_name, self.config.materials.steel_library_name)
        self.assertEqual(machine.copper_material_name, self.config.materials.copper_library_name)

    def test_loss_thermal_and_drive_mapping_preserves_shared_values(self) -> None:
        loss = self.femm_model.loss_from_motor_config(self.config)
        drive = self.femm_model.drive_from_motor_config(self.config)
        self.assertEqual(loss.r_phase_20, self.config.electromagnetic.phase_resistance_20c_ohm)
        self.assertEqual(loss.alpha_cu, self.config.electromagnetic.copper_temp_coeff_per_k)
        self.assertEqual(loss.calibrated, self.config.losses.calibrated)
        self.assertEqual(loss.core_loss_ref_w, self.config.losses.core_loss_ref_w)
        self.assertEqual(loss.magnet_loss_ref_w, self.config.losses.magnet_loss_ref_w)
        self.assertEqual(
            loss.winding_to_case_rth_k_per_w,
            self.config.thermal.winding_to_case_rth_k_per_w,
        )
        self.assertEqual(loss.t_amb, self.config.thermal.ambient_c)
        self.assertEqual(drive.v_dc, self.config.drive.dc_bus_v)
        self.assertEqual(drive.i_max, self.config.drive.current_peak_limit_a)
        self.assertEqual(drive.modulation_limit, self.config.drive.modulation_limit)


if __name__ == "__main__":
    unittest.main()
