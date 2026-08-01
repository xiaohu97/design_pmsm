import math
import os
from pathlib import Path
import sys
from contextlib import contextmanager
from tempfile import TemporaryDirectory
import types
import unittest
from unittest import mock

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")

# These tests exercise geometry and acceptance math only.  Keep them runnable
# in a lightweight unit-test environment without plotting, progress, or COM.
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

import femm_spm_template as motor


class WindingTests(unittest.TestCase):
    def setUp(self):
        self.machine = motor.Machine()

    def test_18s20p_winding_is_closed_balanced_and_120_degree(self):
        layout = motor.winding_layout(self.machine)
        self.assertEqual(len(layout), 18)
        for phase in motor.PHASES:
            phase_sides = [sign for slot_phase, sign in layout if slot_phase == phase]
            self.assertEqual(len(phase_sides), 6)
            self.assertEqual(sum(phase_sides), 0)

        metrics = motor.validate_winding_layout(self.machine)
        self.assertAlmostEqual(metrics["phase_axis_a_deg"], -10.0, places=10)
        self.assertAlmostEqual(metrics["phase_axis_b_deg"], -130.0, places=10)
        self.assertAlmostEqual(metrics["phase_axis_c_deg"], 110.0, places=10)
        self.assertAlmostEqual(
            metrics["fundamental_winding_factor"], 0.9452136366, places=10
        )

    def test_each_single_layer_coil_uses_adjacent_slots(self):
        layout = motor.winding_layout(self.machine)
        for slot in range(0, self.machine.slots, 2):
            phase_1, sign_1 = layout[slot]
            phase_2, sign_2 = layout[slot + 1]
            self.assertEqual(phase_1, phase_2)
            self.assertEqual(sign_1, -sign_2)

    def test_electrical_zero_is_relative_to_phase_a_axis(self):
        offset_deg = math.degrees(motor.electrical_zero_offset_rad(self.machine))
        self.assertAlmostEqual(offset_deg, -100.0, places=10)
        self.assertAlmostEqual(
            math.degrees(motor.electrical_angle_rad(self.machine, 1.0)),
            -110.0,
            places=10,
        )

    def test_dq_abc_round_trip_and_norm(self):
        for theta in np.linspace(-2.0 * math.pi, 2.0 * math.pi, 17):
            id_a, iq_a = 2.5, -7.0
            abc = motor.dq_to_abc(id_a, iq_a, float(theta))
            recovered = motor.abc_to_dq(*abc, float(theta))
            self.assertAlmostEqual(sum(abc), 0.0, places=12)
            self.assertAlmostEqual(recovered[0], id_a, places=12)
            self.assertAlmostEqual(recovered[1], iq_a, places=12)
            self.assertAlmostEqual(
                sum(value * value for value in abc),
                1.5 * (id_a * id_a + iq_a * iq_a),
                places=12,
            )

        q_axis_abc = motor.dq_to_abc(0.0, 1.0, 0.0)
        self.assertAlmostEqual(q_axis_abc[0], 0.0, places=12)
        self.assertLess(q_axis_abc[1], 0.0)
        self.assertGreater(q_axis_abc[2], 0.0)

    def test_magnet_edges_are_unique_and_cover_full_circle(self):
        edges = motor.magnet_edge_angles(self.machine)
        self.assertEqual(len(edges), 40)
        sweeps = [
            (edges[(index + 1) % len(edges)] - edge) % 360.0
            for index, edge in enumerate(edges)
        ]
        self.assertTrue(all(sweep > 0.0 for sweep in sweeps))
        self.assertAlmostEqual(sum(sweeps), 360.0, places=10)
        self.assertEqual(
            {round(value, 10) for value in sweeps},
            {2.7, 15.3},
        )

    def test_airgap_mesh_levels_hold_other_regions_fixed(self):
        levels = motor.AIRGAP_MESH_LEVELS
        self.assertEqual(
            [mesh.airgap_radial_layers for mesh in levels],
            [3, 5, 8, 12],
        )
        self.assertEqual(
            [mesh.airgap_arc_deg for mesh in levels],
            [1.0, 0.5, 0.25, 0.15],
        )
        fixed_fields = (
            "magnet_size_mm",
            "rotor_size_mm",
            "stator_size_mm",
            "coil_size_mm",
            "outer_air_size_mm",
        )
        for field in fixed_fields:
            self.assertEqual(len({getattr(mesh, field) for mesh in levels}), 1)

    def test_convergence_checkpoint_key_covers_mesh_current_and_points(self):
        cache_dir = Path("cache")
        base = motor.torque_sweep_checkpoint_path(
            cache_dir, self.machine, motor.AIRGAP_MESH_LEVELS[0], 5.0, 6
        )
        changed_mesh = motor.torque_sweep_checkpoint_path(
            cache_dir, self.machine, motor.AIRGAP_MESH_LEVELS[1], 5.0, 6
        )
        changed_current = motor.torque_sweep_checkpoint_path(
            cache_dir, self.machine, motor.AIRGAP_MESH_LEVELS[0], 6.0, 6
        )
        changed_points = motor.torque_sweep_checkpoint_path(
            cache_dir, self.machine, motor.AIRGAP_MESH_LEVELS[0], 5.0, 12
        )
        self.assertEqual(len({base, changed_mesh, changed_current, changed_points}), 4)

    def test_case_value_tokens_preserve_fraction_and_sign(self):
        self.assertEqual(motor.case_value_token(600.0), "600")
        self.assertEqual(motor.case_value_token(2.5), "2p5")
        self.assertEqual(motor.case_value_token(-2.5), "m2p5")
        self.assertEqual(
            len(
                {
                    motor.case_value_token(value)
                    for value in (2.0, 2.5, -2.0, -2.5)
                }
            ),
            4,
        )


class FieldSnapshotTests(unittest.TestCase):
    def test_snapshot_uses_explicit_angle_and_records_reproducible_metadata(self):
        machine = motor.Machine()
        mesh = motor.MeshConfig(name="medium")
        expected = motor.OperatingPoint(2.4, 0.032, 0.0047)

        @contextmanager
        def fake_model(*_args, **_kwargs):
            yield

        with TemporaryDirectory() as temporary:
            out_dir = Path(temporary)
            with (
                mock.patch.object(
                    motor, "temporary_femm_model", side_effect=fake_model
                ),
                mock.patch.object(motor, "rotate_rotor") as rotate,
                mock.patch.object(
                    motor, "solve_operating_point", return_value=expected
                ) as solve,
                mock.patch.object(motor, "plot_field_maps") as field_plot,
                mock.patch.object(
                    motor, "plot_airgap_flux_density"
                ) as airgap_plot,
            ):
                actual = motor.run_field_snapshot(
                    machine,
                    out_dir,
                    rpm=600.0,
                    iq_peak_a=2.5,
                    rotor_angle_deg=-1.25,
                    analyses={"field", "airgap"},
                    mesh=mesh,
                    field_radial_points=24,
                    field_angular_points=120,
                    airgap_points=360,
                )

            self.assertEqual(actual, expected)
            rotate.assert_called_once_with(-1.25)
            solve.assert_called_once_with(machine, -1.25, 0.0, 2.5)
            field_plot.assert_called_once_with(
                machine,
                out_dir,
                "rpm600_iq2p5_am1p25",
                False,
                radial_points=24,
                angular_points=120,
            )
            airgap_plot.assert_called_once_with(
                machine,
                out_dir,
                "rpm600_iq2p5_am1p25",
                False,
                num_points=360,
            )

            metadata = (
                out_dir / "field_snapshot_rpm600_iq2p5_am1p25.csv"
            ).read_text(encoding="utf-8")
            self.assertIn("rpm_context,rotor_angle_deg", metadata)
            self.assertIn("600.0,-1.25,0.0,2.5,2.4,0.032,0.0047,medium", metadata)


class AcceptanceMathTests(unittest.TestCase):
    def test_maxwell_airgap_torque_uses_periodic_stress_integral(self):
        machine = motor.Machine()
        radius_mm = machine.r_mag_outer_mm + 0.5 * machine.airgap_mm
        br = np.full(360, 0.8)
        bt = np.full(360, 0.1)
        expected = (
            machine.stack_length_mm
            * 1e-3
            * (radius_mm * 1e-3) ** 2
            * 2.0
            * math.pi
            * 0.08
            / motor.MU0
        )
        self.assertAlmostEqual(
            motor.maxwell_airgap_torque_from_samples(machine, radius_mm, br, bt),
            expected,
            places=12,
        )

    def test_synthetic_surface_pmsm_passes_acceptance(self):
        machine = motor.Machine()
        psi_pm = 0.032
        iq_test = 5.0
        torque_slope = 1.5 * machine.pole_pairs * psi_pm
        samples = []
        for index, angle in enumerate(np.linspace(0.0, 36.0, 12, endpoint=False)):
            cogging = 0.01 * math.sin(2.0 * math.pi * index / 12.0)
            samples.append(
                motor.ElectromagneticSample(
                    rotor_angle_deg=float(angle),
                    theta_e_deg=10.0 + machine.pole_pairs * float(angle),
                    psi_d0_wb=psi_pm,
                    psi_q0_wb=1e-4 * math.sin(2.0 * math.pi * index / 12.0),
                    cogging_torque_nm=cogging,
                    torque_q_pos_nm=cogging + torque_slope * iq_test,
                    torque_q_neg_nm=cogging - torque_slope * iq_test,
                    ld_h=0.00100,
                    lq_h=0.00105,
                )
            )

        metrics = motor.evaluate_electromagnetic_acceptance(
            machine, samples, iq_test
        )
        gated = [
            metric for metric in metrics.values() if metric["criterion"] != "report"
        ]
        self.assertTrue(all(metric["passed"] for metric in gated))
        self.assertAlmostEqual(
            metrics["torque_slope_nm_per_a"]["value"], torque_slope, places=12
        )

    def test_torque_sweep_subtracts_cogging_at_same_angle(self):
        rows = np.array([
            [0.0, -0.02, 2.98],
            [0.5, 0.01, 3.01],
            [1.0, 0.02, 3.02],
            [1.5, -0.01, 2.99],
        ])
        summary = motor.summarize_torque_sweep(rows)
        self.assertAlmostEqual(summary["loaded_mean_nm"], 3.0, places=12)
        self.assertAlmostEqual(summary["loaded_ripple_pp_nm"], 0.0, places=12)
        self.assertAlmostEqual(summary["cogging_mean_nm"], 0.0, places=12)
        self.assertAlmostEqual(summary["cogging_pp_nm"], 0.04, places=12)


if __name__ == "__main__":
    unittest.main()
