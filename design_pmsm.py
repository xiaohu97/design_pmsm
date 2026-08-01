"""Fast, constraint-aware operating-point model for the shared 18-slot/20-pole SPM.

This is a steady-state screening model.  It uses FEMM-derived ``psi_pm``, ``Ld``
and ``Lq`` from :mod:`motor_config`, solves the current circle and voltage
ellipse at every speed, and keeps uncalibrated loss estimates visibly separate
from the electromagnetic power balance.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from motor_config import (
    ConfigError,
    MotorConfig,
    config_to_dict,
    default_config_path,
    load_motor_config,
)


@dataclass(frozen=True)
class CurrentSolution:
    id_peak_a: float
    iq_peak_a: float
    feasible: bool


@dataclass(frozen=True)
class OperatingPoint:
    speed_rpm: float
    requested_shaft_torque_nm: float
    electromagnetic_torque_nm: float
    shaft_torque_nm: float
    feasible: bool
    extrapolated: bool
    loss_model_calibrated: bool
    current_limited: bool
    voltage_limited: bool
    field_weakening: bool
    thermal_limit_exceeded: bool
    limit_reason: str
    id_peak_a: float
    iq_peak_a: float
    phase_current_peak_a: float
    phase_current_rms_a: float
    phase_current_peak_limit_a: float
    current_utilization: float
    vd_peak_v: float
    vq_peak_v: float
    phase_voltage_peak_v: float
    phase_voltage_peak_limit_v: float
    voltage_utilization: float
    phase_resistance_ohm: float
    electrical_frequency_hz: float
    dq_input_power_w: float
    airgap_power_w: float
    shaft_power_w: float
    copper_loss_w: float
    core_loss_w: float
    magnet_loss_w: float
    mechanical_loss_w: float
    inverter_loss_w: float
    total_system_loss_w: float
    dc_input_power_w: float
    efficiency_estimate: float
    dq_power_balance_error_w: float
    system_power_balance_error_w: float
    case_temperature_c: float
    winding_temperature_c: float
    magnet_temperature_c: float


@dataclass(frozen=True)
class AnalysisResult:
    config: MotorConfig
    load_points: tuple[OperatingPoint, ...]
    envelope_points: tuple[OperatingPoint, ...]


def phase_voltage_limit(config: MotorConfig) -> float:
    """Available fundamental phase-voltage peak under the configured convention."""
    return config.drive.modulation_limit * config.drive.dc_bus_v / math.sqrt(3.0)


def electromagnetic_torque_nm(
    config: MotorConfig,
    id_peak_a: float | np.ndarray,
    iq_peak_a: float | np.ndarray,
) -> float | np.ndarray:
    em = config.electromagnetic
    return 1.5 * config.machine.pole_pairs * (
        em.psi_pm_wb * iq_peak_a
        + (em.ld_h - em.lq_h) * id_peak_a * iq_peak_a
    )


def dq_voltage(
    config: MotorConfig,
    speed_rpm: float,
    id_peak_a: float | np.ndarray,
    iq_peak_a: float | np.ndarray,
    phase_resistance_ohm: float,
) -> tuple[float | np.ndarray, float | np.ndarray]:
    em = config.electromagnetic
    omega_m = speed_rpm * 2.0 * math.pi / 60.0
    omega_e = config.machine.pole_pairs * omega_m
    vd = phase_resistance_ohm * id_peak_a - omega_e * em.lq_h * iq_peak_a
    vq = phase_resistance_ohm * iq_peak_a + omega_e * (
        em.ld_h * id_peak_a + em.psi_pm_wb
    )
    return vd, vq


def _id_search_grid(config: MotorConfig) -> np.ndarray:
    return np.linspace(
        -config.drive.current_peak_limit_a,
        0.0,
        config.sweep.current_search_points,
        dtype=float,
    )


def solve_maximum_torque_currents(
    config: MotorConfig,
    speed_rpm: float,
    phase_resistance_ohm: float,
) -> CurrentSolution:
    """Find the maximum positive torque inside the current circle and voltage ellipse."""
    id_values = _id_search_grid(config)
    i_limit = config.drive.current_peak_limit_a
    q_current_max = np.sqrt(np.maximum(0.0, i_limit**2 - id_values**2))

    em = config.electromagnetic
    omega_m = speed_rpm * 2.0 * math.pi / 60.0
    omega_e = config.machine.pole_pairs * omega_m
    rs = phase_resistance_ohm
    v_limit = phase_voltage_limit(config)

    # |v_dq|^2 <= Vmax^2 is quadratic in iq for each selected id.
    a = (omega_e * em.lq_h) ** 2 + rs**2
    b = 2.0 * (
        -rs * id_values * omega_e * em.lq_h
        + rs * omega_e * (em.ld_h * id_values + em.psi_pm_wb)
    )
    c = (
        (rs * id_values) ** 2
        + (omega_e * (em.ld_h * id_values + em.psi_pm_wb)) ** 2
        - v_limit**2
    )
    discriminant = b**2 - 4.0 * a * c
    has_voltage_interval = discriminant >= 0.0
    root = np.sqrt(np.maximum(discriminant, 0.0))
    q_voltage_low = (-b - root) / (2.0 * a)
    q_voltage_high = (-b + root) / (2.0 * a)
    iq_values = np.minimum(q_current_max, q_voltage_high)
    valid = (
        has_voltage_interval
        & (iq_values >= np.maximum(0.0, q_voltage_low) - 1e-10)
        & (iq_values >= 0.0)
    )
    torques = np.asarray(
        electromagnetic_torque_nm(config, id_values, iq_values), dtype=float
    )
    torques[~valid] = -np.inf
    if not np.any(valid):
        # The speed is beyond even the full-current field-weakening capability.
        vd, vq = dq_voltage(config, speed_rpm, id_values, 0.0, rs)
        voltage = np.hypot(vd, vq)
        index = int(np.argmin(voltage))
        return CurrentSolution(float(id_values[index]), 0.0, False)
    index = int(np.argmax(torques))
    return CurrentSolution(float(id_values[index]), float(iq_values[index]), True)


def solve_requested_torque_currents(
    config: MotorConfig,
    speed_rpm: float,
    electromagnetic_torque_request_nm: float,
    phase_resistance_ohm: float,
) -> CurrentSolution:
    """Find the minimum-current feasible point, or clip to the maximum-torque point."""
    id_values = _id_search_grid(config)
    em = config.electromagnetic
    torque_per_iq = 1.5 * config.machine.pole_pairs * (
        em.psi_pm_wb + (em.ld_h - em.lq_h) * id_values
    )
    iq_values = np.divide(
        electromagnetic_torque_request_nm,
        torque_per_iq,
        out=np.full_like(id_values, np.inf),
        where=np.abs(torque_per_iq) > 1e-12,
    )
    current = np.hypot(id_values, iq_values)
    vd, vq = dq_voltage(
        config, speed_rpm, id_values, iq_values, phase_resistance_ohm
    )
    voltage = np.hypot(vd, vq)
    valid = (
        np.isfinite(iq_values)
        & (iq_values >= -1e-12)
        & (current <= config.drive.current_peak_limit_a * (1.0 + 1e-10))
        & (voltage <= phase_voltage_limit(config) * (1.0 + 1e-10))
    )
    if np.any(valid):
        score = np.where(valid, current, np.inf)
        index = int(np.argmin(score))
        return CurrentSolution(float(id_values[index]), float(iq_values[index]), True)
    clipped = solve_maximum_torque_currents(config, speed_rpm, phase_resistance_ohm)
    return CurrentSolution(clipped.id_peak_a, clipped.iq_peak_a, False)


def reference_losses(config: MotorConfig, speed_rpm: float) -> tuple[float, float, float]:
    """Return explicit reference-law core, magnet, and mechanical losses.

    These are trend estimates until ``losses.calibrated`` is set after fitting
    the reference values to measurements or a dedicated transient-loss model.
    """
    loss = config.losses
    electrical_hz = config.machine.pole_pairs * abs(speed_rpm) / 60.0
    core = loss.core_loss_ref_w * (
        electrical_hz / loss.core_loss_ref_electrical_hz
    ) ** loss.core_frequency_exponent
    magnet = loss.magnet_loss_ref_w * (
        electrical_hz / loss.magnet_loss_ref_electrical_hz
    ) ** loss.magnet_frequency_exponent
    mechanical = loss.mechanical_loss_ref_w * (
        abs(speed_rpm) / loss.mechanical_loss_ref_rpm
    ) ** loss.mechanical_speed_exponent
    return core, magnet, mechanical


def _temperatures(
    config: MotorConfig,
    copper_loss_w: float,
    core_loss_w: float,
    magnet_loss_w: float,
    mechanical_loss_w: float,
) -> tuple[float, float, float]:
    thermal = config.thermal
    motor_loss = copper_loss_w + core_loss_w + magnet_loss_w + mechanical_loss_w
    case_c = thermal.ambient_c + motor_loss * thermal.case_to_ambient_rth_k_per_w
    winding_c = case_c + copper_loss_w * thermal.winding_to_case_rth_k_per_w
    magnet_c = case_c + magnet_loss_w * thermal.magnet_to_case_rth_k_per_w
    return case_c, winding_c, magnet_c


def _limit_reason(feasible: bool, current_util: float, voltage_util: float) -> str:
    active: list[str] = []
    if current_util >= 0.995:
        active.append("current")
    if voltage_util >= 0.995:
        active.append("voltage")
    if not active:
        active.append("speed_voltage" if not feasible else "none")
    joined = "+".join(active)
    return f"clipped:{joined}" if not feasible else joined


def evaluate_operating_point(
    config: MotorConfig,
    speed_rpm: float,
    requested_shaft_torque_nm: float,
    *,
    maximum_torque: bool = False,
) -> OperatingPoint:
    """Solve one hot-resistance steady operating point."""
    if speed_rpm < 0.0:
        raise ValueError("This simple motor-mode sweep requires non-negative speed.")
    if requested_shaft_torque_nm < 0.0:
        raise ValueError("This simple motor-mode sweep requires non-negative torque.")

    em = config.electromagnetic
    omega_m = speed_rpm * 2.0 * math.pi / 60.0
    core_loss, magnet_loss, mechanical_loss = reference_losses(config, speed_rpm)
    mechanical_torque = mechanical_loss / omega_m if omega_m > 1e-12 else 0.0
    electromagnetic_request = requested_shaft_torque_nm + mechanical_torque

    winding_c = config.thermal.ambient_c
    solution = CurrentSolution(0.0, 0.0, False)
    case_c = magnet_c = config.thermal.ambient_c
    rs = em.phase_resistance_20c_ohm
    for _ in range(50):
        rs = em.phase_resistance_20c_ohm * (
            1.0 + em.copper_temp_coeff_per_k * (winding_c - 20.0)
        )
        if maximum_torque:
            solution = solve_maximum_torque_currents(config, speed_rpm, rs)
        else:
            solution = solve_requested_torque_currents(
                config, speed_rpm, electromagnetic_request, rs
            )
        copper_loss = 1.5 * rs * (
            solution.id_peak_a**2 + solution.iq_peak_a**2
        )
        new_case_c, new_winding_c, new_magnet_c = _temperatures(
            config, copper_loss, core_loss, magnet_loss, mechanical_loss
        )
        if abs(new_winding_c - winding_c) < 1e-7:
            case_c, winding_c, magnet_c = new_case_c, new_winding_c, new_magnet_c
            break
        # Damping makes the fixed point robust when users enter a high Rth.
        winding_c = 0.5 * winding_c + 0.5 * new_winding_c
        case_c, magnet_c = new_case_c, new_magnet_c
    else:
        raise RuntimeError("Hot-resistance/thermal fixed point did not converge.")

    id_a, iq_a = solution.id_peak_a, solution.iq_peak_a
    vd, vq = dq_voltage(config, speed_rpm, id_a, iq_a, rs)
    voltage = math.hypot(float(vd), float(vq))
    current = math.hypot(id_a, iq_a)
    torque_em = float(electromagnetic_torque_nm(config, id_a, iq_a))
    shaft_torque = torque_em - mechanical_torque
    if abs(shaft_torque) < 1e-12:
        shaft_torque = 0.0
    airgap_power = torque_em * omega_m
    # Keep the sign for the balance.  At a forced speed beyond the feasible
    # motor envelope, a negative value correctly denotes external shaft input.
    shaft_power = shaft_torque * omega_m
    copper_loss = 1.5 * rs * current**2
    dq_input = 1.5 * (float(vd) * id_a + float(vq) * iq_a)
    dq_balance = dq_input - airgap_power - copper_loss

    motor_ac_input = dq_input + core_loss + magnet_loss
    inverter_loss = (
        max(0.0, motor_ac_input) * (1.0 / config.drive.inverter_efficiency - 1.0)
    )
    dc_input = motor_ac_input + inverter_loss
    total_loss = (
        copper_loss + core_loss + magnet_loss + mechanical_loss + inverter_loss
    )
    system_balance = dc_input - shaft_power - total_loss
    efficiency = (
        shaft_power / dc_input
        if shaft_power >= 0.0 and dc_input > 1e-12
        else 0.0
    )
    current_util = current / config.drive.current_peak_limit_a
    voltage_util = voltage / phase_voltage_limit(config)

    requested_feasible = solution.feasible
    if maximum_torque:
        requested_feasible = solution.feasible
        requested_shaft_torque_nm = max(0.0, shaft_torque)
    reason = _limit_reason(requested_feasible, current_util, voltage_util)
    return OperatingPoint(
        speed_rpm=float(speed_rpm),
        requested_shaft_torque_nm=float(requested_shaft_torque_nm),
        electromagnetic_torque_nm=torque_em,
        shaft_torque_nm=shaft_torque,
        feasible=bool(requested_feasible),
        extrapolated=current > em.validated_current_peak_a * (1.0 + 1e-9),
        loss_model_calibrated=config.losses.calibrated,
        current_limited=current_util >= 0.995,
        voltage_limited=voltage_util >= 0.995,
        field_weakening=id_a < -1e-6,
        thermal_limit_exceeded=(
            winding_c > config.thermal.winding_temperature_limit_c
            or magnet_c > config.thermal.magnet_temperature_limit_c
        ),
        limit_reason=reason,
        id_peak_a=id_a,
        iq_peak_a=iq_a,
        phase_current_peak_a=current,
        phase_current_rms_a=current / math.sqrt(2.0),
        phase_current_peak_limit_a=config.drive.current_peak_limit_a,
        current_utilization=current_util,
        vd_peak_v=float(vd),
        vq_peak_v=float(vq),
        phase_voltage_peak_v=voltage,
        phase_voltage_peak_limit_v=phase_voltage_limit(config),
        voltage_utilization=voltage_util,
        phase_resistance_ohm=rs,
        electrical_frequency_hz=config.machine.pole_pairs * speed_rpm / 60.0,
        dq_input_power_w=dq_input,
        airgap_power_w=airgap_power,
        shaft_power_w=shaft_power,
        copper_loss_w=copper_loss,
        core_loss_w=core_loss,
        magnet_loss_w=magnet_loss,
        mechanical_loss_w=mechanical_loss,
        inverter_loss_w=inverter_loss,
        total_system_loss_w=total_loss,
        dc_input_power_w=dc_input,
        efficiency_estimate=efficiency,
        dq_power_balance_error_w=dq_balance,
        system_power_balance_error_w=system_balance,
        case_temperature_c=case_c,
        winding_temperature_c=winding_c,
        magnet_temperature_c=magnet_c,
    )


def run_quick_analysis(config: MotorConfig) -> AnalysisResult:
    config.validate()
    speeds = np.linspace(
        config.sweep.speed_min_rpm,
        config.sweep.speed_max_rpm,
        config.sweep.speed_points,
    )
    load_points: list[OperatingPoint] = []
    envelope_points: list[OperatingPoint] = []
    for speed in speeds:
        envelope_points.append(
            evaluate_operating_point(config, float(speed), 0.0, maximum_torque=True)
        )
        load_points.append(
            evaluate_operating_point(
                config, float(speed), config.sweep.requested_torque_nm
            )
        )
    return AnalysisResult(config, tuple(load_points), tuple(envelope_points))


def _resolved_config_hash(config: MotorConfig) -> str:
    payload = json.dumps(
        config_to_dict(config), sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_summary(result: AnalysisResult) -> dict[str, Any]:
    load = result.load_points
    envelope = result.envelope_points
    low_speed_torque = max(point.shaft_torque_nm for point in envelope[:3])
    base_candidates = [
        point.speed_rpm
        for point in envelope
        if point.shaft_torque_nm < 0.99 * low_speed_torque
    ]
    feasible_speeds = [point.speed_rpm for point in load if point.feasible]
    warnings: list[str] = []
    if not result.config.losses.calibrated:
        warnings.append(
            "Additional loss references are uncalibrated (zero in the supplied config); "
            "efficiency is a copper-only upper bound and thermal results are not design sign-off data."
        )
    if any(point.extrapolated for point in (*load, *envelope)):
        warnings.append(
            "Some currents exceed the "
            f"{result.config.electromagnetic.validated_current_peak_a:g} A FEMM "
            "validation point; those rows are marked extrapolated."
        )
    warnings.append(
        "Temperatures are steady thermal-resistance estimates, not a transient drive-cycle result."
    )
    if any(point.thermal_limit_exceeded for point in load):
        warnings.append("At least one requested-load row exceeds a configured thermal limit.")
    warnings.append(
        "The quick model uses fixed FEMM-derived psi_pm/Ld/Lq and does not model saturation or demagnetization."
    )
    return {
        "model": "constraint-aware steady dq screening model",
        "config_name": result.config.name,
        "resolved_config_sha256": _resolved_config_hash(result.config),
        "parameter_source": result.config.electromagnetic.parameter_source,
        "turns_per_phase": result.config.machine.turns_per_phase,
        "fundamental_winding_factor": (
            result.config.machine.fundamental_winding_factor
        ),
        "torque_constant_nm_per_a_peak": (
            1.5
            * result.config.machine.pole_pairs
            * result.config.electromagnetic.psi_pm_wb
        ),
        "current_convention": "dq and phase peak; phase RMS = peak/sqrt(2)",
        "phase_sequence": result.config.conventions.phase_sequence,
        "electrical_zero_offset_deg": result.config.conventions.electrical_zero_offset_deg,
        "electrical_angle_sign": result.config.conventions.electrical_angle_sign,
        "phase_voltage_peak_limit_v": phase_voltage_limit(result.config),
        "phase_current_peak_limit_a": result.config.drive.current_peak_limit_a,
        "electromagnetic_validated_current_peak_a": (
            result.config.electromagnetic.validated_current_peak_a
        ),
        "requested_shaft_torque_nm": result.config.sweep.requested_torque_nm,
        "efficiency_interpretation": (
            "calibrated estimate"
            if result.config.losses.calibrated
            else "copper-only upper bound"
        ),
        "no_load_id0_voltage_base_speed_rpm": (
            phase_voltage_limit(result.config)
            / (
                result.config.machine.pole_pairs
                * result.config.electromagnetic.psi_pm_wb
            )
            * 60.0
            / (2.0 * math.pi)
        ),
        "estimated_constant_torque_base_speed_rpm": (
            base_candidates[0] if base_candidates else None
        ),
        "low_speed_max_shaft_torque_nm": low_speed_torque,
        "peak_envelope_shaft_power_w": max(point.shaft_power_w for point in envelope),
        "highest_feasible_requested_torque_speed_rpm": (
            max(feasible_speeds) if feasible_speeds else None
        ),
        "max_sweep_winding_temperature_c": max(
            point.winding_temperature_c for point in load
        ),
        "max_feasible_requested_winding_temperature_c": max(
            (point.winding_temperature_c for point in load if point.feasible),
            default=None,
        ),
        "max_sweep_magnet_temperature_c": max(
            point.magnet_temperature_c for point in load
        ),
        "max_feasible_requested_magnet_temperature_c": max(
            (point.magnet_temperature_c for point in load if point.feasible),
            default=None,
        ),
        "max_abs_dq_power_balance_error_w": max(
            abs(point.dq_power_balance_error_w) for point in (*load, *envelope)
        ),
        "max_abs_system_power_balance_error_w": max(
            abs(point.system_power_balance_error_w) for point in (*load, *envelope)
        ),
        "all_requested_points_feasible": all(point.feasible for point in load),
        "any_load_thermal_limit_exceeded": any(
            point.thermal_limit_exceeded for point in load
        ),
        "contains_electromagnetic_extrapolation": any(
            point.extrapolated for point in (*load, *envelope)
        ),
        "warnings": warnings,
        "resolved_config": config_to_dict(result.config),
    }


def save_csv(result: AnalysisResult, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "quick_sweep.csv"
    load_fields = tuple(asdict(result.load_points[0]))
    envelope_fields = tuple(asdict(result.envelope_points[0]))
    fieldnames = [f"load_{name}" for name in load_fields] + [
        f"envelope_{name}" for name in envelope_fields
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for load, envelope in zip(result.load_points, result.envelope_points):
            row = {f"load_{key}": value for key, value in asdict(load).items()}
            row.update(
                {f"envelope_{key}": value for key, value in asdict(envelope).items()}
            )
            writer.writerow(row)
    return path


def save_summary(result: AnalysisResult, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "quick_summary.json"
    path.write_text(
        json.dumps(build_summary(result), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def save_resolved_config(result: AnalysisResult, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "resolved_motor_config.json"
    path.write_text(
        json.dumps(config_to_dict(result.config), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def plot_results(
    result: AnalysisResult,
    out_dir: Path,
    *,
    save_png: bool = True,
    show: bool = False,
) -> Path | None:
    if not show:
        import matplotlib

        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    speed = np.array([point.speed_rpm for point in result.load_points])
    load_torque = np.array([point.shaft_torque_nm for point in result.load_points])
    max_torque = np.array([point.shaft_torque_nm for point in result.envelope_points])
    load_power = np.array([point.shaft_power_w for point in result.load_points])
    max_power = np.array([point.shaft_power_w for point in result.envelope_points])
    efficiency = np.array(
        [point.efficiency_estimate for point in result.load_points]
    ) * 100.0
    id_values = np.array([point.id_peak_a for point in result.load_points])
    iq_values = np.array([point.iq_peak_a for point in result.load_points])
    voltage_util = np.array(
        [point.voltage_utilization for point in result.load_points]
    ) * 100.0
    current_util = np.array(
        [point.current_utilization for point in result.load_points]
    ) * 100.0
    winding_temp = np.array(
        [point.winding_temperature_c for point in result.load_points]
    )
    feasible = np.array([point.feasible for point in result.load_points], dtype=bool)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    ax = axes[0, 0]
    ax.plot(speed, max_torque, label="Max feasible shaft torque")
    ax.plot(speed, load_torque, label="Achieved requested torque")
    ax.axhline(
        result.config.sweep.requested_torque_nm,
        color="gray",
        linestyle="--",
        label="Requested torque",
    )
    if np.any(~feasible):
        ax.scatter(speed[~feasible], load_torque[~feasible], color="red", s=18, label="Clipped")
    ax.set_ylabel("Torque [N m]")
    ax.grid(True)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(speed, max_power, label="Envelope shaft power")
    ax.plot(speed, load_power, label="Load shaft power")
    ax.set_ylabel("Power [W]")
    ax.grid(True)
    ax_eff = ax.twinx()
    efficiency_label = (
        "Estimated load efficiency"
        if result.config.losses.calibrated
        else "Copper-only efficiency upper bound"
    )
    ax_eff.plot(
        speed,
        efficiency,
        color="tab:green",
        linestyle="--",
        label=efficiency_label,
    )
    ax_eff.set_ylabel("Efficiency [%]", color="tab:green")
    handles, labels = ax.get_legend_handles_labels()
    h2, l2 = ax_eff.get_legend_handles_labels()
    ax.legend(handles + h2, labels + l2, fontsize=8)

    ax = axes[1, 0]
    ax.plot(speed, id_values, label="Id peak")
    ax.plot(speed, iq_values, label="Iq peak")
    ax.set_xlabel("Speed [rpm]")
    ax.set_ylabel("Current [A peak]")
    ax.grid(True)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(speed, current_util, label="Current utilization")
    ax.plot(speed, voltage_util, label="Voltage utilization")
    ax.axhline(100.0, color="gray", linestyle="--")
    ax.set_xlabel("Speed [rpm]")
    ax.set_ylabel("Utilization [%]")
    ax.grid(True)
    ax_temp = ax.twinx()
    ax_temp.plot(speed, winding_temp, color="tab:red", linestyle=":", label="Winding temperature")
    ax_temp.set_ylabel("Steady winding temp [degC]", color="tab:red")
    handles, labels = ax.get_legend_handles_labels()
    h2, l2 = ax_temp.get_legend_handles_labels()
    ax.legend(handles + h2, labels + l2, fontsize=8)

    fig.suptitle(f"{result.config.name}: quick constraint-aware screening")
    fig.tight_layout()
    png_path: Path | None = None
    if save_png:
        out_dir.mkdir(parents=True, exist_ok=True)
        png_path = out_dir / "quick_sweep.png"
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return png_path


def _apply_cli_overrides(config: MotorConfig, args: argparse.Namespace) -> MotorConfig:
    drive = config.drive
    sweep = config.sweep
    if args.vdc is not None:
        drive = replace(drive, dc_bus_v=args.vdc)
    if args.current_limit is not None:
        drive = replace(drive, current_peak_limit_a=args.current_limit)
    if args.torque is not None:
        sweep = replace(sweep, requested_torque_nm=args.torque)
    if args.max_rpm is not None:
        sweep = replace(sweep, speed_max_rpm=args.max_rpm)
    if args.points is not None:
        sweep = replace(sweep, speed_points=args.points)
    resolved = replace(config, drive=drive, sweep=sweep)
    resolved.validate()
    return resolved


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fast 18-slot/20-pole PMSM torque-speed screening with current and voltage constraints."
        )
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="Shared SI-unit JSON configuration.",
    )
    parser.add_argument("--out", type=Path, default=Path("output_quick"))
    parser.add_argument("--torque", type=float, help="Override requested shaft torque [N m].")
    parser.add_argument("--max-rpm", type=float, help="Override maximum sweep speed [rpm].")
    parser.add_argument("--points", type=int, help="Override number of sweep speeds.")
    parser.add_argument("--vdc", type=float, help="Override DC bus voltage [V].")
    parser.add_argument(
        "--current-limit", type=float, help="Override dq/phase peak current limit [A]."
    )
    parser.add_argument("--show", action="store_true", help="Show the result figure.")
    parser.add_argument("--no-png", action="store_true", help="Do not save quick_sweep.png.")
    # Compatibility with the original launcher; CSV and PNG are now saved by default.
    parser.add_argument("--save-csv", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--save-png", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--no-show", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--dt", type=float, help=argparse.SUPPRESS)
    parser.add_argument("--t-end", type=float, help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        config = _apply_cli_overrides(load_motor_config(args.config), args)
        result = run_quick_analysis(config)
    except (ConfigError, ValueError, RuntimeError) as exc:
        raise SystemExit(f"Configuration/model error: {exc}") from exc

    if args.dt is not None or args.t_end is not None:
        print("Note: --dt/--t-end are ignored; this version performs a steady-state sweep.")
    csv_path = save_csv(result, args.out)
    summary_path = save_summary(result, args.out)
    config_path = save_resolved_config(result, args.out)
    show_plot = args.show and not args.no_show
    png_path = None
    if not args.no_png or show_plot:
        png_path = plot_results(
            result,
            args.out,
            save_png=not args.no_png,
            show=show_plot,
        )
    summary = build_summary(result)
    print(f"Motor: {config.name} ({config.machine.slots} slots / {2 * config.machine.pole_pairs} poles)")
    print(
        "Limits: "
        f"Vphase_peak={summary['phase_voltage_peak_limit_v']:.3f} V, "
        f"Ipeak={summary['phase_current_peak_limit_a']:.3f} A"
    )
    print(
        "Envelope: "
        f"Tmax_low={summary['low_speed_max_shaft_torque_nm']:.3f} N m, "
        "constant_torque_base_speed~"
        f"{summary['estimated_constant_torque_base_speed_rpm']} rpm, "
        f"Pmax={summary['peak_envelope_shaft_power_w']:.1f} W"
    )
    print(f"Saved CSV: {csv_path.resolve()}")
    print(f"Saved summary: {summary_path.resolve()}")
    print(f"Saved resolved config: {config_path.resolve()}")
    if png_path is not None:
        print(f"Saved figure: {png_path.resolve()}")
    if not config.losses.calibrated:
        print("WARNING: loss/thermal reference values are not calibrated; see quick_summary.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
