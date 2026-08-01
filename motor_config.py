"""Shared, validated configuration for the quick model and FEMM model.

All dimensional values in this file and in the JSON configuration use SI units.
The FEMM adapter is the only place where metres are converted to millimetres.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
from pathlib import Path
from typing import Any, TypeVar


class ConfigError(ValueError):
    """Raised when a motor configuration is incomplete or inconsistent."""


def _require_int(name: str, value: Any) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"{name} must be an integer, not {type(value).__name__}.")


def _require_number(name: str, value: Any) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{name} must be a finite number, not {type(value).__name__}.")
    if not math.isfinite(float(value)):
        raise ConfigError(f"{name} must be finite.")


def _strict_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ConfigError(f"Duplicate JSON object key is not allowed: {key}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ConfigError(f"Non-finite JSON number is not allowed: {value}")


@dataclass(frozen=True)
class MachineConfig:
    pole_pairs: int = 10
    slots: int = 18
    stack_length_m: float = 0.080
    shaft_radius_m: float = 0.020
    rotor_radius_m: float = 0.030
    magnet_thickness_m: float = 0.0015
    airgap_m: float = 0.0005
    stator_outer_radius_m: float = 0.060
    outer_air_radius_m: float = 0.070
    magnet_arc_ratio: float = 0.85
    turns_per_slot: int = 25
    series_coils_per_phase: int = 3
    fundamental_winding_factor: float = 0.9452136366
    winding_layout: str = "18s20p_single_layer_v1"
    connection: str = "wye"
    slot_depth_m: float = 0.008

    @property
    def turns_per_phase(self) -> int:
        return self.turns_per_slot * self.series_coils_per_phase

    def validate(self) -> None:
        for name in ("pole_pairs", "slots", "turns_per_slot", "series_coils_per_phase"):
            _require_int(f"machine.{name}", getattr(self, name))
        _require_number("machine.magnet_arc_ratio", self.magnet_arc_ratio)
        _require_number(
            "machine.fundamental_winding_factor", self.fundamental_winding_factor
        )
        if not isinstance(self.winding_layout, str) or not isinstance(self.connection, str):
            raise ConfigError("machine winding_layout and connection must be strings.")
        if (
            self.pole_pairs <= 0
            or self.slots <= 0
            or self.turns_per_slot <= 0
            or self.series_coils_per_phase <= 0
        ):
            raise ConfigError("pole_pairs, slots, and turns_per_slot must be positive.")
        if (self.slots, self.pole_pairs, self.winding_layout) != (
            18,
            10,
            "18s20p_single_layer_v1",
        ):
            raise ConfigError(
                "This version supports only the validated 18-slot/20-pole single-layer winding."
            )
        if self.series_coils_per_phase != 3 or self.connection != "wye":
            raise ConfigError(
                "The validated winding requires three series coils per phase and a wye connection."
            )
        if not 0.0 < self.fundamental_winding_factor <= 1.0:
            raise ConfigError("machine.fundamental_winding_factor must be in (0, 1].")
        if not 0.0 < self.magnet_arc_ratio < 1.0:
            raise ConfigError("machine.magnet_arc_ratio must be between 0 and 1.")
        positive = {
            "stack_length_m": self.stack_length_m,
            "shaft_radius_m": self.shaft_radius_m,
            "rotor_radius_m": self.rotor_radius_m,
            "magnet_thickness_m": self.magnet_thickness_m,
            "airgap_m": self.airgap_m,
            "stator_outer_radius_m": self.stator_outer_radius_m,
            "outer_air_radius_m": self.outer_air_radius_m,
            "slot_depth_m": self.slot_depth_m,
        }
        for name, value in positive.items():
            _require_number(f"machine.{name}", value)
            if value <= 0.0:
                raise ConfigError(f"machine.{name} must be positive.")
        magnet_outer = self.rotor_radius_m + self.magnet_thickness_m
        stator_inner = magnet_outer + self.airgap_m
        slot_outer = stator_inner + self.slot_depth_m
        radii = (
            self.shaft_radius_m,
            self.rotor_radius_m,
            magnet_outer,
            stator_inner,
            slot_outer,
            self.stator_outer_radius_m,
            self.outer_air_radius_m,
        )
        if any(left >= right for left, right in zip(radii, radii[1:])):
            raise ConfigError(
                "Machine radii must increase in the order shaft, rotor, magnet, "
                "air gap, slot, stator, outer air."
            )


@dataclass(frozen=True)
class ElectromagneticConfig:
    phase_resistance_20c_ohm: float = 0.220
    copper_temp_coeff_per_k: float = 0.00393
    ld_h: float = 0.0009660890764148228
    lq_h: float = 0.0009420675310741955
    psi_pm_wb: float = 0.03201491158003713
    validated_current_peak_a: float = 5.0
    parameter_source: str = "FEMM electromagnetic acceptance, medium mesh"
    source_physical_model_sha256: str = "cae1e23d408cbddc550a4c5cca0629572ad61aa277933f937a2d2c386ade03e0"

    def validate(self) -> None:
        positive = {
            "phase_resistance_20c_ohm": self.phase_resistance_20c_ohm,
            "ld_h": self.ld_h,
            "lq_h": self.lq_h,
            "psi_pm_wb": self.psi_pm_wb,
            "validated_current_peak_a": self.validated_current_peak_a,
        }
        for name, value in positive.items():
            _require_number(f"electromagnetic.{name}", value)
            if value <= 0.0:
                raise ConfigError(f"electromagnetic.{name} must be positive.")
        _require_number(
            "electromagnetic.copper_temp_coeff_per_k",
            self.copper_temp_coeff_per_k,
        )
        if not isinstance(self.parameter_source, str) or not self.parameter_source.strip():
            raise ConfigError("electromagnetic.parameter_source must be a non-empty string.")
        if (
            not isinstance(self.source_physical_model_sha256, str)
            or len(self.source_physical_model_sha256) != 64
        ):
            raise ConfigError(
                "electromagnetic.source_physical_model_sha256 must be a 64-character hash."
            )
        if not 0.0 < self.copper_temp_coeff_per_k < 0.02:
            raise ConfigError(
                "electromagnetic.copper_temp_coeff_per_k must be between 0 and 0.02."
            )


@dataclass(frozen=True)
class MagneticMaterialConfig:
    magnet_br_t: float = 1.22
    magnet_relative_permeability: float = 1.05
    steel_library_name: str = "M-19 Steel"
    copper_library_name: str = "Copper"

    def validate(self) -> None:
        _require_number("materials.magnet_br_t", self.magnet_br_t)
        _require_number(
            "materials.magnet_relative_permeability",
            self.magnet_relative_permeability,
        )
        if self.magnet_br_t <= 0.0 or self.magnet_relative_permeability <= 0.0:
            raise ConfigError("materials magnet Br and relative permeability must be positive.")
        if not isinstance(self.steel_library_name, str) or not isinstance(
            self.copper_library_name, str
        ):
            raise ConfigError("materials library names must be strings.")
        if not self.steel_library_name.strip() or not self.copper_library_name.strip():
            raise ConfigError("materials library names cannot be empty.")


@dataclass(frozen=True)
class ConventionConfig:
    dq_current_basis: str = "peak"
    phase_sequence: str = "A-B-C"
    electrical_zero_offset_deg: float = -100.0
    electrical_angle_sign: int = -1
    positive_iq_torque_sign: int = 1

    def validate(self) -> None:
        _require_int("conventions.electrical_angle_sign", self.electrical_angle_sign)
        _require_int(
            "conventions.positive_iq_torque_sign", self.positive_iq_torque_sign
        )
        _require_number(
            "conventions.electrical_zero_offset_deg",
            self.electrical_zero_offset_deg,
        )
        if not isinstance(self.dq_current_basis, str) or not isinstance(
            self.phase_sequence, str
        ):
            raise ConfigError("convention names must be strings.")
        if self.dq_current_basis != "peak":
            raise ConfigError("Only peak-valued dq currents are supported.")
        if self.phase_sequence != "A-B-C":
            raise ConfigError("Only the validated A-B-C phase sequence is supported.")
        if self.electrical_angle_sign != -1:
            raise ConfigError("FEMM geometry currently requires electrical_angle_sign = -1.")
        if self.positive_iq_torque_sign != 1:
            raise ConfigError("The implemented torque convention requires positive Iq torque.")


def physical_model_fingerprint(
    machine: MachineConfig,
    materials: MagneticMaterialConfig,
    conventions: ConventionConfig,
) -> str:
    payload = {
        "model_variant": "femm_18s20p_single_layer_v2",
        "machine": asdict(machine),
        "materials": asdict(materials),
        "conventions": asdict(conventions),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(canonical).hexdigest()


@dataclass(frozen=True)
class DriveConfig:
    dc_bus_v: float = 48.0
    modulation_limit: float = 0.95
    current_peak_limit_a: float = 15.0
    inverter_efficiency: float = 1.0

    def validate(self) -> None:
        for name in (
            "dc_bus_v",
            "modulation_limit",
            "current_peak_limit_a",
            "inverter_efficiency",
        ):
            _require_number(f"drive.{name}", getattr(self, name))
        if self.dc_bus_v <= 0.0 or self.current_peak_limit_a <= 0.0:
            raise ConfigError("drive.dc_bus_v and current_peak_limit_a must be positive.")
        if not 0.0 < self.modulation_limit <= 1.0:
            raise ConfigError("drive.modulation_limit must be in (0, 1].")
        if not 0.0 < self.inverter_efficiency <= 1.0:
            raise ConfigError("drive.inverter_efficiency must be in (0, 1].")


@dataclass(frozen=True)
class LossConfig:
    calibrated: bool = False
    core_loss_ref_w: float = 0.0
    core_loss_ref_electrical_hz: float = 100.0
    core_frequency_exponent: float = 1.5
    core_flux_density_ref_t: float = 1.2
    magnet_loss_ref_w: float = 0.0
    magnet_loss_ref_electrical_hz: float = 100.0
    magnet_frequency_exponent: float = 2.0
    mechanical_loss_ref_w: float = 0.0
    mechanical_loss_ref_rpm: float = 1000.0
    mechanical_speed_exponent: float = 2.0

    def validate(self) -> None:
        if not isinstance(self.calibrated, bool):
            raise ConfigError("losses.calibrated must be a JSON boolean.")
        nonnegative = {
            "core_loss_ref_w": self.core_loss_ref_w,
            "magnet_loss_ref_w": self.magnet_loss_ref_w,
            "mechanical_loss_ref_w": self.mechanical_loss_ref_w,
        }
        for name, value in nonnegative.items():
            _require_number(f"losses.{name}", value)
            if value < 0.0:
                raise ConfigError(f"losses.{name} must be non-negative.")
        positive = {
            "core_loss_ref_electrical_hz": self.core_loss_ref_electrical_hz,
            "core_frequency_exponent": self.core_frequency_exponent,
            "core_flux_density_ref_t": self.core_flux_density_ref_t,
            "magnet_loss_ref_electrical_hz": self.magnet_loss_ref_electrical_hz,
            "magnet_frequency_exponent": self.magnet_frequency_exponent,
            "mechanical_loss_ref_rpm": self.mechanical_loss_ref_rpm,
            "mechanical_speed_exponent": self.mechanical_speed_exponent,
        }
        for name, value in positive.items():
            _require_number(f"losses.{name}", value)
            if value <= 0.0:
                raise ConfigError(f"losses.{name} must be positive.")


@dataclass(frozen=True)
class ThermalConfig:
    ambient_c: float = 25.0
    winding_to_case_rth_k_per_w: float = 0.35
    magnet_to_case_rth_k_per_w: float = 0.25
    case_to_ambient_rth_k_per_w: float = 0.15
    winding_temperature_limit_c: float = 155.0
    magnet_temperature_limit_c: float = 100.0

    def validate(self) -> None:
        _require_number("thermal.ambient_c", self.ambient_c)
        _require_number(
            "thermal.winding_temperature_limit_c",
            self.winding_temperature_limit_c,
        )
        _require_number(
            "thermal.magnet_temperature_limit_c",
            self.magnet_temperature_limit_c,
        )
        resistances = {
            "winding_to_case_rth_k_per_w": self.winding_to_case_rth_k_per_w,
            "magnet_to_case_rth_k_per_w": self.magnet_to_case_rth_k_per_w,
            "case_to_ambient_rth_k_per_w": self.case_to_ambient_rth_k_per_w,
        }
        for name, value in resistances.items():
            _require_number(f"thermal.{name}", value)
            if value <= 0.0:
                raise ConfigError(f"thermal.{name} must be positive.")
        if self.winding_temperature_limit_c <= self.ambient_c:
            raise ConfigError("thermal.winding_temperature_limit_c must exceed ambient_c.")
        if self.magnet_temperature_limit_c <= self.ambient_c:
            raise ConfigError("thermal.magnet_temperature_limit_c must exceed ambient_c.")


@dataclass(frozen=True)
class SweepConfig:
    speed_min_rpm: float = 0.0
    speed_max_rpm: float = 1500.0
    speed_points: int = 61
    requested_torque_nm: float = 2.0
    current_search_points: int = 2001

    def validate(self) -> None:
        _require_int("sweep.speed_points", self.speed_points)
        _require_int("sweep.current_search_points", self.current_search_points)
        for name in ("speed_min_rpm", "speed_max_rpm", "requested_torque_nm"):
            _require_number(f"sweep.{name}", getattr(self, name))
        if self.speed_min_rpm < 0.0 or self.speed_max_rpm <= self.speed_min_rpm:
            raise ConfigError("sweep speed range must satisfy 0 <= min < max.")
        if self.speed_points < 3:
            raise ConfigError("sweep.speed_points must be at least 3.")
        if self.requested_torque_nm < 0.0:
            raise ConfigError("sweep.requested_torque_nm must be non-negative.")
        if self.current_search_points < 101:
            raise ConfigError("sweep.current_search_points must be at least 101.")


@dataclass(frozen=True)
class MotorConfig:
    schema_version: int = 1
    name: str = "18s20p_spm_demo"
    machine: MachineConfig = MachineConfig()
    electromagnetic: ElectromagneticConfig = ElectromagneticConfig()
    materials: MagneticMaterialConfig = MagneticMaterialConfig()
    conventions: ConventionConfig = ConventionConfig()
    drive: DriveConfig = DriveConfig()
    losses: LossConfig = LossConfig()
    thermal: ThermalConfig = ThermalConfig()
    sweep: SweepConfig = SweepConfig()

    def validate(self) -> None:
        _require_int("schema_version", self.schema_version)
        if self.schema_version != 1:
            raise ConfigError(f"Unsupported schema_version {self.schema_version}; expected 1.")
        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigError("Configuration name cannot be empty.")
        self.machine.validate()
        self.electromagnetic.validate()
        self.materials.validate()
        self.conventions.validate()
        expected_fingerprint = physical_model_fingerprint(
            self.machine, self.materials, self.conventions
        )
        if self.electromagnetic.source_physical_model_sha256 != expected_fingerprint:
            raise ConfigError(
                "The machine/material/convention geometry changed without revalidating "
                "psi_pm/Ld/Lq; source_physical_model_sha256 does not match."
            )
        self.drive.validate()
        self.losses.validate()
        self.thermal.validate()
        self.sweep.validate()


T = TypeVar("T")


def _construct_section(section_type: type[T], raw: Any, section_name: str) -> T:
    if raw is None:
        raise ConfigError(f"Missing required section: {section_name}")
    if not isinstance(raw, dict):
        raise ConfigError(f"{section_name} must be a JSON object.")
    allowed = {field.name for field in fields(section_type)}
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ConfigError(f"Unknown {section_name} field(s): {', '.join(unknown)}")
    missing = sorted(allowed - set(raw))
    if missing:
        raise ConfigError(f"Missing {section_name} field(s): {', '.join(missing)}")
    try:
        return section_type(**raw)
    except TypeError as exc:
        raise ConfigError(f"Invalid {section_name} section: {exc}") from exc


def motor_config_from_dict(raw: dict[str, Any]) -> MotorConfig:
    if not isinstance(raw, dict):
        raise ConfigError("The configuration root must be a JSON object.")
    section_types = {
        "machine": MachineConfig,
        "electromagnetic": ElectromagneticConfig,
        "materials": MagneticMaterialConfig,
        "conventions": ConventionConfig,
        "drive": DriveConfig,
        "losses": LossConfig,
        "thermal": ThermalConfig,
        "sweep": SweepConfig,
    }
    allowed = {"schema_version", "name", *section_types}
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ConfigError(f"Unknown top-level field(s): {', '.join(unknown)}")
    missing = sorted(allowed - set(raw))
    if missing:
        raise ConfigError(f"Missing top-level field(s): {', '.join(missing)}")
    kwargs: dict[str, Any] = {
        "schema_version": raw["schema_version"],
        "name": raw["name"],
    }
    for name, section_type in section_types.items():
        kwargs[name] = _construct_section(section_type, raw[name], name)
    config = MotorConfig(**kwargs)
    config.validate()
    return config


def default_config_path() -> Path:
    return Path(__file__).resolve().with_name("motor_config.json")


def load_motor_config(path: str | Path | None = None) -> MotorConfig:
    config_path = default_config_path() if path is None else Path(path)
    try:
        raw = json.loads(
            config_path.read_text(encoding="utf-8-sig"),
            object_pairs_hook=_strict_json_object,
            parse_constant=_reject_json_constant,
        )
    except FileNotFoundError as exc:
        raise ConfigError(f"Configuration file not found: {config_path}") from exc
    except json.JSONDecodeError as exc:
        raise ConfigError(
            f"Invalid JSON in {config_path} at line {exc.lineno}, column {exc.colno}: {exc.msg}"
        ) from exc
    return motor_config_from_dict(raw)


def config_to_dict(config: MotorConfig) -> dict[str, Any]:
    config.validate()
    return asdict(config)


def canonical_config_sha256(config: MotorConfig) -> str:
    payload = json.dumps(
        config_to_dict(config), sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_motor_config(config: MotorConfig, path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(config_to_dict(config), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return destination
