"""FEMM 2-D SPMSM model and electromagnetic acceptance workflow.

Outputs supported (selected via ``--analysis``):
  basic        torque / loss waveforms and steady thermal-resistance estimate
  field        flux-density cloud and flux-line maps
  airgap       air-gap Bn / Bt distribution
  cogging      zero-current cogging-torque curve
  inductance   centered-incremental Ld / Lq versus rotor position
  validate     dq alignment, torque, and incremental-inductance acceptance
  meshcheck    one-angle electromagnetic check for global-mesh sensitivity
  convergence  isolated gap-mesh, rotor-angle, and torque-integral convergence
  emap         efficiency map (speed by current sweep)
  tncurve      torque-speed envelope
  all          run every analysis above
"""
from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
import re
import sys
import tempfile
import uuid
import multiprocessing
import traceback
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

from motor_config import (
    ConfigError,
    MotorConfig as SharedMotorConfig,
    canonical_config_sha256,
    config_to_dict,
    default_config_path,
    load_motor_config,
    physical_model_fingerprint,
)

try:
    import femm
except ImportError as exc:
    raise SystemExit(
        "pyfemm is not installed. Run: pip install pyfemm, "
        "and install FEMM 4.2 first."
    ) from exc

# ── constants & group IDs ──────────────────────────────────────────
MU0 = 4e-7 * math.pi
ROTOR_GROUP = 2
STATOR_GROUP = 3
COIL_GROUP = 4
AIRGAP_GROUP = 5
PHASES = ("A", "B", "C")

# Validated single-layer 18-slot / 20-pole winding. Each adjacent pair is one
# concentrated coil; every phase has three positive and three negative sides.
WINDING_18S20P = (
    ("A", +1), ("A", -1),
    ("B", -1), ("B", +1),
    ("C", +1), ("C", -1),
    ("C", +1), ("C", -1),
    ("A", -1), ("A", +1),
    ("B", +1), ("B", -1),
    ("B", +1), ("B", -1),
    ("C", -1), ("C", +1),
    ("A", +1), ("A", -1),
)

def get_tmp_fem() -> str:
    """Returns a unique temporary filename for FEMM model to avoid clashes in multiprocessing."""
    return os.path.join(tempfile.gettempdir(), f"femm_spm_tmp_{os.getpid()}_{uuid.uuid4().hex[:6]}.fem")

# ── data classes ───────────────────────────────────────────────────

@dataclass
class Machine:
    pole_pairs: int = 10 
    slots: int = 18 
    stack_length_mm: float = 80.0 
    r_shaft_mm: float = 25.0  # Legacy name: Air-filled bore radius.
    r_rotor_mm: float = 30.0 
    mag_thickness_mm: float = 1.5 
    airgap_mm: float = 0.5 
    r_stator_outer_mm: float = 45.0
    r_air_outer_mm: float = 70.0 
    magnet_arc_ratio: float = 0.85 
    turns_per_slot: int = 25 
    slot_depth_mm: float = 8.0 
    magnet_br_t: float = 1.22
    magnet_relative_permeability: float = 1.05
    steel_material_name: str = "M-19 Steel"
    copper_material_name: str = "Copper"

    def __post_init__(self) -> None:
        if self.pole_pairs <= 0 or self.slots <= 0:
            raise ValueError("pole_pairs and slots must be positive.")
        if not (0.0 < self.magnet_arc_ratio < 1.0):
            raise ValueError("magnet_arc_ratio must be between 0 and 1.")
        if self.magnet_br_t <= 0.0 or self.magnet_relative_permeability <= 0.0:
            raise ValueError("Magnet Br and relative permeability must be positive.")
        radii = (
            self.r_shaft_mm,
            self.r_rotor_mm,
            self.r_mag_outer_mm,
            self.r_stator_inner_mm,
            self.r_slot_outer_mm,
            self.r_stator_outer_mm,
            self.r_air_outer_mm,
        )
        if any(radius <= 0.0 for radius in radii):
            raise ValueError("All radii must be positive.")
        if any(a >= b for a, b in zip(radii, radii[1:])):
            raise ValueError("Machine radii must be strictly increasing.")

    @property
    def r_mag_outer_mm(self):
        return self.r_rotor_mm + self.mag_thickness_mm

    @property
    def r_stator_inner_mm(self):
        return self.r_mag_outer_mm + self.airgap_mm

    @property
    def r_slot_outer_mm(self):
        return self.r_stator_inner_mm + self.slot_depth_mm


@dataclass
class LossThermal:
    r_phase_20: float = 0.22         
    alpha_cu: float = 0.00393        
    calibrated: bool = False
    core_loss_ref_w: float = 0.0
    core_loss_ref_electrical_hz: float = 100.0
    core_frequency_exponent: float = 1.5
    core_flux_density_ref_t: float = 1.2
    magnet_loss_ref_w: float = 0.0
    magnet_loss_ref_electrical_hz: float = 100.0
    magnet_frequency_exponent: float = 2.0
    winding_to_case_rth_k_per_w: float = 0.35
    magnet_to_case_rth_k_per_w: float = 0.25
    case_to_ambient_rth_k_per_w: float = 0.15
    t_amb: float = 25.0

    def __post_init__(self) -> None:
        if self.r_phase_20 <= 0.0 or self.alpha_cu <= 0.0:
            raise ValueError("Copper resistance and temperature coefficient must be positive.")
        if min(
            self.core_loss_ref_w,
            self.magnet_loss_ref_w,
        ) < 0.0:
            raise ValueError("Reference losses must be non-negative.")
        positive = (
            self.core_loss_ref_electrical_hz,
            self.core_frequency_exponent,
            self.core_flux_density_ref_t,
            self.magnet_loss_ref_electrical_hz,
            self.magnet_frequency_exponent,
            self.winding_to_case_rth_k_per_w,
            self.magnet_to_case_rth_k_per_w,
            self.case_to_ambient_rth_k_per_w,
        )
        if any(value <= 0.0 for value in positive):
            raise ValueError("Loss reference scales and thermal resistances must be positive.")


@dataclass
class Drive:
    """Inverter / voltage‑source parameters for T‑N envelope."""
    v_dc: float = 48.0             
    i_max: float = 15.0
    modulation_limit: float = 0.95

    def __post_init__(self) -> None:
        if self.v_dc <= 0.0 or self.i_max <= 0.0:
            raise ValueError("DC voltage and peak-current limit must be positive.")
        if not 0.0 < self.modulation_limit <= 1.0:
            raise ValueError("modulation_limit must be in (0, 1].")


def machine_from_motor_config(config: SharedMotorConfig) -> Machine:
    """Explicit SI-to-mm adapter; FEMM currently supports only 18s/20p."""
    config.validate()
    source = config.machine
    if (source.slots, source.pole_pairs, source.winding_layout) != (
        18,
        10,
        "18s20p_single_layer_v1",
    ):
        raise ConfigError("FEMM supports only the validated 18s/20p winding layout.")
    material = config.materials
    machine = Machine(
        pole_pairs=source.pole_pairs,
        slots=source.slots,
        stack_length_mm=source.stack_length_m * 1000.0,
        r_shaft_mm=source.shaft_radius_m * 1000.0,
        r_rotor_mm=source.rotor_radius_m * 1000.0,
        mag_thickness_mm=source.magnet_thickness_m * 1000.0,
        airgap_mm=source.airgap_m * 1000.0,
        r_stator_outer_mm=source.stator_outer_radius_m * 1000.0,
        r_air_outer_mm=source.outer_air_radius_m * 1000.0,
        magnet_arc_ratio=source.magnet_arc_ratio,
        turns_per_slot=source.turns_per_slot,
        slot_depth_mm=source.slot_depth_m * 1000.0,
        magnet_br_t=material.magnet_br_t,
        magnet_relative_permeability=material.magnet_relative_permeability,
        steel_material_name=material.steel_library_name,
        copper_material_name=material.copper_library_name,
    )
    offset_deg = math.degrees(electrical_zero_offset_rad(machine))
    if not math.isclose(
        offset_deg,
        config.conventions.electrical_zero_offset_deg,
        abs_tol=1e-9,
    ):
        raise ConfigError(
            "Shared Park zero does not match the winding-derived FEMM zero: "
            f"{config.conventions.electrical_zero_offset_deg} vs {offset_deg} deg."
        )
    return machine


def loss_from_motor_config(config: SharedMotorConfig) -> LossThermal:
    em, loss, thermal = config.electromagnetic, config.losses, config.thermal
    return LossThermal(
        r_phase_20=em.phase_resistance_20c_ohm,
        alpha_cu=em.copper_temp_coeff_per_k,
        calibrated=loss.calibrated,
        core_loss_ref_w=loss.core_loss_ref_w,
        core_loss_ref_electrical_hz=loss.core_loss_ref_electrical_hz,
        core_frequency_exponent=loss.core_frequency_exponent,
        core_flux_density_ref_t=loss.core_flux_density_ref_t,
        magnet_loss_ref_w=loss.magnet_loss_ref_w,
        magnet_loss_ref_electrical_hz=loss.magnet_loss_ref_electrical_hz,
        magnet_frequency_exponent=loss.magnet_frequency_exponent,
        winding_to_case_rth_k_per_w=thermal.winding_to_case_rth_k_per_w,
        magnet_to_case_rth_k_per_w=thermal.magnet_to_case_rth_k_per_w,
        case_to_ambient_rth_k_per_w=thermal.case_to_ambient_rth_k_per_w,
        t_amb=thermal.ambient_c,
    )


def drive_from_motor_config(config: SharedMotorConfig) -> Drive:
    return Drive(
        v_dc=config.drive.dc_bus_v,
        i_max=config.drive.current_peak_limit_a,
        modulation_limit=config.drive.modulation_limit,
    )


@dataclass
class SweepCfg:
    points_per_electrical_cycle: int = 72

    def __post_init__(self) -> None:
        if self.points_per_electrical_cycle < 3:
            raise ValueError("points_per_electrical_cycle must be at least 3.")


@dataclass(frozen=True)
class MeshConfig:
    """Explicit mesh controls; sizes are maximum triangle edge lengths."""

    name: str = "medium"
    airgap_radial_layers: int = 5
    airgap_arc_deg: float = 0.5
    magnet_size_mm: float = 0.35
    rotor_size_mm: float = 0.8
    stator_size_mm: float = 1.0
    coil_size_mm: float = 0.7
    outer_air_size_mm: float = 2.0

    def __post_init__(self) -> None:
        if self.airgap_radial_layers < 3:
            raise ValueError("airgap_radial_layers must be at least 3.")
        values = (
            self.airgap_arc_deg,
            self.magnet_size_mm,
            self.rotor_size_mm,
            self.stator_size_mm,
            self.coil_size_mm,
            self.outer_air_size_mm,
        )
        if any(value <= 0.0 for value in values):
            raise ValueError("All mesh controls must be positive.")

    def airgap_size_mm(self, machine: Machine) -> float:
        return machine.airgap_mm / self.airgap_radial_layers


MESH_LEVELS = (
    MeshConfig(
        name="coarse", airgap_radial_layers=3, airgap_arc_deg=1.0,
        magnet_size_mm=0.6, rotor_size_mm=1.4, stator_size_mm=1.8,
        coil_size_mm=1.2, outer_air_size_mm=3.0,
    ),
    MeshConfig(),
    MeshConfig(
        name="fine", airgap_radial_layers=8, airgap_arc_deg=0.25,
        magnet_size_mm=0.22, rotor_size_mm=0.5, stator_size_mm=0.65,
        coil_size_mm=0.45, outer_air_size_mm=1.2,
    ),
)

# Isolated air-gap refinements for convergence: every non-gap size is held at
# the normal medium value so that the observed change can be attributed to the
# nominal radial/arc discretisation of the 0.5 mm air gap.
AIRGAP_MESH_LEVELS = (
    MeshConfig(name="gap3", airgap_radial_layers=3, airgap_arc_deg=1.0),
    MeshConfig(name="gap5", airgap_radial_layers=5, airgap_arc_deg=0.5),
    MeshConfig(name="gap8", airgap_radial_layers=8, airgap_arc_deg=0.25),
    MeshConfig(name="gap12", airgap_radial_layers=12, airgap_arc_deg=0.15),
)

# Angle and torque-integral studies are intentionally held at one fixed mesh.
# The separate four-level gap study then measures the remaining mesh error
# without forcing an already completed angular/integration study to change too.
ANGLE_INTEGRAL_MESH = AIRGAP_MESH_LEVELS[2]


def mesh_config_by_name(name: str) -> MeshConfig:
    for mesh in MESH_LEVELS:
        if mesh.name == name:
            return mesh
    choices = ", ".join(mesh.name for mesh in MESH_LEVELS)
    raise ValueError(f"Unknown mesh level {name!r}; choose one of: {choices}.")


# ── geometry helpers ───────────────────────────────────────────────

def polar_xy(r_mm: float, angle_deg: float) -> tuple[float, float]:
    a = math.radians(angle_deg)
    return r_mm * math.cos(a), r_mm * math.sin(a)


def parse_list(raw) -> list[float]:
    if isinstance(raw, str):
        items = [raw]
    else:
        items = list(raw)
    values = []
    for item in items:
        for token in re.split(r"[\s,]+", str(item).strip()):
            if token:
                values.append(float(token))
    return values


def winding_layout(machine: Machine) -> tuple[tuple[str, int], ...]:
    if (machine.slots, 2 * machine.pole_pairs) != (18, 20):
        raise ValueError(
            "Only the validated single-layer 18-slot / 20-pole winding is "
            f"supported; got {machine.slots} slots / {2 * machine.pole_pairs} poles."
        )
    return WINDING_18S20P


def winding_fundamentals(
    machine: Machine,
    layout: tuple[tuple[str, int], ...] | None = None,
) -> dict[str, complex]:
    sides = winding_layout(machine) if layout is None else layout
    if len(sides) != machine.slots:
        raise ValueError("Winding layout length must equal the slot count.")
    slot_pitch_rad = 2.0 * math.pi / machine.slots
    return {
        phase: sum(
            sign * complex(
                math.cos(machine.pole_pairs * slot * slot_pitch_rad),
                math.sin(machine.pole_pairs * slot * slot_pitch_rad),
            )
            for slot, (slot_phase, sign) in enumerate(sides)
            if slot_phase == phase
        )
        for phase in PHASES
    }


def _wrap_degrees(angle_deg: float) -> float:
    return (angle_deg + 180.0) % 360.0 - 180.0


def validate_winding_layout(machine: Machine) -> dict[str, float]:
    layout = winding_layout(machine)
    fundamentals = winding_fundamentals(machine, layout)
    counts = {phase: sum(1 for item in layout if item[0] == phase) for phase in PHASES}
    net_sides = {
        phase: sum(sign for slot_phase, sign in layout if slot_phase == phase)
        for phase in PHASES
    }
    magnitudes = {phase: abs(fundamentals[phase]) for phase in PHASES}
    axes = {
        phase: math.degrees(math.atan2(fundamentals[phase].imag, fundamentals[phase].real))
        for phase in PHASES
    }
    if len(set(counts.values())) != 1 or any(net_sides.values()):
        raise ValueError("Each phase must have equal slot sides and zero net axial current.")
    reference = magnitudes["A"]
    if reference <= 0.0 or any(abs(value / reference - 1.0) > 1e-9 for value in magnitudes.values()):
        raise ValueError("The three phase fundamental winding magnitudes are not balanced.")
    if abs(_wrap_degrees(axes["B"] - axes["A"] + 120.0)) > 1e-9:
        raise ValueError("Phase B winding axis must lag phase A by 120 electrical degrees.")
    if abs(_wrap_degrees(axes["C"] - axes["A"] - 120.0)) > 1e-9:
        raise ValueError("Phase C winding axis must lead phase A by 120 electrical degrees.")
    return {
        "phase_axis_a_deg": axes["A"],
        "phase_axis_b_deg": axes["B"],
        "phase_axis_c_deg": axes["C"],
        "fundamental_winding_factor": reference / counts["A"],
    }


def electrical_zero_offset_rad(machine: Machine) -> float:
    """Park angle that aligns +d with the no-load PM flux at rotor position zero."""
    phase_axis_a_deg = validate_winding_layout(machine)["phase_axis_a_deg"]
    # FEMM circuit flux linkage follows the conjugate spatial sequence of the
    # coil-side MMF vector and is 90 electrical degrees behind that vector.
    return math.radians(phase_axis_a_deg - 90.0)


def electrical_angle_rad(machine: Machine, rotor_mechanical_deg: float) -> float:
    return (
        -math.radians(machine.pole_pairs * rotor_mechanical_deg)
        + electrical_zero_offset_rad(machine)
    )


def dq_to_abc(id_a: float, iq_a: float, theta_e_rad: float) -> tuple[float, float, float]:
    phase_angles = (
        theta_e_rad,
        theta_e_rad - 2.0 * math.pi / 3.0,
        theta_e_rad + 2.0 * math.pi / 3.0,
    )
    return tuple(
        id_a * math.cos(angle) + iq_a * math.sin(angle)
        for angle in phase_angles
    )


def abc_to_dq(a: float, b: float, c: float, theta_e_rad: float) -> tuple[float, float]:
    phase_angles = (
        theta_e_rad,
        theta_e_rad - 2.0 * math.pi / 3.0,
        theta_e_rad + 2.0 * math.pi / 3.0,
    )
    values = (a, b, c)
    d_value = (2.0 / 3.0) * sum(value * math.cos(angle) for value, angle in zip(values, phase_angles))
    q_value = (2.0 / 3.0) * sum(value * math.sin(angle) for value, angle in zip(values, phase_angles))
    return d_value, q_value


def magnet_edge_angles(machine: Machine) -> tuple[float, ...]:
    pole_pitch = 360.0 / (2 * machine.pole_pairs)
    half_span = 0.5 * pole_pitch * machine.magnet_arc_ratio
    return tuple(sorted({
        round((center + edge) % 360.0, 12)
        for pole in range(2 * machine.pole_pairs)
        for center in (pole * pole_pitch,)
        for edge in (-half_span, half_span)
    }))


def add_grouped_node(x: float, y: float, group: int) -> None:
    femm.mi_addnode(x, y)
    femm.mi_selectnode(x, y)
    femm.mi_setnodeprop("<None>", group)
    femm.mi_clearselected()


def add_segmented_circle(
    radius_mm: float,
    break_angles_deg,
    *,
    group: int = 0,
    boundary: str = "<None>",
    maxseg_deg: float = 1.0,
) -> None:
    """Create one non-overlapping closed circle split at the requested angles."""
    angles = sorted({round(float(angle) % 360.0, 12) for angle in break_angles_deg})
    if len(angles) < 2:
        raise ValueError("A segmented circle requires at least two distinct angles.")
    points = [polar_xy(radius_mm, angle) for angle in angles]
    for x, y in points:
        add_grouped_node(x, y, group)
    for index, angle_start in enumerate(angles):
        angle_end = angles[(index + 1) % len(angles)]
        sweep = (angle_end - angle_start) % 360.0
        if sweep <= 0.0 or sweep > 180.0:
            raise ValueError(f"Invalid circle arc sweep: {sweep} degrees.")
        x1, y1 = points[index]
        x2, y2 = points[(index + 1) % len(points)]
        femm.mi_addarc(x1, y1, x2, y2, sweep, maxseg_deg)
        mx, my = polar_xy(radius_mm, angle_start + 0.5 * sweep)
        femm.mi_selectarcsegment(mx, my)
        femm.mi_setarcsegmentprop(maxseg_deg, boundary, 0, group)
        femm.mi_clearselected()


def add_quarter_arc_circle(
    radius_mm: float,
    group: int = 0,
    boundary: str = "<None>",
    maxseg_deg: float = 1.0,
) -> None:
    add_segmented_circle(
        radius_mm,
        (0.0, 90.0, 180.0, 270.0),
        group=group,
        boundary=boundary,
        maxseg_deg=maxseg_deg,
    )


def add_grouped_segment(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    group: int,
    mesh_size_mm: float,
) -> None:
    femm.mi_addsegment(x1, y1, x2, y2)
    femm.mi_selectsegment(0.5 * (x1 + x2), 0.5 * (y1 + y2))
    femm.mi_setsegmentprop("<None>", mesh_size_mm, 1, 0, group)
    femm.mi_clearselected()


def add_wedge_annulus(
    r_in: float,
    r_out: float,
    a1: float,
    a2: float,
    group: int = 0,
    mesh_size_mm: float = 1.0,
    arc_max_deg: float = 1.0,
) -> tuple[float, float]:
    x1, y1 = polar_xy(r_in, a1)
    x2, y2 = polar_xy(r_in, a2)
    x3, y3 = polar_xy(r_out, a2)
    x4, y4 = polar_xy(r_out, a1)
    for x, y in ((x1, y1), (x2, y2), (x3, y3), (x4, y4)):
        add_grouped_node(x, y, group)
    sweep = (a2 - a1) % 360.0
    if sweep <= 0.0:
        sweep += 360.0
    femm.mi_addarc(x1, y1, x2, y2, sweep, arc_max_deg)
    femm.mi_addarc(x4, y4, x3, y3, sweep, arc_max_deg)
    add_grouped_segment(x2, y2, x3, y3, group=group, mesh_size_mm=mesh_size_mm)
    add_grouped_segment(x1, y1, x4, y4, group=group, mesh_size_mm=mesh_size_mm)

    r_mid = 0.5 * (r_in + r_out)
    a_mid = a1 + sweep * 0.5
    for radius in (r_in, r_out):
        mx, my = polar_xy(radius, a_mid)
        femm.mi_selectarcsegment(mx, my)
        femm.mi_setarcsegmentprop(arc_max_deg, "<None>", 0, group)
        femm.mi_clearselected()
    return polar_xy(r_mid, a_mid)


def set_block(
    material,
    x,
    y,
    circuit="<None>",
    mag_dir=0.0,
    group=0,
    turns=0,
    mesh_size_mm: float | None = None,
) -> None:
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    automesh = 1 if mesh_size_mm is None else 0
    mesh_size = 0.0 if mesh_size_mm is None else mesh_size_mm
    femm.mi_setblockprop(material, automesh, mesh_size, circuit, mag_dir, group, turns)
    femm.mi_clearselected()


# ── materials ──────────────────────────────────────────────────────

def ensure_materials(machine: Machine):
    """Load Air / M‑19 Steel / Copper from FEMM library;
    manually define NdFeB 40 MGOe (not shipped with FEMM).
    """
    for name in ["Air", machine.steel_material_name, machine.copper_material_name]:
        femm.mi_getmaterial(name)
    mu_r = machine.magnet_relative_permeability
    Br = machine.magnet_br_t
    Hc = Br / (MU0 * mu_r)
    femm.mi_addmaterial(
        "NdFeB 40 MGOe", mu_r, mu_r, Hc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    )


# ── model builder ─────────────────────────────────────────────────

def build_model(
    machine: Machine,
    tmp_fem: str,
    mesh: MeshConfig | None = None,
) -> None:
    mesh = MeshConfig() if mesh is None else mesh
    layout = winding_layout(machine)
    validate_winding_layout(machine)
    femm.newdocument(0)
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, machine.stack_length_mm, 30)
    ensure_materials(machine)

    femm.mi_addcircprop("A", 0.0, 1)
    femm.mi_addcircprop("B", 0.0, 1)
    femm.mi_addcircprop("C", 0.0, 1)
    femm.mi_addboundprop("A0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    # Rotor annulus boundaries are segmented once at every magnet edge. The
    # radial separators then form magnet and interpole-air regions without any
    # coincident duplicate arcs.
    edge_angles = magnet_edge_angles(machine)
    add_quarter_arc_circle(machine.r_shaft_mm, group=ROTOR_GROUP)
    add_segmented_circle(
        machine.r_rotor_mm,
        edge_angles,
        group=ROTOR_GROUP,
        maxseg_deg=mesh.airgap_arc_deg,
    )
    add_segmented_circle(
        machine.r_mag_outer_mm,
        edge_angles,
        group=ROTOR_GROUP,
        maxseg_deg=mesh.airgap_arc_deg,
    )
    add_quarter_arc_circle(
        machine.r_stator_inner_mm,
        group=STATOR_GROUP,
        maxseg_deg=mesh.airgap_arc_deg,
    )
    add_quarter_arc_circle(machine.r_stator_outer_mm, group=STATOR_GROUP)
    add_quarter_arc_circle(machine.r_air_outer_mm, boundary="A0")

    for angle in edge_angles:
        x1, y1 = polar_xy(machine.r_rotor_mm, angle)
        x2, y2 = polar_xy(machine.r_mag_outer_mm, angle)
        add_grouped_segment(
            x1, y1, x2, y2,
            group=ROTOR_GROUP,
            mesh_size_mm=mesh.magnet_size_mm,
        )

    set_block("Air", 0.5 * machine.r_shaft_mm, 0.0,
              mesh_size_mm=mesh.outer_air_size_mm)
    set_block(machine.steel_material_name, 0.5 * (machine.r_shaft_mm + machine.r_rotor_mm), 0.0,
              group=ROTOR_GROUP, mesh_size_mm=mesh.rotor_size_mm)
    set_block(
        "Air",
        0.5 * (machine.r_mag_outer_mm + machine.r_stator_inner_mm),
        0.0,
        group=AIRGAP_GROUP,
        mesh_size_mm=mesh.airgap_size_mm(machine),
    )
    set_block(machine.steel_material_name, 0.5 * (machine.r_slot_outer_mm + machine.r_stator_outer_mm),
              12.0, group=STATOR_GROUP, mesh_size_mm=mesh.stator_size_mm)
    set_block("Air", 0.5 * (machine.r_air_outer_mm + machine.r_stator_outer_mm), 0.0,
              mesh_size_mm=mesh.outer_air_size_mm)

    # magnets
    pole_count = 2 * machine.pole_pairs
    pole_pitch = 360.0 / pole_count
    for k in range(pole_count):
        center = k * pole_pitch
        x, y = polar_xy(
            0.5 * (machine.r_rotor_mm + machine.r_mag_outer_mm),
            center,
        )
        mag_dir = center if (k % 2 == 0) else (center + 180.0)
        set_block(
            "NdFeB 40 MGOe", x, y,
            mag_dir=mag_dir,
            group=ROTOR_GROUP,
            mesh_size_mm=mesh.magnet_size_mm,
        )

    # air in rotor gaps between magnets
    r_gap_mid = 0.5 * (machine.r_rotor_mm + machine.r_mag_outer_mm)
    for k in range(pole_count):
        gx, gy = polar_xy(r_gap_mid, (k + 0.5) * pole_pitch)
        set_block("Air", gx, gy, group=ROTOR_GROUP,
                  mesh_size_mm=mesh.magnet_size_mm)

    # stator slots
    slot_pitch = 360.0 / machine.slots
    slot_span = slot_pitch * 0.55
    r_slot_in = machine.r_stator_inner_mm + 0.2
    r_slot_out = machine.r_slot_outer_mm
    for idx, (phase, sign) in enumerate(layout):
        center = idx * slot_pitch
        a1, a2 = center - 0.5 * slot_span, center + 0.5 * slot_span
        x, y = add_wedge_annulus(
            r_slot_in, r_slot_out, a1, a2,
            group=COIL_GROUP,
            mesh_size_mm=mesh.coil_size_mm,
            arc_max_deg=1.0,
        )
        set_block(machine.copper_material_name, x, y, circuit=phase, group=COIL_GROUP,
                  turns=sign * machine.turns_per_slot,
                  mesh_size_mm=mesh.coil_size_mm)

    femm.mi_saveas(tmp_fem)


# ── low‑level helpers ──────────────────────────────────────────────

def apply_phase_currents(ia: float, ib: float, ic: float) -> tuple[float, float, float]:
    femm.mi_modifycircprop("A", 1, ia)
    femm.mi_modifycircprop("B", 1, ib)
    femm.mi_modifycircprop("C", 1, ic)
    return ia, ib, ic


def set_phase_currents(iq_peak: float, theta_e_rad: float) -> tuple[float, float, float]:
    """Compatibility helper: apply positive q-axis peak current."""
    return apply_phase_currents(*dq_to_abc(0.0, iq_peak, theta_e_rad))


def set_dq_currents(
    machine: Machine,
    id_a: float,
    iq_a: float,
    theta_e_rad: float,
) -> tuple[float, float, float]:
    """Apply amplitude-invariant dq peak currents using the validated phase order."""
    winding_layout(machine)
    return apply_phase_currents(*dq_to_abc(id_a, iq_a, theta_e_rad))


def estimate_core_loss(machine: Machine, loss: LossThermal, rpm: float):
    if loss.core_loss_ref_w <= 0.0 or rpm == 0.0:
        return 0.0
    probe_r = 0.5 * (machine.r_slot_outer_mm + machine.r_stator_outer_mm)
    b_max = 0.0
    for deg in [15, 45, 75, 105, 135, 165]:
        bx, by = femm.mo_getb(*polar_xy(probe_r, deg))
        b_max = max(b_max, math.hypot(bx, by))
    f_e = machine.pole_pairs * rpm / 60.0
    return loss.core_loss_ref_w * (
        abs(f_e) / loss.core_loss_ref_electrical_hz
    ) ** loss.core_frequency_exponent * (
        b_max / loss.core_flux_density_ref_t
    ) ** 2


def estimate_magnet_loss(machine: Machine, loss: LossThermal, rpm: float) -> float:
    f_e = machine.pole_pairs * abs(rpm) / 60.0
    return loss.magnet_loss_ref_w * (
        f_e / loss.magnet_loss_ref_electrical_hz
    ) ** loss.magnet_frequency_exponent


def _get_circuit_flux(phase: str):
    """Return flux linkage of *phase* from loaded solution."""
    props = femm.mo_getcircuitproperties(phase)
    # props = (current, voltage, flux_linkage)
    return props[2]


def _close_solution() -> None:
    try:
        femm.mo_close()
    except Exception:
        pass


def _close_model_document() -> None:
    _close_solution()
    try:
        femm.mi_close()
    except Exception:
        pass


def _cleanup_model_files(tmp_fem: str) -> None:
    base = Path(tmp_fem)
    for path in (base, base.with_suffix(".ans")):
        try:
            path.unlink(missing_ok=True)
        except OSError:
            pass


@contextmanager
def temporary_femm_model(machine: Machine, mesh: MeshConfig):
    tmp_fem = get_tmp_fem()
    try:
        _close_model_document()
        build_model(machine, tmp_fem, mesh=mesh)
        yield tmp_fem
    finally:
        _close_model_document()
        _cleanup_model_files(tmp_fem)


def rotate_rotor(delta_mechanical_deg: float) -> None:
    if abs(delta_mechanical_deg) <= 1e-15:
        return
    _close_solution()
    femm.mi_selectgroup(ROTOR_GROUP)
    femm.mi_moverotate(0.0, 0.0, delta_mechanical_deg)
    femm.mi_clearselected()


@dataclass(frozen=True)
class OperatingPoint:
    torque_nm: float
    psi_d_wb: float
    psi_q_wb: float


@dataclass(frozen=True)
class ElectromagneticSample:
    rotor_angle_deg: float
    theta_e_deg: float
    psi_d0_wb: float
    psi_q0_wb: float
    cogging_torque_nm: float
    torque_q_pos_nm: float
    torque_q_neg_nm: float
    ld_h: float
    lq_h: float


@dataclass(frozen=True)
class AcceptanceTolerances:
    max_psi_q_ratio: float = 0.05
    max_torque_symmetry_error: float = 0.10
    max_torque_slope_error: float = 0.20
    max_saliency_ratio: float = 0.20
    min_load_to_cogging_ratio: float = 5.0


def solve_operating_point(
    machine: Machine,
    rotor_angle_deg: float,
    id_a: float,
    iq_a: float,
) -> OperatingPoint:
    _close_solution()
    theta_e = electrical_angle_rad(machine, rotor_angle_deg)
    set_dq_currents(machine, id_a, iq_a, theta_e)
    femm.mi_analyze(1)
    femm.mi_loadsolution()
    femm.mo_groupselectblock(ROTOR_GROUP)
    torque = float(femm.mo_blockintegral(22))
    femm.mo_clearblock()
    psi_d, psi_q = abc_to_dq(
        _get_circuit_flux("A"),
        _get_circuit_flux("B"),
        _get_circuit_flux("C"),
        theta_e,
    )
    return OperatingPoint(torque, psi_d, psi_q)


def run_electromagnetic_samples(
    machine: Machine,
    mesh: MeshConfig,
    *,
    iq_test_a: float,
    delta_current_a: float,
    num_steps: int,
    span_mechanical_deg: float | None = None,
    progress_label: str = "Electromagnetic validation",
    minimum_steps: int = 3,
) -> list[ElectromagneticSample]:
    if iq_test_a <= 0.0 or delta_current_a <= 0.0:
        raise ValueError("iq_test_a and delta_current_a must be positive.")
    if minimum_steps < 1:
        raise ValueError("minimum_steps must be positive.")
    if num_steps < minimum_steps:
        raise ValueError(f"num_steps must be at least {minimum_steps}.")
    if span_mechanical_deg is None:
        span_mechanical_deg = 360.0 / machine.pole_pairs
    if span_mechanical_deg <= 0.0:
        raise ValueError("span_mechanical_deg must be positive.")

    angles = np.linspace(0.0, span_mechanical_deg, num_steps, endpoint=False)
    samples: list[ElectromagneticSample] = []
    previous_angle = 0.0
    with temporary_femm_model(machine, mesh):
        for angle in tqdm(angles, desc=progress_label, leave=False):
            rotor_angle = float(angle)
            rotate_rotor(rotor_angle - previous_angle)
            previous_angle = rotor_angle

            zero = solve_operating_point(machine, rotor_angle, 0.0, 0.0)
            q_pos = solve_operating_point(machine, rotor_angle, 0.0, iq_test_a)
            q_neg = solve_operating_point(machine, rotor_angle, 0.0, -iq_test_a)
            d_delta_pos = solve_operating_point(
                machine, rotor_angle, delta_current_a, 0.0
            )
            d_delta_neg = solve_operating_point(
                machine, rotor_angle, -delta_current_a, 0.0
            )
            if math.isclose(iq_test_a, delta_current_a, rel_tol=0.0, abs_tol=1e-12):
                q_delta_pos, q_delta_neg = q_pos, q_neg
            else:
                q_delta_pos = solve_operating_point(
                    machine, rotor_angle, 0.0, delta_current_a
                )
                q_delta_neg = solve_operating_point(
                    machine, rotor_angle, 0.0, -delta_current_a
                )

            ld_h = (
                d_delta_pos.psi_d_wb - d_delta_neg.psi_d_wb
            ) / (2.0 * delta_current_a)
            lq_h = (
                q_delta_pos.psi_q_wb - q_delta_neg.psi_q_wb
            ) / (2.0 * delta_current_a)
            samples.append(
                ElectromagneticSample(
                    rotor_angle_deg=rotor_angle,
                    theta_e_deg=math.degrees(electrical_angle_rad(machine, rotor_angle)),
                    psi_d0_wb=zero.psi_d_wb,
                    psi_q0_wb=zero.psi_q_wb,
                    cogging_torque_nm=zero.torque_nm,
                    torque_q_pos_nm=q_pos.torque_nm,
                    torque_q_neg_nm=q_neg.torque_nm,
                    ld_h=ld_h,
                    lq_h=lq_h,
                )
            )
    return samples


def evaluate_electromagnetic_acceptance(
    machine: Machine,
    samples: list[ElectromagneticSample],
    iq_test_a: float,
    tolerances: AcceptanceTolerances | None = None,
) -> dict[str, dict[str, float | str | bool]]:
    if not samples:
        raise ValueError("At least one electromagnetic sample is required.")
    tolerances = AcceptanceTolerances() if tolerances is None else tolerances
    psi_d0 = np.array([sample.psi_d0_wb for sample in samples])
    psi_q0 = np.array([sample.psi_q0_wb for sample in samples])
    torque_0 = np.array([sample.cogging_torque_nm for sample in samples])
    torque_pos = np.array([sample.torque_q_pos_nm for sample in samples])
    torque_neg = np.array([sample.torque_q_neg_nm for sample in samples])
    ld_arr = np.array([sample.ld_h for sample in samples])
    lq_arr = np.array([sample.lq_h for sample in samples])

    psi_pm = float(np.mean(psi_d0))
    psi_q_ratio = float(np.sqrt(np.mean(psi_q0**2)) / max(abs(psi_pm), 1e-12))
    loaded_pos = torque_pos - torque_0
    loaded_neg = torque_neg - torque_0
    symmetry_scale = max(
        0.5 * (abs(float(np.mean(loaded_pos))) + abs(float(np.mean(loaded_neg)))),
        1e-12,
    )
    symmetry_error = abs(float(np.mean(loaded_pos + loaded_neg))) / symmetry_scale
    torque_slope = float(np.mean(torque_pos - torque_neg)) / (2.0 * iq_test_a)
    expected_slope = 1.5 * machine.pole_pairs * psi_pm
    slope_error = abs(torque_slope - expected_slope) / max(abs(expected_slope), 1e-12)
    ld_mean = float(np.mean(ld_arr))
    lq_mean = float(np.mean(lq_arr))
    saliency_ratio = abs(ld_mean - lq_mean) / max(
        0.5 * (abs(ld_mean) + abs(lq_mean)), 1e-12
    )
    loaded_mean = 0.5 * float(np.mean(loaded_pos - loaded_neg))
    cogging_pp = float(np.ptp(torque_0))
    load_to_cogging = abs(loaded_mean) / max(cogging_pp, 1e-12)

    return {
        "psi_d0_positive": {
            "value": psi_pm,
            "criterion": ">",
            "limit": 0.0,
            "passed": psi_pm > 0.0,
        },
        "psi_q0_rms_over_psi_d0": {
            "value": psi_q_ratio,
            "criterion": "<=",
            "limit": tolerances.max_psi_q_ratio,
            "passed": psi_q_ratio <= tolerances.max_psi_q_ratio,
        },
        "torque_symmetry_error": {
            "value": symmetry_error,
            "criterion": "<=",
            "limit": tolerances.max_torque_symmetry_error,
            "passed": symmetry_error <= tolerances.max_torque_symmetry_error,
        },
        "torque_slope_relative_error": {
            "value": slope_error,
            "criterion": "<=",
            "limit": tolerances.max_torque_slope_error,
            "passed": slope_error <= tolerances.max_torque_slope_error,
        },
        "ld_lq_saliency_ratio": {
            "value": saliency_ratio,
            "criterion": "<=",
            "limit": tolerances.max_saliency_ratio,
            "passed": saliency_ratio <= tolerances.max_saliency_ratio,
        },
        "load_to_cogging_ratio": {
            "value": load_to_cogging,
            "criterion": ">=",
            "limit": tolerances.min_load_to_cogging_ratio,
            "passed": load_to_cogging >= tolerances.min_load_to_cogging_ratio,
        },
        "psi_pm_wb": {"value": psi_pm, "criterion": "report", "limit": math.nan, "passed": True},
        "torque_slope_nm_per_a": {
            "value": torque_slope, "criterion": "report", "limit": math.nan, "passed": True
        },
        "expected_torque_slope_nm_per_a": {
            "value": expected_slope, "criterion": "report", "limit": math.nan, "passed": True
        },
        "ld_mean_h": {"value": ld_mean, "criterion": "report", "limit": math.nan, "passed": True},
        "lq_mean_h": {"value": lq_mean, "criterion": "report", "limit": math.nan, "passed": True},
        "cogging_peak_to_peak_nm": {
            "value": cogging_pp, "criterion": "report", "limit": math.nan, "passed": True
        },
        "loaded_torque_mean_nm": {
            "value": loaded_mean, "criterion": "report", "limit": math.nan, "passed": True
        },
    }


def save_electromagnetic_acceptance(
    samples: list[ElectromagneticSample],
    metrics: dict[str, dict[str, float | str | bool]],
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_table = np.array([
        (
            sample.rotor_angle_deg,
            sample.theta_e_deg,
            sample.psi_d0_wb,
            sample.psi_q0_wb,
            sample.cogging_torque_nm,
            sample.torque_q_pos_nm,
            sample.torque_q_neg_nm,
            sample.ld_h,
            sample.lq_h,
        )
        for sample in samples
    ])
    np.savetxt(
        out_dir / "electromagnetic_validation_samples.csv",
        sample_table,
        delimiter=",",
        header=(
            "rotor_angle_deg,theta_e_deg,psi_d0_wb,psi_q0_wb,cogging_torque_nm,"
            "torque_q_pos_nm,torque_q_neg_nm,Ld_H,Lq_H"
        ),
        comments="",
    )
    with (out_dir / "electromagnetic_acceptance.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(("metric", "value", "criterion", "limit", "passed"))
        for name, metric in metrics.items():
            writer.writerow((
                name,
                metric["value"],
                metric["criterion"],
                metric["limit"],
                metric["passed"],
            ))
    json_metrics = {}
    for name, metric in metrics.items():
        limit = float(metric["limit"])
        json_metrics[name] = {
            "value": float(metric["value"]),
            "criterion": str(metric["criterion"]),
            "limit": limit if math.isfinite(limit) else None,
            "passed": bool(metric["passed"]),
        }
    (out_dir / "electromagnetic_acceptance_metrics.json").write_text(
        json.dumps(
            {"schema_version": 1, "metrics": json_metrics},
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )


def run_electromagnetic_acceptance(
    machine: Machine,
    out_dir: Path,
    *,
    iq_test_a: float = 5.0,
    delta_current_a: float = 1.0,
    num_steps: int = 12,
    mesh: MeshConfig | None = None,
) -> dict[str, dict[str, float | str | bool]]:
    mesh = MeshConfig() if mesh is None else mesh
    samples = run_electromagnetic_samples(
        machine,
        mesh,
        iq_test_a=iq_test_a,
        delta_current_a=delta_current_a,
        num_steps=num_steps,
    )
    metrics = evaluate_electromagnetic_acceptance(machine, samples, iq_test_a)
    save_electromagnetic_acceptance(samples, metrics, out_dir)
    for name, metric in metrics.items():
        if metric["criterion"] != "report":
            status = "PASS" if metric["passed"] else "FAIL"
            print(f"  {status:4s} {name}: {float(metric['value']):.6g}")
    return metrics


def run_electromagnetic_mesh_check(
    machine: Machine,
    out_dir: Path,
    *,
    iq_test_a: float = 1.0,
    delta_current_a: float = 1.0,
    mesh: MeshConfig | None = None,
) -> dict[str, float | str]:
    """Evaluate one rotor angle for medium/fine global-mesh comparison."""
    mesh = MeshConfig() if mesh is None else mesh
    sample = run_electromagnetic_samples(
        machine,
        mesh,
        iq_test_a=iq_test_a,
        delta_current_a=delta_current_a,
        num_steps=1,
        minimum_steps=1,
        progress_label=f"Global mesh check ({mesh.name})",
    )[0]
    loaded_pos = sample.torque_q_pos_nm - sample.cogging_torque_nm
    loaded_neg = sample.torque_q_neg_nm - sample.cogging_torque_nm
    symmetry_scale = max(0.5 * (abs(loaded_pos) + abs(loaded_neg)), 1e-12)
    torque_slope = (
        sample.torque_q_pos_nm - sample.torque_q_neg_nm
    ) / (2.0 * iq_test_a)
    expected_slope = 1.5 * machine.pole_pairs * sample.psi_d0_wb
    report: dict[str, float | str] = {
        "mesh_level": mesh.name,
        "rotor_angle_deg": sample.rotor_angle_deg,
        "theta_e_deg": sample.theta_e_deg,
        "iq_test_peak_a": iq_test_a,
        "delta_current_peak_a": delta_current_a,
        "psi_d0_wb": sample.psi_d0_wb,
        "psi_q0_wb": sample.psi_q0_wb,
        "psi_q0_abs_over_psi_d0": abs(sample.psi_q0_wb)
        / max(abs(sample.psi_d0_wb), 1e-12),
        "cogging_torque_nm": sample.cogging_torque_nm,
        "loaded_torque_pos_nm": loaded_pos,
        "loaded_torque_neg_nm": loaded_neg,
        "torque_symmetry_error": abs(loaded_pos + loaded_neg) / symmetry_scale,
        "torque_slope_nm_per_a": torque_slope,
        "expected_torque_slope_nm_per_a": expected_slope,
        "torque_slope_relative_error": abs(torque_slope - expected_slope)
        / max(abs(expected_slope), 1e-12),
        "ld_h": sample.ld_h,
        "lq_h": sample.lq_h,
        "ld_lq_saliency_ratio": abs(sample.ld_h - sample.lq_h)
        / max(0.5 * (abs(sample.ld_h) + abs(sample.lq_h)), 1e-12),
    }
    path = out_dir / f"electromagnetic_mesh_check_{mesh.name}.json"
    path.write_text(
        json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"  mesh check written to: {path}")
    return report


def compare_electromagnetic_mesh_checks(
    reference: dict[str, float | str],
    target: dict[str, float | str],
    out_dir: Path,
    *,
    max_relative_change: float = 0.01,
) -> dict[str, object]:
    """Compare two one-angle mesh reports and write explicit pass/fail metrics."""
    if max_relative_change <= 0.0:
        raise ValueError("max_relative_change must be positive.")
    reference_name = str(reference["mesh_level"])
    target_name = str(target["mesh_level"])
    if reference_name == target_name:
        raise ValueError("Mesh-check reference and target levels must differ.")
    for name in (
        "rotor_angle_deg",
        "iq_test_peak_a",
        "delta_current_peak_a",
    ):
        if not math.isclose(
            float(reference[name]),
            float(target[name]),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError(f"Mesh-check operating points differ in {name}.")
    metric_names = (
        "psi_d0_wb",
        "torque_slope_nm_per_a",
        "ld_h",
        "lq_h",
    )
    metrics: dict[str, dict[str, float | bool]] = {}
    for name in metric_names:
        reference_value = float(reference[name])
        target_value = float(target[name])
        relative_change = abs(target_value - reference_value) / max(
            abs(target_value), 1e-12
        )
        metrics[name] = {
            "reference_value": reference_value,
            "target_value": target_value,
            "relative_change": relative_change,
            "limit": max_relative_change,
            "passed": relative_change <= max_relative_change,
        }
    payload: dict[str, object] = {
        "schema_version": 1,
        "reference_mesh": reference_name,
        "target_mesh": target_name,
        "rotor_angle_deg": float(reference["rotor_angle_deg"]),
        "metrics": metrics,
        "all_passed": all(bool(metric["passed"]) for metric in metrics.values()),
    }
    path = out_dir / (
        f"electromagnetic_mesh_comparison_{reference_name}_vs_{target_name}.json"
    )
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"  mesh comparison written to: {path}")
    return payload


# ══════════════════════════════════════════════════════════════════
# 1) BASIC: torque / loss / temperature waveforms (original)
# ══════════════════════════════════════════════════════════════════

def run_one_case(
    machine,
    loss,
    cfg,
    rpm,
    iq_peak,
    progress=False,
    mesh: MeshConfig | None = None,
):
    if rpm <= 0.0:
        raise ValueError("rpm must be positive for a time-domain electrical-cycle case.")
    mesh = MeshConfig() if mesh is None else mesh
    f_e = machine.pole_pairs * rpm / 60.0
    dt = 1.0 / (f_e * cfg.points_per_electrical_cycle)
    dtheta_m_deg = 360.0 / (
        machine.pole_pairs * cfg.points_per_electrical_cycle
    )
    rows = []

    with temporary_femm_model(machine, mesh):
        iterator = range(cfg.points_per_electrical_cycle)
        if progress:
            iterator = tqdm(
                iterator,
                desc=f"FEA (rpm={rpm:n}, Iq={iq_peak:n})",
                leave=False,
            )

        for k in iterator:
            if k > 0:
                rotate_rotor(dtheta_m_deg)
            rotor_angle_deg = k * dtheta_m_deg
            point = solve_operating_point(machine, rotor_angle_deg, 0.0, iq_peak)
            theta_e = electrical_angle_rad(machine, rotor_angle_deg)
            ia, ib, ic = dq_to_abc(0.0, iq_peak, theta_e)
            pfe = estimate_core_loss(machine, loss, rpm)
            ppm = estimate_magnet_loss(machine, loss, rpm)
            temp = loss.t_amb
            for _ in range(100):
                rs = loss.r_phase_20 * (1.0 + loss.alpha_cu * (temp - 20.0))
                pcu = rs * (ia**2 + ib**2 + ic**2)
                ploss = pcu + pfe + ppm
                case_temp = loss.t_amb + ploss * loss.case_to_ambient_rth_k_per_w
                new_temp = (
                    case_temp + pcu * loss.winding_to_case_rth_k_per_w
                )
                if abs(new_temp - temp) < 1e-7:
                    temp = new_temp
                    break
                temp = 0.5 * temp + 0.5 * new_temp
            else:
                raise RuntimeError(
                    "FEMM steady thermal-resistance iteration did not converge."
                )
            rows.append(
                (
                    k * dt,
                    theta_e,
                    ia,
                    ib,
                    ic,
                    point.torque_nm,
                    pcu,
                    pfe,
                    ppm,
                    ploss,
                    temp,
                )
            )
    return np.array(rows, dtype=float)


def save_case_outputs(rows, out_dir, case_tag, show):
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"femm_waveforms_{case_tag}.csv"
    np.savetxt(csv_path, rows, delimiter=",",
               header="t_s,theta_e_rad,ia_a,ib_a,ic_a,torque_nm,"
                      "pcu_w,pfe_w,ppm_w,ploss_w,winding_temp_steady_c", comments="")

    sample_angle_deg = np.linspace(0.0, 360.0, len(rows) + 1)
    torque = rows[:, 5]
    pcu, pfe, ppm, ploss, temp = rows[:, 6], rows[:, 7], rows[:, 8], rows[:, 9], rows[:, 10]

    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    avg_t = np.mean(torque)
    axs[0].plot(sample_angle_deg, np.append(torque, torque[0]), label="FEMM torque")
    axs[0].axhline(avg_t, color="gray", ls="--", label=f"Avg = {avg_t:.3f} N·m")
    axs[0].set_ylabel("Torque [N·m]"); axs[0].grid(True); axs[0].legend()

    axs[1].plot(sample_angle_deg, np.append(pcu, pcu[0]), label="Copper loss")
    axs[1].plot(sample_angle_deg, np.append(pfe, pfe[0]), label="Core loss")
    axs[1].plot(sample_angle_deg, np.append(ppm, ppm[0]), label="Magnet loss")
    axs[1].plot(sample_angle_deg, np.append(ploss, ploss[0]), "k", label="Total loss")
    axs[1].set_ylabel("Loss [W]"); axs[1].grid(True); axs[1].legend()

    axs[2].plot(
        sample_angle_deg,
        np.append(temp, temp[0]),
        "r",
        label="Steady winding estimate at each rotor position",
    )
    axs[2].set_xlabel("Electrical-cycle sample angle [deg]")
    axs[2].set_ylabel("Temperature [°C]")
    axs[2].grid(True); axs[2].legend()
    fig.tight_layout()
    png_path = out_dir / f"femm_waveforms_{case_tag}.png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig) if not show else plt.show()
    return csv_path, png_path


# ══════════════════════════════════════════════════════════════════
# 2) FIELD MAPS: flux‑density cloud + flux‑line contour
# ══════════════════════════════════════════════════════════════════

def plot_field_maps(
    machine,
    out_dir,
    case_tag,
    show=False,
    *,
    radial_points: int = 48,
    angular_points: int = 180,
):
    """Sample |B| and A on a polar grid from loaded solution; produce two PNGs."""
    if radial_points < 8 or angular_points < 36:
        raise ValueError(
            "Field-map sampling requires at least 8 radial and 36 angular points."
        )
    r_min = machine.r_shaft_mm + 1.0
    r_max = machine.r_stator_outer_mm - 1.0
    nr, ntheta = radial_points, angular_points
    r_arr = np.linspace(r_min, r_max, nr)
    theta_arr = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    R, TH = np.meshgrid(r_arr, theta_arr)
    B_mag = np.zeros_like(R)
    A_val = np.zeros_like(R)

    with tqdm(total=ntheta, desc="Field Maps", leave=False) as pbar:
        for i in range(ntheta):
            for j in range(nr):
                x = R[i, j] * math.cos(TH[i, j])
                y = R[i, j] * math.sin(TH[i, j])
                bx, by = femm.mo_getb(x, y)
                B_mag[i, j] = math.hypot(bx, by)
                A_val[i, j] = femm.mo_geta(x, y)
            pbar.update(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    theta_closed = np.append(theta_arr, 2.0 * np.pi)
    TH_closed, R_closed = np.meshgrid(theta_closed, r_arr, indexing="ij")
    B_closed = np.vstack((B_mag, B_mag[0]))
    A_closed = np.vstack((A_val, A_val[0]))

    # ── flux density cloud ──
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))
    c = ax.pcolormesh(
        TH_closed,
        R_closed,
        B_closed,
        shading="auto",
        cmap="viridis",
    )
    fig.colorbar(c, ax=ax, label="|B| [T]", pad=0.1)
    ax.set_title(f"Flux Density  ({case_tag})", va="bottom", fontsize=13)
    ax.set_yticklabels([])
    fig.savefig(out_dir / f"field_density_{case_tag}.png", dpi=150, bbox_inches="tight")
    plt.close(fig) if not show else plt.show()

    # ── flux lines (iso‑A contour) ──
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))
    if float(np.ptp(A_closed)) > 1e-15:
        levels = np.linspace(A_closed.min(), A_closed.max(), 40)
        ax.contour(
            TH_closed,
            R_closed,
            A_closed,
            levels=levels,
            colors="k",
            linewidths=0.5,
        )
    ax.set_title(f"Flux Lines  ({case_tag})", va="bottom", fontsize=13)
    ax.set_yticklabels([])
    fig.savefig(out_dir / f"field_lines_{case_tag}.png", dpi=150, bbox_inches="tight")
    plt.close(fig) if not show else plt.show()


# ══════════════════════════════════════════════════════════════════
# 3) AIRGAP: radial / tangential flux density along air‑gap
# ══════════════════════════════════════════════════════════════════

def plot_airgap_flux_density(machine, out_dir, case_tag, show=False, num_points=360):
    """Extract Bn / Bt along airgap midline from loaded solution."""
    if num_points < 36:
        raise ValueError("Air-gap sampling requires at least 36 angular points.")
    r_gap = machine.r_mag_outer_mm + machine.airgap_mm / 2.0
    angles = np.linspace(0, 360, num_points, endpoint=False)
    bn_arr, bt_arr = [], []
    for deg in tqdm(angles, desc="Airgap Flux", leave=False):
        x, y = polar_xy(r_gap, deg)
        bx, by = femm.mo_getb(x, y)
        rad = math.radians(deg)
        bn = bx * math.cos(rad) + by * math.sin(rad)   # radial
        bt = -bx * math.sin(rad) + by * math.cos(rad)   # tangential
        bn_arr.append(bn)
        bt_arr.append(bt)
    bn_arr, bt_arr = np.array(bn_arr), np.array(bt_arr)

    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"airgap_B_{case_tag}.csv"
    np.savetxt(csv_path, np.column_stack([angles, bn_arr, bt_arr]),
               delimiter=",", header="angle_deg,Bn_T,Bt_T", comments="")

    fig, ax = plt.subplots(figsize=(10, 4))
    plot_angles = np.append(angles, 360.0)
    ax.plot(plot_angles, np.append(bn_arr, bn_arr[0]), label="Bn (radial)")
    ax.plot(
        plot_angles,
        np.append(bt_arr, bt_arr[0]),
        label="Bt (tangential)",
        alpha=0.7,
    )
    ax.set_xlabel("Mechanical angle [°]"); ax.set_ylabel("Flux density [T]")
    ax.set_title(f"Airgap Flux Density  ({case_tag})")
    ax.grid(True); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"airgap_B_{case_tag}.png", dpi=150)
    plt.close(fig) if not show else plt.show()
    return csv_path, out_dir / f"airgap_B_{case_tag}.png"


def case_value_token(value: float) -> str:
    """Return a readable, collision-free token for a numeric case value."""
    return f"{value:g}".replace("-", "m").replace(".", "p")


def run_field_snapshot(
    machine: Machine,
    out_dir: Path,
    *,
    rpm: float,
    iq_peak_a: float,
    rotor_angle_deg: float,
    analyses: set[str],
    mesh: MeshConfig,
    show: bool = False,
    field_radial_points: int = 48,
    field_angular_points: int = 180,
    airgap_points: int = 360,
) -> OperatingPoint:
    """Solve one explicit rotor-angle snapshot for reproducible field plots."""
    if not analyses & {"field", "airgap"}:
        raise ValueError("A field snapshot requires field and/or airgap analysis.")
    with temporary_femm_model(machine, mesh):
        rotate_rotor(rotor_angle_deg)
        point = solve_operating_point(
            machine, rotor_angle_deg, 0.0, iq_peak_a
        )
        case_tag = (
            f"rpm{case_value_token(rpm)}_"
            f"iq{case_value_token(iq_peak_a)}_"
            f"a{case_value_token(rotor_angle_deg)}"
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        metadata_path = out_dir / f"field_snapshot_{case_tag}.csv"
        with metadata_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                (
                    "rpm_context",
                    "rotor_angle_deg",
                    "id_peak_a",
                    "iq_peak_a",
                    "torque_nm",
                    "psi_d_wb",
                    "psi_q_wb",
                    "mesh",
                )
            )
            writer.writerow(
                (
                    rpm,
                    rotor_angle_deg,
                    0.0,
                    iq_peak_a,
                    point.torque_nm,
                    point.psi_d_wb,
                    point.psi_q_wb,
                    mesh.name,
                )
            )
        if "field" in analyses:
            plot_field_maps(
                machine,
                out_dir,
                case_tag,
                show,
                radial_points=field_radial_points,
                angular_points=field_angular_points,
            )
        if "airgap" in analyses:
            plot_airgap_flux_density(
                machine,
                out_dir,
                case_tag,
                show,
                num_points=airgap_points,
            )
    return point


# ══════════════════════════════════════════════════════════════════
# 4) COGGING TORQUE
# ══════════════════════════════════════════════════════════════════

def run_cogging_torque(
    machine,
    out_dir,
    num_steps=72,
    show=False,
    mesh: MeshConfig | None = None,
):
    """Rotate through one true cogging period and record zero-current torque."""
    if num_steps < 3:
        raise ValueError("num_steps must be at least 3.")
    mesh = MeshConfig() if mesh is None else mesh
    period_deg = cogging_period_deg(machine)
    step_deg = period_deg / num_steps
    angles, torques = [], []
    previous_angle = 0.0
    with temporary_femm_model(machine, mesh):
        for k in tqdm(range(num_steps), desc="Cogging Torque"):
            angle = k * step_deg
            rotate_rotor(angle - previous_angle)
            previous_angle = angle
            point = solve_operating_point(machine, angle, 0.0, 0.0)
            angles.append(angle)
            torques.append(point.torque_nm)
    angles, torques = np.array(angles), np.array(torques)

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_dir / "cogging_torque.csv",
               np.column_stack([angles, torques]),
               delimiter=",", header="rotor_angle_deg,cogging_torque_nm", comments="")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(angles, torques, "b-o", markersize=3)
    ax.axhline(0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("Rotor angle [°]"); ax.set_ylabel("Cogging torque [N·m]")
    ax.set_title(f"Cogging Torque (period = {period_deg:.3f}° mechanical)")
    ax.grid(True); fig.tight_layout()
    fig.savefig(out_dir / "cogging_torque.png", dpi=150)
    plt.close(fig) if not show else plt.show()
    print(f"  Cogging torque peak‑to‑peak: {torques.max() - torques.min():.4f} N·m")


# ══════════════════════════════════════════════════════════════════
# 5) INDUCTANCE: Ld / Lq vs rotor position
# ══════════════════════════════════════════════════════════════════

def run_inductance_analysis(
    machine,
    out_dir,
    test_current=1.0,
    num_steps=37,
    show=False,
    mesh: MeshConfig | None = None,
):
    """Compute incremental Ld/Lq with centered +/- current perturbations."""
    if test_current <= 0.0:
        raise ValueError("test_current must be positive.")
    mesh = MeshConfig() if mesh is None else mesh
    samples = run_electromagnetic_samples(
        machine,
        mesh,
        iq_test_a=test_current,
        delta_current_a=test_current,
        num_steps=num_steps,
        progress_label="Incremental inductance",
    )
    pos_arr = np.array([sample.rotor_angle_deg for sample in samples])
    theta_e_arr = np.array([sample.theta_e_deg for sample in samples])
    psi_pm_arr = np.array([sample.psi_d0_wb for sample in samples])
    psi_q0_arr = np.array([sample.psi_q0_wb for sample in samples])
    ld_arr = np.array([sample.ld_h for sample in samples])
    lq_arr = np.array([sample.lq_h for sample in samples])

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        out_dir / "inductance.csv",
        np.column_stack([
            pos_arr,
            theta_e_arr,
            psi_pm_arr,
            psi_q0_arr,
            ld_arr * 1e3,
            lq_arr * 1e3,
        ]),
        delimiter=",",
        header=(
            "rotor_mechanical_angle_deg,electrical_angle_deg,psi_d0_wb,psi_q0_wb,"
            "Ld_incremental_mH,Lq_incremental_mH"
        ),
        comments="",
    )

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(theta_e_arr, ld_arr * 1e3, "b-o", markersize=3, label="Ld incremental")
    ax.plot(theta_e_arr, lq_arr * 1e3, "r-s", markersize=3, label="Lq incremental")
    ax.set_xlabel("Electrical position [deg]")
    ax.set_ylabel("Incremental inductance [mH]")
    ax.set_title(f"d/q Incremental Inductance (delta current = {test_current} A)")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(out_dir / "inductance.png", dpi=150)
    plt.close(fig) if not show else plt.show()
    print(
        f"  Ld avg = {np.mean(ld_arr)*1e3:.3f} mH, "
        f"Lq avg = {np.mean(lq_arr)*1e3:.3f} mH"
    )
    return float(np.mean(ld_arr)), float(np.mean(lq_arr))


def cogging_period_deg(machine: Machine) -> float:
    return 360.0 / math.lcm(machine.slots, 2 * machine.pole_pairs)


def summarize_torque_sweep(rows: np.ndarray) -> dict[str, float]:
    if rows.ndim != 2 or rows.shape[1] != 3 or len(rows) < 3:
        raise ValueError("Torque sweep rows must be an N x 3 array with N >= 3.")
    cogging = rows[:, 1]
    loaded = rows[:, 2] - cogging
    return {
        "loaded_mean_nm": float(np.mean(loaded)),
        "loaded_ripple_pp_nm": float(np.ptp(loaded)),
        "cogging_mean_nm": float(np.mean(cogging)),
        "cogging_pp_nm": float(np.ptp(cogging)),
    }


def maxwell_airgap_torque_from_samples(
    machine: Machine,
    radius_mm: float,
    br_t: np.ndarray,
    bt_t: np.ndarray,
) -> float:
    """Integrate planar Maxwell shear stress on a circular air-gap contour."""
    br = np.asarray(br_t, dtype=float)
    bt = np.asarray(bt_t, dtype=float)
    if br.shape != bt.shape or br.ndim != 1 or len(br) < 3:
        raise ValueError("br_t and bt_t must be equal-length 1-D arrays with N >= 3.")
    if radius_mm <= 0.0:
        raise ValueError("radius_mm must be positive.")
    radius_m = radius_mm * 1e-3
    depth_m = machine.stack_length_mm * 1e-3
    return float(
        depth_m * radius_m**2 * (2.0 * math.pi / MU0) * np.mean(br * bt)
    )


def sample_airgap_maxwell_torques(
    machine: Machine,
    radius_fraction: float,
    point_levels: tuple[int, ...],
) -> dict[int, float]:
    """Evaluate nested periodic quadratures from one finest air-gap field sample."""
    if not 0.0 < radius_fraction < 1.0:
        raise ValueError("radius_fraction must lie strictly inside the air gap.")
    if (
        len(point_levels) < 2
        or tuple(sorted(set(point_levels))) != point_levels
        or any(points < 3 for points in point_levels)
    ):
        raise ValueError("point_levels must contain at least two increasing counts >= 3.")
    finest_points = point_levels[-1]
    if any(finest_points % points for points in point_levels):
        raise ValueError("Each torque-integral point count must divide the finest count.")

    radius_mm = machine.r_mag_outer_mm + radius_fraction * machine.airgap_mm
    angles = 2.0 * math.pi * np.arange(finest_points) / finest_points
    br = np.empty(finest_points, dtype=float)
    bt = np.empty(finest_points, dtype=float)
    femm.mo_smooth("on")
    for index, angle in enumerate(angles):
        x = radius_mm * math.cos(float(angle))
        y = radius_mm * math.sin(float(angle))
        bx, by = femm.mo_getb(x, y)
        br[index] = float(bx) * math.cos(float(angle)) + float(by) * math.sin(float(angle))
        bt[index] = -float(bx) * math.sin(float(angle)) + float(by) * math.cos(float(angle))

    torques = {}
    for points in point_levels:
        stride = finest_points // points
        torques[points] = maxwell_airgap_torque_from_samples(
            machine,
            radius_mm,
            br[::stride],
            bt[::stride],
        )
    return torques


def run_torque_integral_convergence(
    machine: Machine,
    out_dir: Path,
    *,
    iq_test_a: float,
    mesh: MeshConfig,
    radius_fractions: tuple[float, ...] = (0.25, 0.50, 0.75),
    point_levels: tuple[int, ...] = (180, 360, 720),
) -> tuple[
    dict[str, dict[str, float | str | bool]],
    dict[str, float],
]:
    """Converge Maxwell torque quadrature/path and compare it with WST torque."""
    if iq_test_a <= 0.0:
        raise ValueError("iq_test_a must be positive.")
    if (
        len(radius_fractions) < 3
        or tuple(sorted(set(radius_fractions))) != radius_fractions
        or any(not 0.0 < fraction < 1.0 for fraction in radius_fractions)
    ):
        raise ValueError("radius_fractions must contain at least three increasing gap positions.")

    sampled: dict[str, dict[float, dict[int, float]]] = {}
    wst_torque: dict[str, float] = {}
    with temporary_femm_model(machine, mesh):
        for state, iq_a in (("zero", 0.0), ("loaded", iq_test_a)):
            point = solve_operating_point(machine, 0.0, 0.0, iq_a)
            wst_torque[state] = point.torque_nm
            sampled[state] = {
                fraction: sample_airgap_maxwell_torques(
                    machine, fraction, point_levels
                )
                for fraction in radius_fractions
            }

    incremental = {
        (fraction, points): (
            sampled["loaded"][fraction][points]
            - sampled["zero"][fraction][points]
        )
        for fraction in radius_fractions
        for points in point_levels
    }
    wst_increment = wst_torque["loaded"] - wst_torque["zero"]
    middle_fraction = min(radius_fractions, key=lambda value: abs(value - 0.5))
    finest_points = point_levels[-1]
    quadrature_previous = incremental[(middle_fraction, point_levels[-2])]
    quadrature_reference = incremental[(middle_fraction, finest_points)]
    quadrature_delta = _relative_change(quadrature_previous, quadrature_reference)
    finest_radius_values = np.array(
        [incremental[(fraction, finest_points)] for fraction in radius_fractions]
    )
    radius_reference = float(np.mean(finest_radius_values))
    radius_spread_nm = float(np.ptp(finest_radius_values))
    radius_spread = radius_spread_nm / max(abs(radius_reference), 1e-12)
    wst_difference = _relative_change(quadrature_reference, wst_increment)

    acceptance = {
        "integral_quadrature_relative_change": {
            "value": quadrature_delta, "criterion": "<=", "limit": 0.02,
            "passed": quadrature_delta <= 0.02,
        },
        "integral_radius_relative_spread": {
            "value": radius_spread, "criterion": "<=", "limit": 0.05,
            "passed": radius_spread <= 0.05,
        },
        "maxwell_vs_wst_relative_difference": {
            "value": wst_difference, "criterion": "<=", "limit": 0.10,
            "passed": wst_difference <= 0.10,
        },
    }
    diagnostics = {
        "quadrature_absolute_delta_nm": abs(quadrature_previous - quadrature_reference),
        "radius_spread_nm": radius_spread_nm,
        "wst_increment_nm": wst_increment,
        "maxwell_increment_nm": quadrature_reference,
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "torque_integral_convergence.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow((
            "mesh", "airgap_fraction", "radius_mm", "points",
            "cogging_maxwell_nm", "loaded_total_maxwell_nm",
            "loaded_increment_maxwell_nm", "loaded_increment_wst_nm",
        ))
        for fraction in radius_fractions:
            radius_mm = machine.r_mag_outer_mm + fraction * machine.airgap_mm
            for points in point_levels:
                writer.writerow((
                    mesh.name,
                    fraction,
                    radius_mm,
                    points,
                    sampled["zero"][fraction][points],
                    sampled["loaded"][fraction][points],
                    incremental[(fraction, points)],
                    wst_increment,
                ))
    return acceptance, diagnostics


def run_torque_sweep(
    machine: Machine,
    mesh: MeshConfig,
    *,
    iq_test_a: float,
    num_points: int,
    span_mechanical_deg: float | None = None,
    progress_label: str = "Torque sweep",
) -> tuple[np.ndarray, dict[str, float]]:
    if iq_test_a <= 0.0:
        raise ValueError("iq_test_a must be positive.")
    if num_points < 3:
        raise ValueError("num_points must be at least 3.")
    span = cogging_period_deg(machine) if span_mechanical_deg is None else span_mechanical_deg
    if span <= 0.0:
        raise ValueError("span_mechanical_deg must be positive.")
    angles = np.linspace(0.0, span, num_points, endpoint=False)
    rows = []
    previous_angle = 0.0
    with temporary_femm_model(machine, mesh):
        for angle in tqdm(angles, desc=progress_label, leave=False):
            rotor_angle = float(angle)
            rotate_rotor(rotor_angle - previous_angle)
            previous_angle = rotor_angle
            zero = solve_operating_point(machine, rotor_angle, 0.0, 0.0)
            loaded = solve_operating_point(machine, rotor_angle, 0.0, iq_test_a)
            rows.append((rotor_angle, zero.torque_nm, loaded.torque_nm))
    table = np.array(rows, dtype=float)
    return table, summarize_torque_sweep(table)


def _relative_change(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1e-12)


def torque_sweep_checkpoint_path(
    cache_dir: Path,
    machine: Machine,
    mesh: MeshConfig,
    iq_test_a: float,
    num_points: int,
) -> Path:
    payload = {
        "machine": asdict(machine),
        "mesh": asdict(mesh),
        "iq_test_a": float(iq_test_a),
        "num_points": int(num_points),
        "span_mechanical_deg": cogging_period_deg(machine),
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:16]
    return cache_dir / f"torque_{mesh.name}_n{num_points}_{digest}.csv"


def load_torque_sweep_checkpoint(path: Path, expected_points: int) -> np.ndarray:
    table = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if table.shape != (expected_points, 3) or not np.all(np.isfinite(table)):
        raise ValueError(f"Invalid torque convergence checkpoint: {path}")
    return np.asarray(table, dtype=float)


def run_torque_convergence(
    machine: Machine,
    out_dir: Path,
    *,
    iq_test_a: float = 5.0,
    mesh_angle_points: int = 8,
    angle_point_levels: tuple[int, ...] = (6, 12, 24),
) -> dict[str, dict[str, float | str | bool]]:
    if mesh_angle_points < 3 or any(points < 3 for points in angle_point_levels):
        raise ValueError("All convergence point counts must be at least 3.")
    if tuple(sorted(set(angle_point_levels))) != angle_point_levels:
        raise ValueError("angle_point_levels must be strictly increasing and unique.")

    cache: dict[tuple[str, int], tuple[np.ndarray, dict[str, float]]] = {}
    checkpoint_dir = out_dir / ".convergence_cache"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def get_sweep(mesh: MeshConfig, points: int, label: str):
        key = (mesh.name, points)
        if key not in cache:
            checkpoint = torque_sweep_checkpoint_path(
                checkpoint_dir, machine, mesh, iq_test_a, points
            )
            if checkpoint.exists():
                table = load_torque_sweep_checkpoint(checkpoint, points)
                cache[key] = (table, summarize_torque_sweep(table))
                print(f"  Reusing convergence checkpoint: {checkpoint.name}")
            else:
                table, summary = run_torque_sweep(
                    machine,
                    mesh,
                    iq_test_a=iq_test_a,
                    num_points=points,
                    progress_label=label,
                )
                temporary_checkpoint = checkpoint.with_suffix(".tmp")
                np.savetxt(
                    temporary_checkpoint,
                    table,
                    delimiter=",",
                    header=(
                        "rotor_angle_deg,cogging_torque_nm,"
                        "loaded_total_torque_nm"
                    ),
                    comments="",
                )
                temporary_checkpoint.replace(checkpoint)
                cache[key] = (table, summary)
        return cache[key]

    studies: list[dict[str, float | str]] = []
    raw_rows: list[tuple[str, str, int, float, float, float]] = []
    for mesh in AIRGAP_MESH_LEVELS:
        table, summary = get_sweep(
            mesh, mesh_angle_points, f"Mesh convergence ({mesh.name})"
        )
        studies.append({
            "study": "mesh",
            "level": mesh.name,
            "points": mesh_angle_points,
            "airgap_radial_layers": mesh.airgap_radial_layers,
            "airgap_arc_deg": mesh.airgap_arc_deg,
            **summary,
        })
        raw_rows.extend(
            ("mesh", mesh.name, mesh_angle_points, *row)
            for row in table
        )

    converged_gap_mesh = ANGLE_INTEGRAL_MESH
    for points in angle_point_levels:
        level = f"n{points}"
        table, summary = get_sweep(
            converged_gap_mesh, points, f"Angle convergence ({points} points)"
        )
        studies.append({
            "study": "angle",
            "level": level,
            "points": points,
            "airgap_radial_layers": converged_gap_mesh.airgap_radial_layers,
            "airgap_arc_deg": converged_gap_mesh.airgap_arc_deg,
            **summary,
        })
        raw_rows.extend(
            ("angle", level, points, *row)
            for row in table
        )

    mesh_results = [row for row in studies if row["study"] == "mesh"]
    angle_results = [row for row in studies if row["study"] == "angle"]
    mesh_reference = mesh_results[-1]
    angle_reference = angle_results[-1]
    for rows, reference in ((mesh_results, mesh_reference), (angle_results, angle_reference)):
        for row in rows:
            row["loaded_mean_error"] = _relative_change(
                float(row["loaded_mean_nm"]), float(reference["loaded_mean_nm"])
            )
            row["cogging_pp_error"] = _relative_change(
                float(row["cogging_pp_nm"]), float(reference["cogging_pp_nm"])
            )

    out_dir.mkdir(parents=True, exist_ok=True)
    integral_acceptance, integral_diagnostics = run_torque_integral_convergence(
        machine,
        out_dir,
        iq_test_a=iq_test_a,
        mesh=converged_gap_mesh,
    )

    mesh_delta = _relative_change(
        float(mesh_results[-2]["loaded_mean_nm"]),
        float(mesh_results[-1]["loaded_mean_nm"]),
    )
    angle_delta = _relative_change(
        float(angle_results[-2]["loaded_mean_nm"]),
        float(angle_results[-1]["loaded_mean_nm"]),
    )
    mesh_cogging_delta = _relative_change(
        float(mesh_results[-2]["cogging_pp_nm"]),
        float(mesh_results[-1]["cogging_pp_nm"]),
    )
    angle_cogging_delta = _relative_change(
        float(angle_results[-2]["cogging_pp_nm"]),
        float(angle_results[-1]["cogging_pp_nm"]),
    )
    loaded_signal = abs(float(mesh_reference["loaded_mean_nm"]))
    cogging_pp = abs(float(mesh_reference["cogging_pp_nm"]))
    remesh_bias = abs(float(mesh_reference["cogging_mean_nm"]))
    absolute_mesh_delta = abs(
        float(mesh_results[-1]["loaded_mean_nm"])
        - float(mesh_results[-2]["loaded_mean_nm"])
    )
    absolute_angle_delta = abs(
        float(angle_results[-1]["loaded_mean_nm"])
        - float(angle_results[-2]["loaded_mean_nm"])
    )
    numerical_floor = max(
        remesh_bias,
        absolute_mesh_delta,
        absolute_angle_delta,
        integral_diagnostics["quadrature_absolute_delta_nm"],
        integral_diagnostics["radius_spread_nm"],
        1e-12,
    )
    load_to_cogging = loaded_signal / max(cogging_pp, 1e-12)
    load_to_numerical_floor = loaded_signal / numerical_floor

    acceptance = {
        "mesh_loaded_torque_relative_change": {
            "value": mesh_delta, "criterion": "<=", "limit": 0.02,
            "passed": mesh_delta <= 0.02,
        },
        "angle_loaded_torque_relative_change": {
            "value": angle_delta, "criterion": "<=", "limit": 0.02,
            "passed": angle_delta <= 0.02,
        },
        "mesh_cogging_pp_relative_change": {
            "value": mesh_cogging_delta, "criterion": "<=", "limit": 0.10,
            "passed": mesh_cogging_delta <= 0.10,
        },
        "angle_cogging_pp_relative_change": {
            "value": angle_cogging_delta, "criterion": "<=", "limit": 0.10,
            "passed": angle_cogging_delta <= 0.10,
        },
        "load_to_cogging_ratio": {
            "value": load_to_cogging, "criterion": ">=", "limit": 5.0,
            "passed": load_to_cogging >= 5.0,
        },
        "load_to_numerical_floor_ratio": {
            "value": load_to_numerical_floor, "criterion": ">=", "limit": 20.0,
            "passed": load_to_numerical_floor >= 20.0,
        },
        **integral_acceptance,
    }

    with (out_dir / "torque_convergence_samples.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow((
            "study", "level", "points", "rotor_angle_deg",
            "cogging_torque_nm", "loaded_total_torque_nm",
        ))
        writer.writerows(raw_rows)
    with (out_dir / "torque_convergence_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        fieldnames = (
            "study", "level", "points", "airgap_radial_layers", "airgap_arc_deg",
            "loaded_mean_nm", "loaded_ripple_pp_nm", "cogging_mean_nm",
            "cogging_pp_nm", "loaded_mean_error", "cogging_pp_error",
        )
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(studies)
    with (out_dir / "torque_convergence_acceptance.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(("metric", "value", "criterion", "limit", "passed"))
        for name, metric in acceptance.items():
            writer.writerow((
                name, metric["value"], metric["criterion"], metric["limit"], metric["passed"]
            ))

    for name, metric in acceptance.items():
        status = "PASS" if metric["passed"] else "FAIL"
        print(f"  {status:4s} {name}: {float(metric['value']):.6g}")
    return acceptance


# ══════════════════════════════════════════════════════════════════
# 6) EFFICIENCY MAP & PARALLEL WORKERS
# ══════════════════════════════════════════════════════════════════

def _worker_emap(args):
    """Multiprocessing worker function for pure performance sweeping."""
    machine, loss, rpm, iq, pts, mesh = args
    import pythoncom
    import time, random
    pythoncom.CoInitialize()
    time.sleep(random.uniform(0.2, 2.0)) # stagger COM startup to prevent hangs
    try:
        import femm
        femm.openfemm(1)
        cfg = SweepCfg(points_per_electrical_cycle=pts)
        rows = run_one_case(machine, loss, cfg, rpm, iq, progress=False, mesh=mesh)
        t_avg = float(np.mean(rows[:, 5]))
        ploss_avg = float(np.mean(rows[:, 9]))
        omega = rpm * 2 * math.pi / 60.0
        p_mech = t_avg * omega
        eta = (
            p_mech / (p_mech + ploss_avg)
            if p_mech >= 0.0 and (p_mech + ploss_avg) > 0.0
            else math.nan
        )
        return (rpm, iq, t_avg, ploss_avg, eta)
    except Exception as exc:
        traceback.print_exc()
        raise exc
    finally:
        try:
            femm.closefemm()
        except:
            pass


def _worker_basic(args):
    """Multiprocessing worker function for basic waveform generation."""
    (
        machine,
        loss,
        rpm,
        iq,
        pts,
        out_dir,
        show,
        analyses,
        mesh,
        snapshot_angle_deg,
        field_radial_points,
        field_angular_points,
        airgap_points,
    ) = args
    import pythoncom
    import time, random
    pythoncom.CoInitialize()
    time.sleep(random.uniform(0.2, 2.0)) # stagger COM startup to prevent hangs
    try:
        import femm
        femm.openfemm(1)
        cfg = SweepCfg(points_per_electrical_cycle=pts)
        case_tag = f"rpm{case_value_token(rpm)}_iq{case_value_token(iq)}"
        summary = None
        if "basic" in analyses:
            rows = run_one_case(
                machine, loss, cfg, rpm, iq, progress=False, mesh=mesh
            )
            save_case_outputs(rows, out_dir, case_tag, show)
            summary = (
                rpm,
                iq,
                float(np.mean(rows[:, 5])),
                float(np.mean(rows[:, 9])),
                float(rows[-1, 10]),
            )
        if analyses & {"field", "airgap"}:
            run_field_snapshot(
                machine,
                out_dir,
                rpm=rpm,
                iq_peak_a=iq,
                rotor_angle_deg=snapshot_angle_deg,
                analyses=analyses,
                mesh=mesh,
                show=show,
                field_radial_points=field_radial_points,
                field_angular_points=field_angular_points,
                airgap_points=airgap_points,
            )
        return summary
    except Exception as exc:
        traceback.print_exc()
        raise exc
    finally:
        try:
            femm.closefemm()
        except:
            pass

def run_efficiency_map(machine, loss, rpm_list, iq_list, out_dir,
                       points_per_cycle=12, show=False, workers=1,
                       mesh: MeshConfig | None = None):
    """Sweep rpm × iq, plot η contour."""
    mesh = MeshConfig() if mesh is None else mesh
    results = []   
    tasks = [
        (machine, loss, rpm, iq, points_per_cycle, mesh)
        for rpm in rpm_list
        for iq in iq_list
    ]

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_req = {executor.submit(_worker_emap, t): t for t in tasks}
            for future in tqdm(as_completed(future_to_req), total=len(tasks), desc="Emap Sweep"):
                results.append(future.result())
    else:
        # Sequential
        for t in tqdm(tasks, desc="Emap Sweep"):
            rpm, iq = t[2], t[3]
            cfg = SweepCfg(points_per_electrical_cycle=t[4])
            rows = run_one_case(
                machine, loss, cfg, rpm, iq, progress=True, mesh=mesh
            )
            t_avg = float(np.mean(rows[:, 5]))
            ploss_avg = float(np.mean(rows[:, 9]))
            omega = rpm * 2 * math.pi / 60.0
            p_mech = t_avg * omega
            eta = (
                p_mech / (p_mech + ploss_avg)
                if p_mech >= 0.0 and (p_mech + ploss_avg) > 0.0
                else math.nan
            )
            results.append((rpm, iq, t_avg, ploss_avg, eta))

    data = np.array(results, dtype=float)
    out_dir.mkdir(parents=True, exist_ok=True)
    efficiency_name = (
        "efficiency_estimate"
        if loss.calibrated
        else "copper_only_efficiency_upper_bound"
    )
    np.savetxt(
        out_dir / "efficiency_map.csv",
        data,
        delimiter=",",
        header=f"rpm,iq_peak_a,avg_torque_nm,avg_loss_w,{efficiency_name}",
        comments="",
    )

    # Match mesh shapes
    rpm_uq = np.unique(data[:, 0])
    iq_uq = np.unique(data[:, 1])
    n_rpm, n_iq = len(rpm_uq), len(iq_uq)
    
    # Sort data for meshgrid
    data = data[np.lexsort((data[:, 1], data[:, 0]))]
    RPM = data[:, 0].reshape(n_rpm, n_iq)
    TQ  = data[:, 2].reshape(n_rpm, n_iq)
    ETA = data[:, 4].reshape(n_rpm, n_iq) * 100 

    fig, ax = plt.subplots(figsize=(10, 6))
    levels = np.arange(0, 101, 5)
    cs = ax.contourf(RPM, TQ, ETA, levels=levels, cmap="RdYlGn")
    colorbar_label = (
        "Estimated efficiency [%]"
        if loss.calibrated
        else "Copper-only efficiency upper bound [%]"
    )
    fig.colorbar(cs, ax=ax, label=colorbar_label)
    ax.contour(RPM, TQ, ETA, levels=levels, colors="k", linewidths=0.3)
    ax.set_xlabel("Speed [rpm]"); ax.set_ylabel("Avg Torque [N·m]")
    ax.set_title(
        "Efficiency Map" if loss.calibrated else "Copper-only Efficiency Upper Bound"
    )
    ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(out_dir / "efficiency_map.png", dpi=150)
    plt.close(fig) if not show else plt.show()
    return data


# ══════════════════════════════════════════════════════════════════
# 7) TORQUE‑SPEED CURVE
# ══════════════════════════════════════════════════════════════════

def plot_torque_speed_curve(
    machine,
    drive,
    emap_data,
    out_dir,
    show=False,
    shared_config: SharedMotorConfig | None = None,
):
    """Plot raw FEA points against the shared constraint-aware dq envelope."""
    if shared_config is None:
        shared_config = load_motor_config()
    out_dir.mkdir(parents=True, exist_ok=True)
    rpm_vals = np.unique(emap_data[:, 0])
    t_max_fem = []
    for rpm in rpm_vals:
        mask = emap_data[:, 0] == rpm
        t_max_fem.append(np.max(emap_data[mask, 2]))
    t_max_fem = np.array(t_max_fem)

    from design_pmsm import evaluate_operating_point

    rpm_dense = np.linspace(0.0, rpm_vals.max() * 1.5, 200)
    envelope_points = [
        evaluate_operating_point(
            shared_config, float(rpm), 0.0, maximum_torque=True
        )
        for rpm in rpm_dense
    ]
    t_envelope = np.maximum(
        0.0, np.array([point.shaft_torque_nm for point in envelope_points])
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.plot(
        rpm_vals,
        t_max_fem,
        "bs-",
        label="Raw FEA max (not voltage filtered)",
        markersize=5,
    )
    ax1.plot(rpm_dense, t_envelope, "r--", label="Shared voltage/current envelope")
    ax1.set_ylabel("Torque [N·m]"); ax1.set_title("Torque–Speed Curve")
    ax1.legend(); ax1.grid(True)

    p_fem = t_max_fem * rpm_vals * 2 * math.pi / 60.0
    p_env = np.array(
        [max(0.0, point.shaft_power_w) for point in envelope_points]
    )
    ax2.plot(rpm_vals, p_fem, "bs-", label="Raw FEA max power", markersize=5)
    ax2.plot(rpm_dense, p_env, "r--", label="Shared constrained envelope")
    ax2.set_xlabel("Speed [rpm]"); ax2.set_ylabel("Power [W]")
    ax2.set_title("Power–Speed Curve")
    ax2.legend(); ax2.grid(True)

    fig.tight_layout()
    fig.savefig(out_dir / "torque_speed_curve.png", dpi=150)
    plt.close(fig) if not show else plt.show()


# ══════════════════════════════════════════════════════════════════
# CLI + MAIN
# ══════════════════════════════════════════════════════════════════

ALL_ANALYSES = {
    "basic", "field", "airgap", "cogging", "inductance", "validate",
    "meshcheck", "convergence", "emap", "tncurve",
}


def parse_args():
    p = argparse.ArgumentParser(
        description="FEMM 2D SPMSM template – comprehensive motor analysis.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--analysis", nargs="+", default=["basic"],
                   help=(
                       "Options: basic, field, airgap, cogging, inductance, validate, "
                       "meshcheck, convergence, emap, tncurve, all."
                   ))
    p.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="Shared SI-unit motor JSON used by both quick and FEMM models.",
    )
    p.add_argument("--rpm-list", nargs="+", default=["300", "600", "900"],
                   help="Speed list (rpm) for basic analysis.")
    p.add_argument("--iq-list", nargs="+", default=["5", "10", "15"],
                   help="Current peak list (A) for basic analysis.")
    p.add_argument("--points", type=int, default=72, help="Samples per electrical cycle.")
    p.add_argument("--out", type=Path, default=Path("output_femm"), help="Output folder.")
    p.add_argument("--show", action="store_true", help="Show matplotlib figures.")
    p.add_argument("--show-femm", action="store_true", help="Show FEMM GUI.")
    p.add_argument("--workers", type=int, default=1, help="Number of parallel FEMM workers (default 1).")
    p.add_argument(
        "--mesh-level",
        choices=[mesh.name for mesh in MESH_LEVELS],
        default="medium",
        help="Explicit FEMM mesh level for normal analyses and validation.",
    )
    p.add_argument(
        "--snapshot-angle-deg",
        type=float,
        default=0.0,
        help="Explicit mechanical rotor angle for field/air-gap snapshots.",
    )
    p.add_argument(
        "--field-radial-points",
        type=int,
        default=48,
        help="Radial samples in a field-map snapshot (minimum 8).",
    )
    p.add_argument(
        "--field-angular-points",
        type=int,
        default=180,
        help="Angular samples in a field-map snapshot (minimum 36).",
    )
    p.add_argument(
        "--airgap-points",
        type=int,
        default=360,
        help="Angular samples around the air-gap midline (minimum 36).",
    )
    
    p.add_argument("--emap-rpm", nargs="+", default=["300,600,900,1200"])
    p.add_argument("--emap-iq", nargs="+", default=["2.5,5,10,15"])
    p.add_argument("--emap-pts", type=int, default=12)
    p.add_argument("--cogging-steps", type=int, default=72)
    p.add_argument("--ind-steps", type=int, default=37)
    p.add_argument("--ind-current", type=float, default=1.0)
    p.add_argument("--validation-iq", type=float, default=5.0)
    p.add_argument("--validation-delta-current", type=float, default=1.0)
    p.add_argument("--validation-steps", type=int, default=12)
    p.add_argument(
        "--mesh-check-reference-level",
        choices=[mesh.name for mesh in MESH_LEVELS],
        default="medium",
        help="Baseline global mesh used by the one-angle mesh comparison.",
    )
    p.add_argument(
        "--mesh-check-level",
        choices=[mesh.name for mesh in MESH_LEVELS],
        default="fine",
        help="Target global mesh used by the one-angle mesh comparison.",
    )
    p.add_argument("--convergence-iq", type=float, default=5.0)
    p.add_argument("--convergence-mesh-points", type=int, default=8)
    p.add_argument("--convergence-angle-points", nargs="+", default=["6,12,24"])
    return p.parse_args()


def open_femm_or_raise(show_femm: bool):
    try:
        femm.openfemm(0 if show_femm else 1)
    except Exception as exc:
        raise SystemExit(f"Failed to connect FEMM COM server. Ensure FEMM 4.2 is installed and registered. {exc}") from exc


def write_resolved_femm_config(
    out_dir: Path,
    source_path: Path,
    config: SharedMotorConfig,
    machine: Machine,
    loss: LossThermal,
    drive: Drive,
) -> Path:
    shared_dict = config_to_dict(config)
    physical_hash = physical_model_fingerprint(
        config.machine,
        config.materials,
        config.conventions,
    )
    payload = {
        "source_path": str(source_path.resolve()),
        "shared_config_sha256": canonical_config_sha256(config),
        "physical_model_sha256": physical_hash,
        "shared_config": shared_dict,
        "resolved_femm_machine": asdict(machine),
        "resolved_femm_loss_thermal": asdict(loss),
        "resolved_femm_drive": asdict(drive),
    }
    path = out_dir / "resolved_femm_config.json"
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _snapshot_output_files(out_dir: Path) -> dict[str, tuple[int, int]]:
    snapshot: dict[str, tuple[int, int]] = {}
    for path in out_dir.rglob("*"):
        if not path.is_file() or path.name == "femm_run_manifest.json":
            continue
        stat = path.stat()
        snapshot[path.relative_to(out_dir).as_posix()] = (
            stat.st_mtime_ns,
            stat.st_size,
        )
    return snapshot


def _jsonable_arg(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_jsonable_arg(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable_arg(item) for key, item in value.items()}
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def write_femm_run_manifest(
    out_dir: Path,
    *,
    before: dict[str, tuple[int, int]],
    started_at_utc: str,
    analyses: set[str],
    args: argparse.Namespace,
    config: SharedMotorConfig,
) -> Path:
    """Append provenance for files created or replaced by one successful run."""
    changed_files = []
    for path in sorted(out_dir.rglob("*"), key=lambda item: item.as_posix()):
        if not path.is_file() or path.name == "femm_run_manifest.json":
            continue
        relative = path.relative_to(out_dir).as_posix()
        stat = path.stat()
        state = (stat.st_mtime_ns, stat.st_size)
        if before.get(relative) == state:
            continue
        changed_files.append(
            {
                "path": relative,
                "sha256": _file_sha256(path),
                "size_bytes": stat.st_size,
            }
        )
    if not changed_files:
        raise RuntimeError("Successful FEMM run did not create or replace any output files.")

    manifest_path = out_dir / "femm_run_manifest.json"
    if manifest_path.is_file():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"Invalid FEMM run manifest {manifest_path}: {exc}") from exc
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            raise RuntimeError(f"Unsupported FEMM run manifest: {manifest_path}")
        runs = payload.get("runs")
        if not isinstance(runs, list):
            raise RuntimeError(f"Invalid FEMM run list in {manifest_path}")
    else:
        payload = {"schema_version": 1, "runs": []}
        runs = payload["runs"]

    runs.append(
        {
            "started_at_utc": started_at_utc,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "analyses": sorted(analyses),
            "arguments": {
                key: _jsonable_arg(value)
                for key, value in sorted(vars(args).items())
            },
            "command": [str(item) for item in sys.argv],
            "shared_config_sha256": canonical_config_sha256(config),
            "physical_model_sha256": physical_model_fingerprint(
                config.machine,
                config.materials,
                config.conventions,
            ),
            "files": changed_files,
        }
    )
    temporary = manifest_path.with_name(".femm_run_manifest.tmp.json")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(manifest_path)
    return manifest_path


def main():
    args = parse_args()
    # multiprocessing entry point guard handled naturally
    raw_analyses = []
    for item in args.analysis:
        raw_analyses.extend(
            value
            for value in re.split(r"[\s,]+", item.strip().lower())
            if value
        )
    unknown_analyses = set(raw_analyses) - ALL_ANALYSES - {"all"}
    if unknown_analyses:
        choices = ", ".join(sorted(ALL_ANALYSES | {"all"}))
        unknown = ", ".join(sorted(unknown_analyses))
        raise SystemExit(
            f"Unknown analysis name(s): {unknown}. Valid choices: {choices}."
        )
    if "all" in raw_analyses:
        analyses = ALL_ANALYSES.copy()
    else:
        analyses = set(raw_analyses)
    if not analyses:
        raise SystemExit("At least one analysis must be selected.")
    if "tncurve" in analyses:
        analyses.add("emap")

    rpm_list = parse_list(args.rpm_list)
    iq_list = parse_list(args.iq_list)
    try:
        shared_config = load_motor_config(args.config)
        machine = machine_from_motor_config(shared_config)
        loss = loss_from_motor_config(shared_config)
        drive = drive_from_motor_config(shared_config)
    except (ConfigError, ValueError) as exc:
        raise SystemExit(f"Shared motor configuration error: {exc}") from exc
    if args.points < 3:
        raise SystemExit("--points must be at least 3.")
    cfg = SweepCfg(points_per_electrical_cycle=args.points)
    mesh = mesh_config_by_name(args.mesh_level)
    mesh_check_reference = mesh_config_by_name(args.mesh_check_reference_level)
    mesh_check_target = mesh_config_by_name(args.mesh_check_level)
    if "meshcheck" in analyses and mesh_check_reference.name == mesh_check_target.name:
        raise SystemExit(
            "--mesh-check-reference-level and --mesh-check-level must differ."
        )
    convergence_angle_points = tuple(
        int(value) for value in parse_list(args.convergence_angle_points)
    )
    if args.field_radial_points < 8:
        raise SystemExit("--field-radial-points must be at least 8.")
    if args.field_angular_points < 36:
        raise SystemExit("--field-angular-points must be at least 36.")
    if args.airgap_points < 36:
        raise SystemExit("--airgap-points must be at least 36.")

    args.out.mkdir(parents=True, exist_ok=True)
    files_before_run = _snapshot_output_files(args.out)
    run_started_at_utc = datetime.now(timezone.utc).isoformat()
    resolved_path = write_resolved_femm_config(
        args.out, args.config, shared_config, machine, loss, drive
    )
    
    # We delay opening FEMM globally until we are sure we are not running purely in parallel
    print(f"Selected analyses: {', '.join(sorted(analyses))}")
    print(
        f"Motor config: {shared_config.name} "
        f"({machine.slots} slots / {2 * machine.pole_pairs} poles)"
    )
    print(f"Resolved config: {resolved_path}")
    print(f"Workers allocated: {args.workers}\n")
    requested_currents: list[float] = []
    if analyses & {"basic", "field", "airgap"}:
        requested_currents.extend(abs(value) for value in iq_list)
    if analyses & {"validate", "meshcheck"}:
        requested_currents.extend(
            (abs(args.validation_iq), abs(args.validation_delta_current))
        )
    if "convergence" in analyses:
        requested_currents.append(abs(args.convergence_iq))
    if "emap" in analyses:
        requested_currents.extend(abs(value) for value in parse_list(args.emap_iq))
    validated_current = shared_config.electromagnetic.validated_current_peak_a
    if requested_currents and max(requested_currents) > validated_current + 1e-12:
        print(
            "WARNING: requested current exceeds the electromagnetic parameter "
            f"validation point ({validated_current:g} A peak); results are extrapolated.\n"
        )
    if analyses & {"basic", "field", "airgap", "emap"}:
        print(
            "NOTE: FEMM rpm/Iq cases impose current directly and are not filtered "
            "by the inverter voltage ellipse; use the quick sweep for reachability.\n"
        )
    if not loss.calibrated and analyses & {"basic", "emap", "tncurve"}:
        print(
            "WARNING: additional losses are uncalibrated; efficiency/temperature "
            "outputs are screening estimates only.\n"
        )

    femm_is_open = False
    try:
        # Sequential/Parallel routing
        if analyses & {"basic", "field", "airgap"}:
            summary_rows = []
            
            if args.workers > 1:
                tasks = [
                    (
                        machine,
                        loss,
                        r,
                        i,
                        args.points,
                        args.out,
                        args.show,
                        analyses,
                        mesh,
                        args.snapshot_angle_deg,
                        args.field_radial_points,
                        args.field_angular_points,
                        args.airgap_points,
                    )
                    for r in rpm_list
                    for i in iq_list
                ]

                with ProcessPoolExecutor(max_workers=args.workers) as executor:
                    for future in tqdm(as_completed([executor.submit(_worker_basic, t) for t in tasks]), 
                                       total=len(tasks), desc="Basic/Field Sweep"):
                        summary = future.result()
                        if summary is not None:
                            summary_rows.append(summary)
            else:
                open_femm_or_raise(args.show_femm)
                femm_is_open = True
                for rpm in tqdm(rpm_list, desc="RPM List"):
                    for iq in tqdm(iq_list, desc="Iq List", leave=False):
                        if "basic" in analyses:
                            case_tag = (
                                f"rpm{case_value_token(rpm)}_"
                                f"iq{case_value_token(iq)}"
                            )
                            rows = run_one_case(
                                machine,
                                loss,
                                cfg,
                                rpm,
                                iq,
                                progress=True,
                                mesh=mesh,
                            )
                            save_case_outputs(rows, args.out, case_tag, args.show)
                            avg_torque = float(np.mean(rows[:, 5]))
                            avg_loss = float(np.mean(rows[:, 9]))
                            final_temp = float(rows[-1, 10])
                            summary_rows.append(
                                (rpm, iq, avg_torque, avg_loss, final_temp)
                            )
                        if analyses & {"field", "airgap"}:
                            run_field_snapshot(
                                machine,
                                args.out,
                                rpm=rpm,
                                iq_peak_a=iq,
                                rotor_angle_deg=args.snapshot_angle_deg,
                                analyses=analyses,
                                mesh=mesh,
                                show=args.show,
                                field_radial_points=args.field_radial_points,
                                field_angular_points=args.field_angular_points,
                                airgap_points=args.airgap_points,
                            )

            if "basic" in analyses and summary_rows:
                summary = np.array(summary_rows, dtype=float)
                np.savetxt(args.out / "femm_summary.csv", summary, delimiter=",",
                           header=(
                               "rpm,iq_peak_a,avg_torque_nm,avg_total_loss_w,"
                               "estimated_steady_winding_temp_c"
                           ),
                           comments="")
                print("\nBasic summary written to:", args.out / "femm_summary.csv")

        # Fallback to main FEMM instance for remaining sequential parts
        sequential_analyses = {
            "cogging", "inductance", "validate", "meshcheck", "convergence"
        }
        needs_main_femm = bool(analyses & sequential_analyses) or (
            args.workers == 1 and "emap" in analyses
        )
        if needs_main_femm and not femm_is_open:
            open_femm_or_raise(args.show_femm)
            femm_is_open = True
            
        if "cogging" in analyses:
            print("\n── Cogging torque analysis ──")
            run_cogging_torque(
                machine, args.out, num_steps=args.cogging_steps,
                show=args.show, mesh=mesh,
            )

        if "inductance" in analyses:
            print("\n── Inductance analysis ──")
            run_inductance_analysis(machine, args.out, test_current=args.ind_current,
                                    num_steps=args.ind_steps, show=args.show, mesh=mesh)

        if "validate" in analyses:
            print("\n── Electromagnetic acceptance ──")
            run_electromagnetic_acceptance(
                machine,
                args.out,
                iq_test_a=args.validation_iq,
                delta_current_a=args.validation_delta_current,
                num_steps=args.validation_steps,
                mesh=mesh,
            )

        if "meshcheck" in analyses:
            print("\n── Global mesh sensitivity check ──")
            reference_report = run_electromagnetic_mesh_check(
                machine,
                args.out,
                iq_test_a=args.validation_iq,
                delta_current_a=args.validation_delta_current,
                mesh=mesh_check_reference,
            )
            target_report = run_electromagnetic_mesh_check(
                machine,
                args.out,
                iq_test_a=args.validation_iq,
                delta_current_a=args.validation_delta_current,
                mesh=mesh_check_target,
            )
            comparison = compare_electromagnetic_mesh_checks(
                reference_report,
                target_report,
                args.out,
            )
            status = "PASS" if comparison["all_passed"] else "FAIL"
            print(
                f"  {status} {mesh_check_reference.name}/"
                f"{mesh_check_target.name} global-mesh comparison"
            )

        if "convergence" in analyses:
            print("\n── Torque convergence ──")
            run_torque_convergence(
                machine,
                args.out,
                iq_test_a=args.convergence_iq,
                mesh_angle_points=args.convergence_mesh_points,
                angle_point_levels=convergence_angle_points,
            )

        emap_data = None
        if "emap" in analyses:
            print("\n── Efficiency map ──")
            emap_rpm = parse_list(args.emap_rpm)
            emap_iq = parse_list(args.emap_iq)
            emap_data = run_efficiency_map(machine, loss, emap_rpm, emap_iq,
                                           args.out, points_per_cycle=args.emap_pts,
                                           show=args.show, workers=args.workers,
                                           mesh=mesh)

        if "tncurve" in analyses and emap_data is not None:
            print("\n── Torque‑speed curve ──")
            plot_torque_speed_curve(
                machine,
                drive,
                emap_data,
                args.out,
                args.show,
                shared_config=shared_config,
            )

    finally:
        try:
            femm.closefemm()
        except Exception:
            pass

    manifest_path = write_femm_run_manifest(
        args.out,
        before=files_before_run,
        started_at_utc=run_started_at_utc,
        analyses=analyses,
        args=args,
        config=shared_config,
    )
    print(f"Run manifest: {manifest_path}")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
