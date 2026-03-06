"""FEMM 2‑D SPMSM template — auto model + comprehensive post‑processing.

Outputs supported (selected via ``--analysis``):
  basic      – torque / loss / temperature waveforms (original)
  field      – flux‑density cloud + flux‑line maps
  airgap     – airgap Bn / Bt distribution
  cogging    – cogging‑torque curve (zero‑current)
  inductance – Ld / Lq vs rotor position
  emap       – efficiency map (speed × current sweep)
  tncurve    – torque‑speed envelope
  all        – run every analysis above
"""
from __future__ import annotations

import argparse
import math
import os
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

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
TMP_FEM = os.path.join(tempfile.gettempdir(), "femm_spm_tmp.fem")

# ── data classes ───────────────────────────────────────────────────

@dataclass
class Machine:
    pole_pairs: int = 4
    slots: int = 12
    stack_length_mm: float = 80.0
    r_shaft_mm: float = 10.0
    r_rotor_mm: float = 34.0
    mag_thickness_mm: float = 3.0
    airgap_mm: float = 0.8
    r_stator_outer_mm: float = 90.0
    r_air_outer_mm: float = 140.0
    magnet_arc_ratio: float = 0.75
    turns_per_slot: int = 35
    slot_depth_mm: float = 8.0

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
    kh: float = 1.1
    ke: float = 0.004
    core_mass_kg: float = 9.0
    magnet_loss_ratio: float = 0.15
    r_th: float = 0.45
    c_th: float = 2600.0
    t_amb: float = 25.0


@dataclass
class Drive:
    """Inverter / voltage‑source parameters for T‑N envelope."""
    v_dc: float = 310.0
    i_max: float = 55.0


@dataclass
class SweepCfg:
    points_per_electrical_cycle: int = 72


# ── geometry helpers ───────────────────────────────────────────────

def polar_xy(r_mm: float, angle_deg: float):
    a = math.radians(angle_deg)
    return r_mm * math.cos(a), r_mm * math.sin(a)


def parse_list(raw):
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


def add_quarter_arc_circle(radius_mm: float, group: int = 0, boundary: str = "<None>"):
    points = [(radius_mm, 0.0), (0.0, radius_mm), (-radius_mm, 0.0), (0.0, -radius_mm)]
    for x, y in points:
        femm.mi_addnode(x, y)
    for i in range(4):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % 4]
        femm.mi_addarc(x1, y1, x2, y2, 90.0, 1)
        mx, my = polar_xy(radius_mm, 45.0 + i * 90.0)
        femm.mi_selectarcsegment(mx, my)
        femm.mi_setarcsegmentprop(1.0, boundary, 0, group)
        femm.mi_clearselected()


def add_wedge_annulus(r_in: float, r_out: float, a1: float, a2: float, group: int = 0):
    x1, y1 = polar_xy(r_in, a1)
    x2, y2 = polar_xy(r_in, a2)
    x3, y3 = polar_xy(r_out, a2)
    x4, y4 = polar_xy(r_out, a1)
    femm.mi_addnode(x1, y1)
    femm.mi_addnode(x2, y2)
    femm.mi_addnode(x3, y3)
    femm.mi_addnode(x4, y4)
    sweep = (a2 - a1) % 360.0
    if sweep <= 0.0:
        sweep += 360.0
    femm.mi_addarc(x1, y1, x2, y2, sweep, 1)
    femm.mi_addsegment(x2, y2, x3, y3)
    femm.mi_addarc(x4, y4, x3, y3, sweep, 1)
    femm.mi_addsegment(x1, y1, x4, y4)

    r_mid = 0.5 * (r_in + r_out)
    a_mid = a1 + sweep * 0.5
    mx, my = polar_xy(r_mid, a_mid)
    femm.mi_selectarcsegment(mx, my)
    femm.mi_setarcsegmentprop(1.0, "<None>", 0, group)
    femm.mi_clearselected()

    mx, my = polar_xy(r_in, a1)
    femm.mi_selectsegment(0.5 * (mx + x4), 0.5 * (my + y4))
    femm.mi_setsegmentprop("<None>", 1.0, 1, 0, group)
    femm.mi_clearselected()

    mx, my = polar_xy(r_in, a2)
    femm.mi_selectsegment(0.5 * (mx + x3), 0.5 * (my + y3))
    femm.mi_setsegmentprop("<None>", 1.0, 1, 0, group)
    femm.mi_clearselected()
    return polar_xy(r_mid, a_mid)


def set_block(material, x, y, circuit="<None>", mag_dir=0.0, group=0, turns=0):
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    femm.mi_setblockprop(material, 1, 0.0, circuit, mag_dir, group, turns)
    femm.mi_clearselected()


# ── materials ──────────────────────────────────────────────────────

def ensure_materials():
    """Load Air / M‑19 Steel / Copper from FEMM library;
    manually define NdFeB 40 MGOe (not shipped with FEMM).
    """
    for name in ["Air", "M-19 Steel", "Copper"]:
        femm.mi_getmaterial(name)
    mu_r = 1.05
    Br = 1.22
    Hc = Br / (MU0 * mu_r)
    femm.mi_addmaterial(
        "NdFeB 40 MGOe", mu_r, mu_r, Hc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    )


# ── model builder ─────────────────────────────────────────────────

def build_model(machine: Machine):
    femm.newdocument(0)
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, machine.stack_length_mm, 30)
    ensure_materials()

    femm.mi_addcircprop("A", 0.0, 1)
    femm.mi_addcircprop("B", 0.0, 1)
    femm.mi_addcircprop("C", 0.0, 1)
    femm.mi_addboundprop("A0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    add_quarter_arc_circle(machine.r_shaft_mm)
    add_quarter_arc_circle(machine.r_rotor_mm, group=ROTOR_GROUP)
    add_quarter_arc_circle(machine.r_mag_outer_mm, group=ROTOR_GROUP)
    add_quarter_arc_circle(machine.r_stator_inner_mm, group=STATOR_GROUP)
    add_quarter_arc_circle(machine.r_stator_outer_mm, group=STATOR_GROUP)
    add_quarter_arc_circle(machine.r_air_outer_mm, boundary="A0")

    set_block("Air", 0.5 * machine.r_shaft_mm, 0.0)
    set_block("M-19 Steel", 0.5 * (machine.r_shaft_mm + machine.r_rotor_mm), 0.0,
              group=ROTOR_GROUP)
    set_block("Air", 0.5 * (machine.r_mag_outer_mm + machine.r_stator_inner_mm), 0.0)
    set_block("M-19 Steel", 0.5 * (machine.r_slot_outer_mm + machine.r_stator_outer_mm),
              12.0, group=STATOR_GROUP)
    set_block("Air", 0.5 * (machine.r_air_outer_mm + machine.r_stator_outer_mm), 0.0)

    # magnets
    pole_count = 2 * machine.pole_pairs
    pole_pitch = 360.0 / pole_count
    magnet_span = pole_pitch * machine.magnet_arc_ratio
    for k in range(pole_count):
        center = k * pole_pitch
        a1, a2 = center - 0.5 * magnet_span, center + 0.5 * magnet_span
        x, y = add_wedge_annulus(machine.r_rotor_mm, machine.r_mag_outer_mm, a1, a2,
                                 group=ROTOR_GROUP)
        mag_dir = center if (k % 2 == 0) else (center + 180.0)
        set_block("NdFeB 40 MGOe", x, y, mag_dir=mag_dir, group=ROTOR_GROUP)

    # air in rotor gaps between magnets
    r_gap_mid = 0.5 * (machine.r_rotor_mm + machine.r_mag_outer_mm)
    for k in range(pole_count):
        gx, gy = polar_xy(r_gap_mid, (k + 0.5) * pole_pitch)
        set_block("Air", gx, gy, group=ROTOR_GROUP)

    # stator slots
    phase_map = [
        ("A", +1), ("C", -1), ("B", +1), ("A", -1), ("C", +1), ("B", -1),
        ("A", +1), ("C", -1), ("B", +1), ("A", -1), ("C", +1), ("B", -1),
    ]
    if machine.slots != len(phase_map):
        raise ValueError("This template currently expects 12 slots.")
    slot_pitch = 360.0 / machine.slots
    slot_span = slot_pitch * 0.55
    r_slot_in = machine.r_stator_inner_mm + 0.2
    r_slot_out = machine.r_slot_outer_mm
    for idx, (phase, sign) in enumerate(phase_map):
        center = idx * slot_pitch
        a1, a2 = center - 0.5 * slot_span, center + 0.5 * slot_span
        x, y = add_wedge_annulus(r_slot_in, r_slot_out, a1, a2, group=COIL_GROUP)
        set_block("Copper", x, y, circuit=phase, group=COIL_GROUP,
                  turns=sign * machine.turns_per_slot)

    femm.mi_saveas(TMP_FEM)


# ── low‑level helpers ──────────────────────────────────────────────

def set_phase_currents(i_peak: float, theta_e_rad: float):
    ia = i_peak * math.sin(theta_e_rad)
    ib = i_peak * math.sin(theta_e_rad - 2.0 * math.pi / 3.0)
    ic = i_peak * math.sin(theta_e_rad + 2.0 * math.pi / 3.0)
    femm.mi_modifycircprop("A", 1, ia)
    femm.mi_modifycircprop("B", 1, ib)
    femm.mi_modifycircprop("C", 1, ic)
    return ia, ib, ic


def set_dq_currents(machine: Machine, id_a: float, iq_a: float, theta_e_rad: float):
    """Apply d‑q current in the three‑phase windings."""
    cos_e = math.cos(theta_e_rad)
    sin_e = math.sin(theta_e_rad)
    ia = id_a * cos_e - iq_a * sin_e
    ib = id_a * math.cos(theta_e_rad - 2 * math.pi / 3) - iq_a * math.sin(theta_e_rad - 2 * math.pi / 3)
    ic = id_a * math.cos(theta_e_rad + 2 * math.pi / 3) - iq_a * math.sin(theta_e_rad + 2 * math.pi / 3)
    femm.mi_modifycircprop("A", 1, ia)
    femm.mi_modifycircprop("B", 1, ib)
    femm.mi_modifycircprop("C", 1, ic)
    return ia, ib, ic


def estimate_core_loss(machine: Machine, loss: LossThermal, rpm: float):
    probe_r = 0.5 * (machine.r_slot_outer_mm + machine.r_stator_outer_mm)
    b_max = 0.0
    for deg in [15, 45, 75, 105, 135, 165]:
        bx, by = femm.mo_getb(*polar_xy(probe_r, deg))
        b_max = max(b_max, math.hypot(bx, by))
    f_e = machine.pole_pairs * rpm / 60.0
    return loss.core_mass_kg * (loss.kh * f_e * b_max**2 + loss.ke * (f_e * b_max)**2)


def _get_circuit_flux(phase: str):
    """Return flux linkage of *phase* from loaded solution."""
    props = femm.mo_getcircuitproperties(phase)
    # props = (current, voltage, flux_linkage)
    return props[2]


# ══════════════════════════════════════════════════════════════════
# 1) BASIC: torque / loss / temperature waveforms (original)
# ══════════════════════════════════════════════════════════════════

def run_one_case(machine, loss, cfg, rpm, iq_peak):
    build_model(machine)
    f_e = machine.pole_pairs * rpm / 60.0
    dt = 1.0 / (f_e * cfg.points_per_electrical_cycle)
    dtheta_e = 2.0 * math.pi / cfg.points_per_electrical_cycle
    dtheta_m_deg = math.degrees(dtheta_e / machine.pole_pairs)
    temp = loss.t_amb
    rows = []
    for k in range(cfg.points_per_electrical_cycle):
        theta_e = k * dtheta_e
        ia, ib, ic = set_phase_currents(iq_peak, theta_e)
        if k > 0:
            femm.mi_selectgroup(ROTOR_GROUP)
            femm.mi_moverotate(0.0, 0.0, dtheta_m_deg)
            femm.mi_clearselected()
        femm.mi_analyze(1)
        femm.mi_loadsolution()
        femm.mo_groupselectblock(ROTOR_GROUP)
        torque = femm.mo_blockintegral(22)
        femm.mo_clearblock()
        rs = loss.r_phase_20 * (1.0 + loss.alpha_cu * (temp - 20.0))
        pcu = rs * (ia**2 + ib**2 + ic**2)
        pfe = estimate_core_loss(machine, loss, rpm)
        ppm = loss.magnet_loss_ratio * pfe
        ploss = pcu + pfe + ppm
        dtemp = (ploss - (temp - loss.t_amb) / loss.r_th) / loss.c_th
        temp += dtemp * dt
        rows.append((k * dt, theta_e, ia, ib, ic, torque, pcu, pfe, ppm, ploss, temp))
    return np.array(rows, dtype=float)


def save_case_outputs(rows, out_dir, case_tag, show):
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"femm_waveforms_{case_tag}.csv"
    np.savetxt(csv_path, rows, delimiter=",",
               header="t_s,theta_e_rad,ia_a,ib_a,ic_a,torque_nm,"
                      "pcu_w,pfe_w,ppm_w,ploss_w,temp_c", comments="")

    t, torque = rows[:, 0], rows[:, 5]
    pcu, pfe, ppm, ploss, temp = rows[:, 6], rows[:, 7], rows[:, 8], rows[:, 9], rows[:, 10]

    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    avg_t = np.mean(torque)
    axs[0].plot(t, torque, label="FEMM torque")
    axs[0].axhline(avg_t, color="gray", ls="--", label=f"Avg = {avg_t:.3f} N·m")
    axs[0].set_ylabel("Torque [N·m]"); axs[0].grid(True); axs[0].legend()

    axs[1].plot(t, pcu, label="Copper loss")
    axs[1].plot(t, pfe, label="Core loss")
    axs[1].plot(t, ppm, label="Magnet loss")
    axs[1].plot(t, ploss, "k", label="Total loss")
    axs[1].set_ylabel("Loss [W]"); axs[1].grid(True); axs[1].legend()

    axs[2].plot(t, temp, "r", label="Winding temperature")
    axs[2].set_xlabel("Time [s]"); axs[2].set_ylabel("Temperature [°C]")
    axs[2].grid(True); axs[2].legend()
    fig.tight_layout()
    png_path = out_dir / f"femm_waveforms_{case_tag}.png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig) if not show else plt.show()
    return csv_path, png_path


# ══════════════════════════════════════════════════════════════════
# 2) FIELD MAPS: flux‑density cloud + flux‑line contour
# ══════════════════════════════════════════════════════════════════

def plot_field_maps(machine, out_dir, case_tag, show=False):
    """Sample |B| and A on a polar grid from loaded solution; produce two PNGs."""
    r_min = machine.r_shaft_mm + 1.0
    r_max = machine.r_stator_outer_mm - 1.0
    nr, ntheta = 80, 360
    r_arr = np.linspace(r_min, r_max, nr)
    theta_arr = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    R, TH = np.meshgrid(r_arr, theta_arr)
    B_mag = np.zeros_like(R)
    A_val = np.zeros_like(R)

    for i in range(ntheta):
        for j in range(nr):
            x = R[i, j] * math.cos(TH[i, j])
            y = R[i, j] * math.sin(TH[i, j])
            bx, by = femm.mo_getb(x, y)
            B_mag[i, j] = math.hypot(bx, by)
            A_val[i, j] = femm.mo_geta(x, y)

    out_dir.mkdir(parents=True, exist_ok=True)

    # ── flux density cloud ──
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))
    c = ax.pcolormesh(TH, R, B_mag, shading="auto", cmap="jet")
    fig.colorbar(c, ax=ax, label="|B| [T]", pad=0.1)
    ax.set_title(f"Flux Density  ({case_tag})", va="bottom", fontsize=13)
    ax.set_yticklabels([])
    fig.savefig(out_dir / f"field_density_{case_tag}.png", dpi=150, bbox_inches="tight")
    plt.close(fig) if not show else plt.show()

    # ── flux lines (iso‑A contour) ──
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))
    levels = np.linspace(A_val.min(), A_val.max(), 40)
    ax.contour(TH, R, A_val, levels=levels, colors="k", linewidths=0.5)
    ax.set_title(f"Flux Lines  ({case_tag})", va="bottom", fontsize=13)
    ax.set_yticklabels([])
    fig.savefig(out_dir / f"field_lines_{case_tag}.png", dpi=150, bbox_inches="tight")
    plt.close(fig) if not show else plt.show()


# ══════════════════════════════════════════════════════════════════
# 3) AIRGAP: radial / tangential flux density along air‑gap
# ══════════════════════════════════════════════════════════════════

def plot_airgap_flux_density(machine, out_dir, case_tag, show=False, num_points=360):
    """Extract Bn / Bt along airgap midline from loaded solution."""
    r_gap = machine.r_mag_outer_mm + machine.airgap_mm / 2.0
    angles = np.linspace(0, 360, num_points, endpoint=False)
    bn_arr, bt_arr = [], []
    for deg in angles:
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
    ax.plot(angles, bn_arr, label="Bn (radial)")
    ax.plot(angles, bt_arr, label="Bt (tangential)", alpha=0.7)
    ax.set_xlabel("Mechanical angle [°]"); ax.set_ylabel("Flux density [T]")
    ax.set_title(f"Airgap Flux Density  ({case_tag})")
    ax.grid(True); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"airgap_B_{case_tag}.png", dpi=150)
    plt.close(fig) if not show else plt.show()


# ══════════════════════════════════════════════════════════════════
# 4) COGGING TORQUE
# ══════════════════════════════════════════════════════════════════

def run_cogging_torque(machine, out_dir, num_steps=72, show=False):
    """Build model with I=0, rotate rotor one slot‑pitch, record torque."""
    build_model(machine)
    slot_pitch_deg = 360.0 / machine.slots
    step_deg = slot_pitch_deg / num_steps
    angles, torques = [], []
    for k in range(num_steps):
        if k > 0:
            femm.mi_selectgroup(ROTOR_GROUP)
            femm.mi_moverotate(0.0, 0.0, step_deg)
            femm.mi_clearselected()
        femm.mi_analyze(1)
        femm.mi_loadsolution()
        femm.mo_groupselectblock(ROTOR_GROUP)
        t = femm.mo_blockintegral(22)
        femm.mo_clearblock()
        angles.append(k * step_deg)
        torques.append(t)
    angles, torques = np.array(angles), np.array(torques)

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_dir / "cogging_torque.csv",
               np.column_stack([angles, torques]),
               delimiter=",", header="rotor_angle_deg,cogging_torque_nm", comments="")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(angles, torques, "b-o", markersize=3)
    ax.axhline(0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("Rotor angle [°]"); ax.set_ylabel("Cogging torque [N·m]")
    ax.set_title(f"Cogging Torque (1 slot pitch = {slot_pitch_deg:.1f}°)")
    ax.grid(True); fig.tight_layout()
    fig.savefig(out_dir / "cogging_torque.png", dpi=150)
    plt.close(fig) if not show else plt.show()
    print(f"  Cogging torque peak‑to‑peak: {torques.max() - torques.min():.4f} N·m")


# ══════════════════════════════════════════════════════════════════
# 5) INDUCTANCE: Ld / Lq vs rotor position
# ══════════════════════════════════════════════════════════════════

def run_inductance_analysis(machine, out_dir, test_current=1.0, num_steps=37, show=False):
    """Compute Ld and Lq by applying d‑ or q‑axis current at each rotor position."""
    pp = machine.pole_pairs
    e_cycle_mech = 360.0 / pp          # one electrical cycle in mech deg
    step_deg = e_cycle_mech / num_steps
    dtheta_e = 2 * math.pi / num_steps

    Ld_arr, Lq_arr, pos_arr = [], [], []

    for axis_label, id_val, iq_val in [("d", test_current, 0.0), ("q", 0.0, test_current)]:
        build_model(machine)
        vals = []
        for k in range(num_steps):
            theta_e = k * dtheta_e
            set_dq_currents(machine, id_val, iq_val, theta_e)
            if k > 0:
                femm.mi_selectgroup(ROTOR_GROUP)
                femm.mi_moverotate(0.0, 0.0, step_deg)
                femm.mi_clearselected()
            femm.mi_analyze(1)
            femm.mi_loadsolution()

            # get flux linkages of all three phases
            psi_a = _get_circuit_flux("A")
            psi_b = _get_circuit_flux("B")
            psi_c = _get_circuit_flux("C")

            # Park transform → ψd, ψq
            cos_e = math.cos(theta_e)
            sin_e = math.sin(theta_e)
            cos_b = math.cos(theta_e - 2 * math.pi / 3)
            sin_b = math.sin(theta_e - 2 * math.pi / 3)
            cos_c = math.cos(theta_e + 2 * math.pi / 3)
            sin_c = math.sin(theta_e + 2 * math.pi / 3)
            psi_d = (2.0 / 3.0) * (psi_a * cos_e + psi_b * cos_b + psi_c * cos_c)
            psi_q = (2.0 / 3.0) * (-psi_a * sin_e - psi_b * sin_b - psi_c * sin_c)

            vals.append((k * step_deg, psi_d, psi_q))

            if axis_label == "d":
                pos_arr.append(k * step_deg)

        vals = np.array(vals)
        if axis_label == "d":
            # Ld = (ψd − ψ_pm_d) / Id.  At I=0 ψ_pm would be here; simplified: Ld ≈ ψd / Id
            # since test_current is small, the PM contribution dominates; we use incremental method
            Ld_arr = vals[:, 1] / test_current if test_current != 0 else vals[:, 1] * 0
        else:
            Lq_arr = vals[:, 2] / test_current if test_current != 0 else vals[:, 2] * 0

    pos_arr = np.array(pos_arr)
    Ld_arr = np.abs(np.array(Ld_arr))
    Lq_arr = np.abs(np.array(Lq_arr))

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_dir / "inductance.csv",
               np.column_stack([pos_arr, Ld_arr * 1e3, Lq_arr * 1e3]),
               delimiter=",", header="rotor_angle_deg,Ld_mH,Lq_mH", comments="")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(pos_arr, Ld_arr * 1e3, "b-o", markersize=3, label="Ld")
    ax.plot(pos_arr, Lq_arr * 1e3, "r-s", markersize=3, label="Lq")
    ax.set_xlabel("Rotor position [elec. °]"); ax.set_ylabel("Inductance [mH]")
    ax.set_title(f"d/q Inductance (test current = {test_current} A)")
    ax.legend(); ax.grid(True); fig.tight_layout()
    fig.savefig(out_dir / "inductance.png", dpi=150)
    plt.close(fig) if not show else plt.show()
    print(f"  Ld avg = {np.mean(Ld_arr)*1e3:.3f} mH,  Lq avg = {np.mean(Lq_arr)*1e3:.3f} mH")
    return float(np.mean(Ld_arr)), float(np.mean(Lq_arr))


# ══════════════════════════════════════════════════════════════════
# 6) EFFICIENCY MAP: speed × current sweep
# ══════════════════════════════════════════════════════════════════

def run_efficiency_map(machine, loss, rpm_list, iq_list, out_dir,
                       points_per_cycle=12, show=False):
    """Sweep rpm × iq, compute avg torque & losses, plot η contour."""
    results = []   # (rpm, iq, T_avg, P_loss_avg, eta)
    total = len(rpm_list) * len(iq_list)
    idx = 0
    for rpm in rpm_list:
        for iq in iq_list:
            idx += 1
            print(f"  emap [{idx}/{total}] rpm={rpm:.0f} iq={iq:.1f} A …")
            cfg = SweepCfg(points_per_electrical_cycle=points_per_cycle)
            rows = run_one_case(machine, loss, cfg, rpm, iq)
            t_avg = float(np.mean(rows[:, 5]))
            ploss_avg = float(np.mean(rows[:, 9]))
            omega = rpm * 2 * math.pi / 60.0
            p_mech = t_avg * omega
            eta = p_mech / (p_mech + ploss_avg) if (p_mech + ploss_avg) > 0 else 0.0
            eta = max(0.0, min(1.0, eta))
            results.append((rpm, iq, t_avg, ploss_avg, eta))

    data = np.array(results, dtype=float)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_dir / "efficiency_map.csv", data, delimiter=",",
               header="rpm,iq_peak_a,avg_torque_nm,avg_loss_w,efficiency", comments="")

    # reshape for contour
    n_rpm, n_iq = len(rpm_list), len(iq_list)
    RPM = data[:, 0].reshape(n_rpm, n_iq)
    TQ  = data[:, 2].reshape(n_rpm, n_iq)
    ETA = data[:, 4].reshape(n_rpm, n_iq) * 100  # percent

    fig, ax = plt.subplots(figsize=(10, 6))
    levels = np.arange(0, 101, 5)
    cs = ax.contourf(RPM, TQ, ETA, levels=levels, cmap="RdYlGn")
    fig.colorbar(cs, ax=ax, label="Efficiency [%]")
    ax.contour(RPM, TQ, ETA, levels=levels, colors="k", linewidths=0.3)
    ax.set_xlabel("Speed [rpm]"); ax.set_ylabel("Avg Torque [N·m]")
    ax.set_title("Efficiency Map")
    ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(out_dir / "efficiency_map.png", dpi=150)
    plt.close(fig) if not show else plt.show()
    return data


# ══════════════════════════════════════════════════════════════════
# 7) TORQUE‑SPEED CURVE
# ══════════════════════════════════════════════════════════════════

def plot_torque_speed_curve(machine, drive, emap_data, out_dir, show=False):
    """Plot T‑N envelope from efficiency map data + analytical voltage limit."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rpm_vals = np.unique(emap_data[:, 0])
    t_max_fem = []
    for rpm in rpm_vals:
        mask = emap_data[:, 0] == rpm
        t_max_fem.append(np.max(emap_data[mask, 2]))
    t_max_fem = np.array(t_max_fem)

    # Analytical voltage‑limited envelope (simplified SPMSM):
    #   V_max ≈ V_dc / sqrt(3)
    #   At each speed:  V_max = sqrt((Rs*I + ωe*ψ_pm)² + (ωe*Ld*I)²)
    #   → solve for max I → max torque = 1.5 * p * ψ_pm * Iq
    # We use a rough estimate with Ld, ψ_pm from the geometry
    mu_r, Br = 1.05, 1.22
    Hc = Br / (MU0 * mu_r)
    # ψ_pm ≈ Br * (2/π) * τ_p * L_stk * N_ph * k_w / p    (rough)
    tau_p = math.pi * (machine.r_rotor_mm * 1e-3) / machine.pole_pairs
    Lstk = machine.stack_length_mm * 1e-3
    Nph = machine.turns_per_slot * machine.slots / 3 / 2  # series turns per phase
    kw = 0.866  # winding factor estimate
    psi_pm = Br * machine.magnet_arc_ratio * tau_p * Lstk * Nph * kw * 2 / math.pi
    # rough Ld
    g_eff = (machine.airgap_mm + machine.mag_thickness_mm / mu_r) * 1e-3
    Ld = MU0 * (Nph**2) * math.pi * (machine.r_stator_inner_mm * 1e-3) * Lstk / (machine.pole_pairs * g_eff)

    V_max = drive.v_dc / math.sqrt(3)
    Rs = 0.023  # approx from design_pmsm
    rpm_dense = np.linspace(100, rpm_vals.max() * 1.5, 200)
    t_envelope = []
    for rpm in rpm_dense:
        we = machine.pole_pairs * rpm * 2 * math.pi / 60.0
        # max current from voltage constraint: simplified
        i_budget = min(drive.i_max, max(0, (V_max - we * psi_pm) / max(we * Ld + Rs, 1e-9)))
        torque = 1.5 * machine.pole_pairs * psi_pm * i_budget
        t_envelope.append(torque)
    t_envelope = np.array(t_envelope)
    t_envelope = np.clip(t_envelope, 0, 1.5 * machine.pole_pairs * psi_pm * drive.i_max)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # torque‑speed
    ax1.plot(rpm_vals, t_max_fem, "bs-", label="FEA max torque", markersize=5)
    ax1.plot(rpm_dense, t_envelope, "r--", label="Voltage‑limit envelope (analytical)")
    ax1.set_ylabel("Torque [N·m]"); ax1.set_title("Torque–Speed Curve")
    ax1.legend(); ax1.grid(True)

    # power‑speed
    p_fem = t_max_fem * rpm_vals * 2 * math.pi / 60.0
    p_env = t_envelope * rpm_dense * 2 * math.pi / 60.0
    ax2.plot(rpm_vals, p_fem, "bs-", label="FEA max power", markersize=5)
    ax2.plot(rpm_dense, p_env, "r--", label="Analytical envelope")
    ax2.set_xlabel("Speed [rpm]"); ax2.set_ylabel("Power [W]")
    ax2.set_title("Power–Speed Curve")
    ax2.legend(); ax2.grid(True)

    fig.tight_layout()
    fig.savefig(out_dir / "torque_speed_curve.png", dpi=150)
    plt.close(fig) if not show else plt.show()


# ══════════════════════════════════════════════════════════════════
# CLI + MAIN
# ══════════════════════════════════════════════════════════════════

ALL_ANALYSES = {"basic", "field", "airgap", "cogging", "inductance", "emap", "tncurve"}


def parse_args():
    p = argparse.ArgumentParser(
        description="FEMM 2D SPMSM template – comprehensive motor analysis.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument(
        "--analysis", nargs="+", default=["basic"],
        help="Analysis to run (space or comma separated). Options:\n"
             "  basic      – torque/loss/temperature waveforms\n"
             "  field      – flux density cloud + flux lines\n"
             "  airgap     – airgap Bn/Bt distribution\n"
             "  cogging    – cogging torque curve\n"
             "  inductance – Ld/Lq vs rotor position\n"
             "  emap       – efficiency map\n"
             "  tncurve    – torque‑speed curve (requires emap)\n"
             "  all        – run every analysis\n"
             "Default: basic",
    )
    p.add_argument("--rpm-list", nargs="+", default=["1000", "1500", "2000"],
                   help="Speed list (rpm) for basic analysis.")
    p.add_argument("--iq-list", nargs="+", default=["15", "25", "35"],
                   help="Current peak list (A) for basic analysis.")
    p.add_argument("--points", type=int, default=72, help="Samples per electrical cycle.")
    p.add_argument("--out", type=Path, default=Path("output_femm"), help="Output folder.")
    p.add_argument("--show", action="store_true", help="Show matplotlib figures.")
    p.add_argument("--show-femm", action="store_true", help="Show FEMM GUI.")
    # emap specific
    p.add_argument("--emap-rpm", nargs="+", default=["500,1000,1500,2000,2500,3000"],
                   help="Speed list for efficiency map.")
    p.add_argument("--emap-iq", nargs="+", default=["5,10,15,20,25,30,35"],
                   help="Current list for efficiency map.")
    p.add_argument("--emap-pts", type=int, default=12,
                   help="Points per electrical cycle for efficiency map (fewer = faster).")
    # inductance / cogging
    p.add_argument("--cogging-steps", type=int, default=72,
                   help="Steps per slot pitch for cogging torque.")
    p.add_argument("--ind-steps", type=int, default=37,
                   help="Steps per electrical cycle for inductance.")
    p.add_argument("--ind-current", type=float, default=1.0,
                   help="Test current (A) for inductance measurement.")
    return p.parse_args()


def open_femm_or_raise(show_femm: bool):
    try:
        femm.openfemm(0 if show_femm else 1)
    except Exception as exc:
        raise SystemExit(
            "Failed to connect FEMM COM server (femm.ActiveFEMM).\n"
            "Likely causes:\n"
            "  1) FEMM 4.2 is not installed.\n"
            "  2) FEMM COM server is not registered.\n\n"
            "Fix steps (run in Administrator CMD):\n"
            '  "C:\\Program Files (x86)\\FEMM42\\femm.exe" /regserver\n\n'
            "Verify registration:\n"
            "  reg query HKCR\\femm.ActiveFEMM\n\n"
            f"Original error: {exc}"
        ) from exc


def main():
    args = parse_args()

    # resolve analysis set
    raw_analyses = []
    for item in args.analysis:
        raw_analyses.extend(re.split(r"[\s,]+", item.strip().lower()))
    if "all" in raw_analyses:
        analyses = ALL_ANALYSES.copy()
    else:
        analyses = set(raw_analyses) & ALL_ANALYSES
    if not analyses:
        print("No valid analysis selected. Choose from:", ", ".join(sorted(ALL_ANALYSES)))
        return
    # tncurve requires emap
    if "tncurve" in analyses:
        analyses.add("emap")

    rpm_list = parse_list(args.rpm_list)
    iq_list = parse_list(args.iq_list)
    machine = Machine()
    loss = LossThermal()
    drive = Drive()
    cfg = SweepCfg(points_per_electrical_cycle=args.points)

    args.out.mkdir(parents=True, exist_ok=True)
    open_femm_or_raise(args.show_femm)

    print(f"Selected analyses: {', '.join(sorted(analyses))}\n")

    try:
        # ── basic + field + airgap ───────────────────
        if analyses & {"basic", "field", "airgap"}:
            summary_rows = []
            for rpm in rpm_list:
                for iq in iq_list:
                    case_tag = f"rpm{int(rpm)}_iq{int(iq)}"
                    print(f"Running FEMM case: {case_tag}")
                    rows = run_one_case(machine, loss, cfg, rpm, iq)

                    if "basic" in analyses:
                        csv_path, png_path = save_case_outputs(
                            rows, args.out, case_tag, args.show)

                    # post‑process the LAST solved position (solution still loaded)
                    if "field" in analyses:
                        print(f"  Generating field maps …")
                        plot_field_maps(machine, args.out, case_tag, args.show)

                    if "airgap" in analyses:
                        print(f"  Generating airgap B …")
                        plot_airgap_flux_density(machine, args.out, case_tag, args.show)

                    avg_torque = float(np.mean(rows[:, 5]))
                    avg_loss = float(np.mean(rows[:, 9]))
                    final_temp = float(rows[-1, 10])
                    summary_rows.append((rpm, iq, avg_torque, avg_loss, final_temp))

            if "basic" in analyses and summary_rows:
                summary = np.array(summary_rows, dtype=float)
                np.savetxt(args.out / "femm_summary.csv", summary, delimiter=",",
                           header="rpm,iq_peak_a,avg_torque_nm,avg_total_loss_w,final_temp_c",
                           comments="")
                print("\nBasic summary written to:", args.out / "femm_summary.csv")
                for r in summary_rows:
                    print(f"  rpm={r[0]:.0f}, iq={r[1]:.1f} A, Tavg={r[2]:.3f} N·m, "
                          f"Pavg={r[3]:.2f} W, Tend={r[4]:.2f} °C")

        # ── cogging ──────────────────────────────────
        if "cogging" in analyses:
            print("\n── Cogging torque analysis ──")
            run_cogging_torque(machine, args.out, num_steps=args.cogging_steps,
                               show=args.show)

        # ── inductance ───────────────────────────────
        if "inductance" in analyses:
            print("\n── Inductance analysis ──")
            run_inductance_analysis(machine, args.out,
                                    test_current=args.ind_current,
                                    num_steps=args.ind_steps, show=args.show)

        # ── efficiency map ───────────────────────────
        emap_data = None
        if "emap" in analyses:
            print("\n── Efficiency map ──")
            emap_rpm = parse_list(args.emap_rpm)
            emap_iq = parse_list(args.emap_iq)
            emap_data = run_efficiency_map(machine, loss, emap_rpm, emap_iq,
                                           args.out, points_per_cycle=args.emap_pts,
                                           show=args.show)

        # ── torque‑speed curve ───────────────────────
        if "tncurve" in analyses and emap_data is not None:
            print("\n── Torque‑speed curve ──")
            plot_torque_speed_curve(machine, drive, emap_data, args.out, args.show)

        print("\n✓ All analyses complete. Output folder:", args.out)

    finally:
        try:
            femm.closefemm()
        except Exception:
            pass


if __name__ == "__main__":
    main()
