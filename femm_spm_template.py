from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import femm
except ImportError as exc:
    raise SystemExit(
        "pyfemm is not installed. Run: pip install pyfemm, and install FEMM 4.2 first."
    ) from exc


ROTOR_GROUP = 2
STATOR_GROUP = 3
COIL_GROUP = 4


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
class SweepCfg:
    points_per_electrical_cycle: int = 72


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


def set_block(material: str, x: float, y: float, circuit: str = "<None>", mag_dir: float = 0.0, group: int = 0, turns: int = 0):
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    femm.mi_setblockprop(material, 1, 0.0, circuit, mag_dir, group, turns)
    femm.mi_clearselected()


def ensure_materials():
    for name in ["Air", "M-19 Steel", "NdFeB 40 MGOe", "Copper"]:
        try:
            femm.mi_getmaterial(name)
        except Exception:
            pass


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
    set_block(
        "M-19 Steel",
        0.5 * (machine.r_shaft_mm + machine.r_rotor_mm),
        0.0,
        group=ROTOR_GROUP,
    )
    set_block(
        "Air",
        0.5 * (machine.r_mag_outer_mm + machine.r_stator_inner_mm),
        0.0,
    )
    set_block(
        "M-19 Steel",
        0.5 * (machine.r_slot_outer_mm + machine.r_stator_outer_mm),
        12.0,
        group=STATOR_GROUP,
    )
    set_block("Air", 0.5 * (machine.r_air_outer_mm + machine.r_stator_outer_mm), 0.0)

    pole_count = 2 * machine.pole_pairs
    pole_pitch = 360.0 / pole_count
    magnet_span = pole_pitch * machine.magnet_arc_ratio

    for k in range(pole_count):
        center = k * pole_pitch
        a1 = center - 0.5 * magnet_span
        a2 = center + 0.5 * magnet_span
        x, y = add_wedge_annulus(
            machine.r_rotor_mm,
            machine.r_mag_outer_mm,
            a1,
            a2,
            group=ROTOR_GROUP,
        )
        mag_dir = center if (k % 2 == 0) else (center + 180.0)
        set_block("NdFeB 40 MGOe", x, y, mag_dir=mag_dir, group=ROTOR_GROUP)

    phase_map = [
        ("A", +1),
        ("C", -1),
        ("B", +1),
        ("A", -1),
        ("C", +1),
        ("B", -1),
        ("A", +1),
        ("C", -1),
        ("B", +1),
        ("A", -1),
        ("C", +1),
        ("B", -1),
    ]
    if machine.slots != len(phase_map):
        raise ValueError("This template currently expects 12 slots.")

    slot_pitch = 360.0 / machine.slots
    slot_span = slot_pitch * 0.55
    r_slot_in = machine.r_stator_inner_mm + 0.2
    r_slot_out = machine.r_slot_outer_mm

    for idx, (phase, sign) in enumerate(phase_map):
        center = idx * slot_pitch
        a1 = center - 0.5 * slot_span
        a2 = center + 0.5 * slot_span
        x, y = add_wedge_annulus(r_slot_in, r_slot_out, a1, a2, group=COIL_GROUP)
        set_block(
            "Copper",
            x,
            y,
            circuit=phase,
            group=COIL_GROUP,
            turns=sign * machine.turns_per_slot,
        )


def set_phase_currents(i_peak: float, theta_e_rad: float):
    ia = i_peak * math.sin(theta_e_rad)
    ib = i_peak * math.sin(theta_e_rad - 2.0 * math.pi / 3.0)
    ic = i_peak * math.sin(theta_e_rad + 2.0 * math.pi / 3.0)

    femm.mi_modifycircprop("A", 1, ia)
    femm.mi_modifycircprop("B", 1, ib)
    femm.mi_modifycircprop("C", 1, ic)
    return ia, ib, ic


def estimate_core_loss(machine: Machine, loss: LossThermal, rpm: float):
    probe_r = 0.5 * (machine.r_slot_outer_mm + machine.r_stator_outer_mm)
    b_max = 0.0
    for deg in [15, 45, 75, 105, 135, 165]:
        x, y = polar_xy(probe_r, deg)
        bx, by = femm.mo_getb(x, y)
        b = math.sqrt(bx * bx + by * by)
        b_max = max(b_max, b)

    f_e = machine.pole_pairs * rpm / 60.0
    pfe = loss.core_mass_kg * (loss.kh * f_e * b_max * b_max + loss.ke * (f_e * b_max) ** 2)
    return pfe


def run_one_case(machine: Machine, loss: LossThermal, cfg: SweepCfg, rpm: float, iq_peak: float):
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
        pcu = rs * (ia * ia + ib * ib + ic * ic)
        pfe = estimate_core_loss(machine, loss, rpm)
        ppm = loss.magnet_loss_ratio * pfe
        ploss = pcu + pfe + ppm

        dtemp = (ploss - (temp - loss.t_amb) / loss.r_th) / loss.c_th
        temp += dtemp * dt

        rows.append((k * dt, theta_e, ia, ib, ic, torque, pcu, pfe, ppm, ploss, temp))

    return np.array(rows, dtype=float)


def save_case_outputs(rows: np.ndarray, out_dir: Path, case_tag: str, show: bool):
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / f"femm_waveforms_{case_tag}.csv"
    np.savetxt(
        csv_path,
        rows,
        delimiter=",",
        header=(
            "t_s,theta_e_rad,ia_a,ib_a,ic_a,torque_nm,"
            "pcu_w,pfe_w,ppm_w,ploss_w,temp_c"
        ),
        comments="",
    )

    t = rows[:, 0]
    torque = rows[:, 5]
    pcu = rows[:, 6]
    pfe = rows[:, 7]
    ppm = rows[:, 8]
    ploss = rows[:, 9]
    temp = rows[:, 10]

    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    axs[0].plot(t, torque, label="FEMM torque")
    axs[0].set_ylabel("Torque [N.m]")
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(t, pcu, label="Copper loss")
    axs[1].plot(t, pfe, label="Core loss")
    axs[1].plot(t, ppm, label="Magnet loss")
    axs[1].plot(t, ploss, "k", label="Total loss")
    axs[1].set_ylabel("Loss [W]")
    axs[1].grid(True)
    axs[1].legend()

    axs[2].plot(t, temp, "r", label="Winding temperature")
    axs[2].set_xlabel("Time [s]")
    axs[2].set_ylabel("Temperature [degC]")
    axs[2].grid(True)
    axs[2].legend()

    fig.tight_layout()
    png_path = out_dir / f"femm_waveforms_{case_tag}.png"
    fig.savefig(png_path, dpi=150)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return csv_path, png_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="FEMM 2D SPMSM template: auto model + speed/load sweep + torque/loss waveforms."
    )
    parser.add_argument(
        "--rpm-list",
        nargs="+",
        default=["1000", "1500", "2000"],
        help="Speed list in rpm, supports comma or space separated values.",
    )
    parser.add_argument(
        "--iq-list",
        nargs="+",
        default=["15", "25", "35"],
        help="Current-peak list in A, supports comma or space separated values.",
    )
    parser.add_argument("--points", type=int, default=72, help="Samples per electrical cycle.")
    parser.add_argument("--out", type=Path, default=Path("output_femm"), help="Output folder.")
    parser.add_argument("--show", action="store_true", help="Show matplotlib figures.")
    parser.add_argument("--show-femm", action="store_true", help="Show FEMM GUI during solve.")
    return parser.parse_args()



def open_femm_or_raise(show_femm: bool):
    try:
        femm.openfemm(0 if show_femm else 1)
    except Exception as exc:
        msg = str(exc)
        raise SystemExit(
            "Failed to connect FEMM COM server (femm.ActiveFEMM).\n"
            "Likely causes:\n"
            "  1) FEMM 4.2 is not installed.\n"
            "  2) FEMM COM server is not registered.\n\n"
            "Fix steps (run in Administrator CMD):\n"
            "  \"C:\\Program Files (x86)\\FEMM42\\femm.exe\" /regserver\n\n"
            "Verify registration:\n"
            "  reg query HKCR\\femm.ActiveFEMM\n\n"
            f"Original error: {msg}"
        ) from exc
def main():
    args = parse_args()
    rpm_list = parse_list(args.rpm_list)
    iq_list = parse_list(args.iq_list)

    machine = Machine()
    loss = LossThermal()
    cfg = SweepCfg(points_per_electrical_cycle=args.points)

    args.out.mkdir(parents=True, exist_ok=True)
    open_femm_or_raise(args.show_femm)

    summary_rows = []
    try:
        for rpm in rpm_list:
            for iq in iq_list:
                case_tag = f"rpm{int(rpm)}_iq{int(iq)}"
                print(f"Running FEMM case: {case_tag}")
                rows = run_one_case(machine, loss, cfg, rpm, iq)
                csv_path, png_path = save_case_outputs(rows, args.out, case_tag, args.show)
                avg_torque = float(np.mean(rows[:, 5]))
                avg_loss = float(np.mean(rows[:, 9]))
                final_temp = float(rows[-1, 10])
                summary_rows.append((rpm, iq, avg_torque, avg_loss, final_temp, csv_path.name, png_path.name))

        summary = np.array([[r[0], r[1], r[2], r[3], r[4]] for r in summary_rows], dtype=float)
        np.savetxt(
            args.out / "femm_summary.csv",
            summary,
            delimiter=",",
            header="rpm,iq_peak_a,avg_torque_nm,avg_total_loss_w,final_temp_c",
            comments="",
        )

        print("Done. Summary written to:", args.out / "femm_summary.csv")
        for row in summary_rows:
            print(
                f"  rpm={row[0]:.0f}, iq={row[1]:.1f} A, "
                f"Tavg={row[2]:.3f} N.m, Pavg={row[3]:.2f} W, "
                f"Tend={row[4]:.2f} C, csv={row[5]}, png={row[6]}"
            )
    finally:
        try:
            femm.closefemm()
        except Exception:
            pass


if __name__ == "__main__":
    main()


