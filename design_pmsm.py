from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

MU0 = 4e-7 * np.pi  # 真空磁导率


@dataclass
class Geometry:
    # 几何相关参数（SPMSM 简化模型）
    pole_pairs: int = 4
    d_in: float = 0.11
    d_out: float = 0.18
    stack_len: float = 0.08
    air_gap: float = 0.0008
    mag_thick: float = 0.003
    mag_arc_ratio: float = 0.75
    turns_phase: int = 70
    mean_turn_len: float = 0.42
    a_cond_phase: float = 22e-6
    core_mass: float = 9.0
    magnet_mass: float = 1.2


@dataclass
class Material:
    # 材料与经验损耗参数
    br: float = 1.22
    mu_r_mag: float = 1.05
    rho_cu_20: float = 1.724e-8
    alpha_cu: float = 0.00393
    kh: float = 0.90
    ke: float = 0.0030
    b_core_base: float = 1.10
    k_pm: float = 0.15


@dataclass
class Thermal:
    # 一阶热网络参数
    t_amb: float = 25.0
    r_th: float = 0.45
    c_th: float = 2600.0


@dataclass
class Drive:
    # 逆变器/控制器约束
    v_dc: float = 310.0
    m_max: float = 0.98
    i_max: float = 55.0


@dataclass
class SimConfig:
    dt: float = 1e-4
    t_end: float = 4.0


def estimate_em_params(geo: Geometry, mat: Material):
    """由尺寸与材料粗估 R、Ld/Lq、永磁磁链（工程初值）。"""
    tau_p = np.pi * geo.d_in / (2.0 * geo.pole_pairs)
    pole_arc = geo.mag_arc_ratio * tau_p
    b_gap = mat.br * geo.mag_thick / (geo.mag_thick + mat.mu_r_mag * geo.air_gap)
    phi_p = b_gap * pole_arc * geo.stack_len

    winding_factor = 0.95
    psi_pm = winding_factor * geo.turns_phase * phi_p

    g_eff = geo.air_gap + geo.mag_thick / mat.mu_r_mag
    ld = MU0 * geo.turns_phase**2 * geo.stack_len * tau_p / (2.0 * g_eff)
    lq = ld  # SPMSM 近似 Ld=Lq

    r20 = mat.rho_cu_20 * geo.turns_phase * geo.mean_turn_len / geo.a_cond_phase
    return r20, ld, lq, psi_pm


def speed_profile_rpm(time_s: float):
    # 示例速度工况：1500rpm 基值 + 低频扰动
    return 1500.0 + 300.0 * np.sin(2.0 * np.pi * 0.5 * time_s)


def torque_cmd_profile(time_s: float):
    # 示例转矩工况：1s 后阶跃，同时叠加小幅波动
    base = 18.0 if time_s < 1.0 else 28.0
    return base + 4.0 * np.sin(2.0 * np.pi * 2.0 * time_s)


def run_simulation(
    geo: Geometry,
    mat: Material,
    thermal: Thermal,
    drive: Drive,
    cfg: SimConfig,
):
    r20, ld, lq, psi_pm = estimate_em_params(geo, mat)

    n = int(np.floor(cfg.t_end / cfg.dt)) + 1
    ts = np.arange(n) * cfg.dt

    # 结果缓存
    te_cmd_arr = np.zeros(n)
    te_arr = np.zeros(n)
    pcu_arr = np.zeros(n)
    pfe_arr = np.zeros(n)
    ppm_arr = np.zeros(n)
    ploss_arr = np.zeros(n)
    temp_arr = np.zeros(n)
    rs_arr = np.zeros(n)

    # 状态量：dq 电流与绕组温度
    i_d = 0.0
    i_q = 0.0
    temp = thermal.t_amb

    # 简化电流环：按带宽设置 P 控增益
    bw_i = 1200.0
    kp_d = bw_i * ld
    kp_q = bw_i * lq

    for k, t in enumerate(ts):
        rpm = speed_profile_rpm(t)
        w_m = rpm * 2.0 * np.pi / 60.0
        w_e = geo.pole_pairs * w_m

        # MTPA 简化：SPMSM 默认 id=0，通过目标转矩反算 iq*
        te_cmd = torque_cmd_profile(t)
        i_d_ref = 0.0
        i_q_ref = te_cmd / (1.5 * geo.pole_pairs * max(psi_pm, 1e-6))

        # 电流限幅
        i_ref_mag = np.hypot(i_d_ref, i_q_ref)
        if i_ref_mag > drive.i_max:
            scale = drive.i_max / i_ref_mag
            i_d_ref *= scale
            i_q_ref *= scale

        # 铜阻随温度变化
        r_s = r20 * (1.0 + mat.alpha_cu * (temp - 20.0))

        # dq 电压前馈 + P 控制
        v_d_ff = r_s * i_d - w_e * lq * i_q
        v_q_ff = r_s * i_q + w_e * (ld * i_d + psi_pm)

        v_d = v_d_ff + kp_d * (i_d_ref - i_d)
        v_q = v_q_ff + kp_q * (i_q_ref - i_q)

        # 受母线电压限制（SVPWM 线性区近似）
        v_max = drive.m_max * drive.v_dc / np.sqrt(3.0)
        v_mag = np.hypot(v_d, v_q)
        if v_mag > v_max:
            scale = v_max / v_mag
            v_d *= scale
            v_q *= scale

        # 电流状态方程离散积分
        di_d = (v_d - r_s * i_d + w_e * lq * i_q) / ld
        di_q = (v_q - r_s * i_q - w_e * (ld * i_d + psi_pm)) / lq
        i_d += di_d * cfg.dt
        i_q += di_q * cfg.dt

        # 电磁转矩
        te = 1.5 * geo.pole_pairs * (psi_pm * i_q + (ld - lq) * i_d * i_q)

        # 铜耗（3 相）
        i_phase_rms = np.hypot(i_d, i_q) / np.sqrt(2.0)
        pcu = 3.0 * (i_phase_rms**2) * r_s

        # 铁耗/磁钢损耗（经验公式）
        f_e = abs(w_e) / (2.0 * np.pi)
        b_pk = mat.b_core_base + 0.002 * abs(i_q)
        pfe = geo.core_mass * (mat.kh * f_e * b_pk**2 + mat.ke * (f_e * b_pk) ** 2)
        ppm = geo.magnet_mass * mat.k_pm * (f_e * b_pk) ** 2

        # 一阶热网络：Cth*dT/dt = Ploss - (T-Tamb)/Rth
        ploss = pcu + pfe + ppm
        dtemp = (ploss - (temp - thermal.t_amb) / thermal.r_th) / thermal.c_th
        temp += dtemp * cfg.dt

        te_cmd_arr[k] = te_cmd
        te_arr[k] = te
        pcu_arr[k] = pcu
        pfe_arr[k] = pfe
        ppm_arr[k] = ppm
        ploss_arr[k] = ploss
        temp_arr[k] = temp
        rs_arr[k] = r_s

    return {
        "t": ts,
        "te_cmd": te_cmd_arr,
        "te": te_arr,
        "pcu": pcu_arr,
        "pfe": pfe_arr,
        "ppm": ppm_arr,
        "ploss": ploss_arr,
        "temp": temp_arr,
        "rs": rs_arr,
        "r20": r20,
        "ld": ld,
        "lq": lq,
        "psi_pm": psi_pm,
    }


def save_csv(results, out_dir: Path):
    # 导出时域波形到 CSV，便于后处理
    out_dir.mkdir(parents=True, exist_ok=True)
    table = np.column_stack(
        [
            results["t"],
            results["te_cmd"],
            results["te"],
            results["pcu"],
            results["pfe"],
            results["ppm"],
            results["ploss"],
            results["temp"],
            results["rs"],
        ]
    )
    header = (
        "t_s,torque_cmd_nm,torque_em_nm,copper_loss_w,core_loss_w,"
        "magnet_loss_w,total_loss_w,temp_c,phase_res_ohm"
    )
    np.savetxt(out_dir / "pmsm_waveforms.csv", table, delimiter=",", header=header, comments="")


def plot_results(results, save_png: bool, out_dir: Path, show: bool):
    # 3 张子图：转矩、损耗、温度/相电阻
    t = results["t"]
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axs[0].plot(t, results["te_cmd"], "--", label="Torque cmd")
    axs[0].plot(t, results["te"], label="Electromagnetic torque")
    axs[0].set_ylabel("Torque [N.m]")
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(t, results["pcu"], label="Copper loss")
    axs[1].plot(t, results["pfe"], label="Core loss")
    axs[1].plot(t, results["ppm"], label="Magnet loss")
    axs[1].plot(t, results["ploss"], "k", linewidth=1.2, label="Total loss")
    axs[1].set_ylabel("Power loss [W]")
    axs[1].grid(True)
    axs[1].legend()

    ax2 = axs[2].twinx()
    axs[2].plot(t, results["temp"], "r", label="Winding temperature")
    ax2.plot(t, results["rs"], "b--", label="Phase resistance")
    axs[2].set_ylabel("Temperature [degC]", color="r")
    ax2.set_ylabel("Rs [ohm]", color="b")
    axs[2].set_xlabel("Time [s]")
    axs[2].grid(True)

    h1, l1 = axs[2].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axs[2].legend(h1 + h2, l1 + l2, loc="upper left")

    fig.tight_layout()

    if save_png:
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "pmsm_waveforms.png", dpi=150)

    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_args():
    # 命令行参数：方便扫母线电压、仿真时长和输出目录
    parser = argparse.ArgumentParser(description="SPMSM dq model with torque/loss/temperature waveforms.")
    parser.add_argument("--vdc", type=float, default=310.0, help="DC bus voltage [V].")
    parser.add_argument("--dt", type=float, default=1e-4, help="Simulation time step [s].")
    parser.add_argument("--t-end", type=float, default=4.0, help="Simulation end time [s].")
    parser.add_argument("--out", type=Path, default=Path("output"), help="Output directory.")
    parser.add_argument("--save-csv", action="store_true", help="Save waveform CSV.")
    parser.add_argument("--save-png", action="store_true", help="Save waveform PNG.")
    parser.add_argument("--no-show", action="store_true", help="Disable interactive plot window.")
    return parser.parse_args()


def main():
    args = parse_args()

    # 这里可替换成你的真实电机参数
    geo = Geometry()
    mat = Material()
    thermal = Thermal()
    drive = Drive(v_dc=args.vdc)
    cfg = SimConfig(dt=args.dt, t_end=args.t_end)

    results = run_simulation(geo, mat, thermal, drive, cfg)

    print(
        "Estimated params: "
        f"R20={results['r20']:.4f} ohm, "
        f"Ld={results['ld']*1e3:.3f} mH, "
        f"Lq={results['lq']*1e3:.3f} mH, "
        f"psi_pm={results['psi_pm']:.4f} Wb"
    )

    if args.save_csv:
        save_csv(results, args.out)
        print(f"Saved CSV: {args.out / 'pmsm_waveforms.csv'}")

    plot_results(results, save_png=args.save_png, out_dir=args.out, show=not args.no_show)
    if args.save_png:
        print(f"Saved figure: {args.out / 'pmsm_waveforms.png'}")


if __name__ == "__main__":
    main()
