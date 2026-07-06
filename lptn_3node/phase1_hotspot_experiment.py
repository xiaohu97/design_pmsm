"""
Phase-1 热点子节点实验: 绕组拆成 热点(C_h) + 本体(C_b)
======================================================
物理动机: 堵转时电流集中在某一相, 铜损集中于局部导体; 传感器测到的是热点温度。
模型:
  dTh/dt = ( γ·P(Th)         - (Th-Tb)/R_h )                    / C_h
  dTb/dt = ( (1-γ)·P(Tb) + (Th-Tb)/R_h - (Tb-T3)/R13 )          / C_b
  P(T)   = a_eff·(1+α·(T-T_ref))·I²      (各自用本地温度做电阻修正)
测量: T1 传感器 = Th;  T3 实测作边界;  每段 Th0=Tb0=实测T1[0]

对比: 单节点联合拟合 (基线) vs 热点模型 (Seg1单独 / 全段联合)
用法: python phase1_hotspot_experiment.py <CSV1> [CSV2] ...
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import (CONFIG, load_multiple, extract_lock_segments,
                              p1_sim_segment, fit_phase1_segments)

# 参数: C_h, C_b, R_h, R13, a_eff, gamma
HS_NAMES = ['C_h', 'C_b', 'R_h', 'R13', 'a_eff', 'gamma']
HS_BOUNDS = [(1, 200), (20, 800), (0.05, 5), (0.05, 10), (0.1, 20), (0.05, 1.0)]
OUTPUT_PLOT = 'phase1_hotspot_compare.png'
MAX_SUBSTEP = 25


def hs_sim_segment(seg, p, alpha_cu, T_ref):
    """两节点绕组前向 Euler (快模态自动细分子步), 返回 (Th, Tb)"""
    C_h, C_b, R_h, R13, a_eff, gamma = p
    I2 = seg['I2']; T3 = seg['T3']; dt_arr = seg['dt']
    N = len(seg['t'])
    Th = np.empty(N); Tb = np.empty(N)
    Th[0] = Tb[0] = seg['T1'][0]
    inv_Ch = 1.0/C_h; inv_Cb = 1.0/C_b
    inv_Rh = 1.0/R_h; inv_R13 = 1.0/R13
    tau_fast = C_h * R_h          # 最快模态量级
    for k in range(N - 1):
        dt = dt_arr[k]
        m = int(dt / (0.4 * tau_fast)) + 1
        if m > MAX_SUBSTEP:
            return None, None     # 太刚性, 判为无效参数
        h = dt / m
        th, tb = Th[k], Tb[k]
        for _ in range(m):
            Ph = gamma * a_eff * (1.0 + alpha_cu*(th - T_ref)) * I2[k]
            Pb = (1.0-gamma) * a_eff * (1.0 + alpha_cu*(tb - T_ref)) * I2[k]
            q_int = (th - tb) * inv_Rh
            dth = (Ph - q_int) * inv_Ch
            dtb = (Pb + q_int - (tb - T3[k]) * inv_R13) * inv_Cb
            th += dth * h; tb += dtb * h
        Th[k+1] = th; Tb[k+1] = tb
    return Th, Tb


def hs_cost(p, segments, alpha_cu, T_ref):
    sq = dsq = 0.0
    n = dn = 0
    try:
        for seg in segments:
            Th, _ = hs_sim_segment(seg, p, alpha_cu, T_ref)
            if Th is None or np.any(np.isnan(Th)) or np.any(np.abs(Th) > 500):
                return 1e6
            sq += np.sum((Th - seg['T1'])**2); n += len(Th)
            dTh = np.diff(Th) / seg['dt_safe']
            dsq += np.sum((dTh - seg['dT1_meas'])**2); dn += len(dTh)
        return sq/n + 0.5*dsq/dn
    except Exception:
        return 1e6


def fit_hotspot(segments, config, seeds=(42, 123)):
    alpha_cu = config.get('alpha_cu', 0.00393)
    T_ref = config.get('T_ref', 25.0)
    best_c, best_x = 1e10, None
    for seed in seeds:
        print(f'    DE-HS (seed={seed})...', end=' ', flush=True)
        t0 = time.time()
        r = differential_evolution(hs_cost, HS_BOUNDS,
            args=(segments, alpha_cu, T_ref), seed=seed,
            maxiter=250, popsize=15, tol=1e-9, polish=False, workers=1,
            mutation=(0.5, 1.5), recombination=0.9)
        r2 = minimize(hs_cost, r.x, args=(segments, alpha_cu, T_ref),
            method='L-BFGS-B', bounds=HS_BOUNDS,
            options={'maxiter': 400, 'ftol': 1e-14})
        print(f'cost={r2.fun:.4f} ({time.time()-t0:.0f}s)')
        if r2.fun < best_c:
            best_c, best_x = r2.fun, r2.x
    return best_x, best_c


def seg_rmse(T_sim, T_meas):
    return np.sqrt(np.mean((T_sim - T_meas)**2))


def print_params(p, tag):
    C_h, C_b, R_h, R13, a_eff, gamma = p
    print(f'  [{tag}] ' + ', '.join(f'{n}={v:.4g}' for n, v in zip(HS_NAMES, p)))
    print(f'         tau_hot=C_h·R_h={C_h*R_h:.1f}s, tau_bulk=C_b·R13={C_b*R13:.1f}s, '
          f'a_eff·R13={a_eff*R13:.3f} K/(W·A²)... γ={gamma:.2f}')
    for v, (lo, hi), n in zip(p, HS_BOUNDS, HS_NAMES):
        m = 1e-3 * (hi - lo)
        if abs(v-lo) < m: print(f'         ⚠ {n} 碰下界 {lo}')
        if abs(v-hi) < m: print(f'         ⚠ {n} 碰上界 {hi}')


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    alpha_cu = CONFIG.get('alpha_cu', 0.00393)
    T_ref = CONFIG.get('T_ref', 25.0)

    data = load_multiple(sys.argv[1:], CONFIG)
    segments = extract_lock_segments(data, CONFIG)
    print(f'{len(segments)} 个堵转段')

    # ── 基线: 单节点联合拟合 ──
    print('\n[基线] 单节点联合拟合:')
    x1, _ = fit_phase1_segments(segments, CONFIG, verbose=True)
    base_sims = [p1_sim_segment(s, *x1, alpha_cu, T_ref) for s in segments]
    base_rmse = [seg_rmse(sim, s['T1']) for sim, s in zip(base_sims, segments)]
    print('  ' + ', '.join(f'Seg{i+1} RMSE={r:.2f}°C' for i, r in enumerate(base_rmse)))

    # ── 热点模型: Seg1 单独 ──
    print('\n[热点模型] Seg 1 单独拟合:')
    xh1, _ = fit_hotspot([segments[0]], CONFIG)
    print_params(xh1, 'Seg1-own')

    # ── 热点模型: 全段联合 ──
    print('\n[热点模型] 全段联合拟合:')
    xhj, _ = fit_hotspot(segments, CONFIG)
    print_params(xhj, 'joint')

    # ── 汇总 ──
    print('\n' + '='*72)
    print(f"{'':14s}" + ''.join(f'{"Seg"+str(i+1)+" RMSE[°C]":>16s}' for i in range(len(segments))))
    print(f"{'单节点联合':14s}" + ''.join(f'{r:16.2f}' for r in base_rmse))
    hs_own = [hs_sim_segment(segments[0], xh1, alpha_cu, T_ref)[0]]
    r_own = seg_rmse(hs_own[0], segments[0]['T1'])
    print(f"{'热点 Seg1单独':14s}{r_own:16.2f}" + ''.join(f'{"-":>16s}' for _ in segments[1:]))
    hs_joint = [hs_sim_segment(s, xhj, alpha_cu, T_ref) for s in segments]
    r_joint = [seg_rmse(sim[0], s['T1']) for sim, s in zip(hs_joint, segments)]
    print(f"{'热点 联合':14s}" + ''.join(f'{r:16.2f}' for r in r_joint))

    # ── 图: 每段一行 ──
    n = len(segments)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3.6*n), squeeze=False)
    for i, seg in enumerate(segments):
        ax = axes[i][0]
        tt = seg['t'] - seg['t'][0]
        ax.plot(tt, seg['T1'], color='#6b7280', lw=1.4, label='T1 measured')
        ax.plot(tt, base_sims[i], color='#3477c9', lw=1.6, ls='--',
                label=f'single-node (RMSE={base_rmse[i]:.2f}°C)')
        Th, Tb = hs_joint[i]
        ax.plot(tt, Th, color='#d97706', lw=1.8, ls='-.',
                label=f'hotspot Th (RMSE={r_joint[i]:.2f}°C)')
        ax.plot(tt, Tb, color='#d97706', lw=1.0, ls=':', alpha=.7,
                label='hotspot Tb (bulk, unmeasured)')
        if i == 0:
            ax.plot(tt, hs_own[0], color='#0e7a5f', lw=1.2, ls=(0, (3, 1, 1, 1)),
                    alpha=.9, label=f'hotspot Seg1-own (RMSE={r_own:.2f}°C)')
        ax.set_title(f"Seg {i+1}: t={seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s, "
                     f"I_rms={np.sqrt(np.mean(seg['I2'])):.1f}A", fontsize=10)
        ax.set_ylabel('T (°C)'); ax.grid(True, alpha=.25)
        ax.legend(fontsize=8.5, loc='best')
    axes[-1][0].set_xlabel('t within segment (s)')
    fig.suptitle('Locked-rotor: single-node vs hotspot sub-node winding model',
                 fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(OUTPUT_PLOT, dpi=150, bbox_inches='tight')
    print(f'\nPlot: {OUTPUT_PLOT}')


if __name__ == '__main__':
    main()
