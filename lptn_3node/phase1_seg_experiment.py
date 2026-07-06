"""
Phase-1 堵转分段辨识对比实验
============================
目的: 区分堵转段拟合不准的原因
  1) 联合辨识: 所有堵转段共用一组 (C1, a_eff, R13) — 分段独立仿真(已修bug)
  2) 分段辨识: 每个有效堵转段单独辨识一组 (C1, a_eff, R13)
判读:
  - 各段参数差异大  → 模型结构问题(局部热点/铜损系数不恒定), 单组参数本就拟合不了
  - 各段参数接近    → 原问题主要来自旧版拼接bug/代价权重, 修bug即可解决

用法: python phase1_seg_experiment.py <CSV1> [CSV2] ...
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import (CONFIG, load_multiple, extract_lock_segments,
                              p1_sim_segment, fit_phase1_segments)

MIN_PTS = 30       # 单段辨识最少点数
MIN_RANGE = 5.0    # 单段辨识最小温升范围 [°C], 低于此激励不足

P1_BOUNDS = [(5, 500), (0.1, 20), (0.1, 50)]  # 与 fit_phase1_segments 内一致
OUTPUT_PLOT = 'phase1_seg_compare.png'


def bound_flags(x):
    flags = []
    for v, (lo, hi), name in zip(x, P1_BOUNDS, ['C1', 'a_eff', 'R13']):
        m = 1e-3 * (hi - lo)
        if abs(v - lo) < m: flags.append(f'{name}碰下界')
        elif abs(v - hi) < m: flags.append(f'{name}碰上界')
    return ','.join(flags) if flags else '-'


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)

    alpha_cu = CONFIG.get('alpha_cu', 0.00393)
    T_ref = CONFIG.get('T_ref', 25.0)

    print('Loading:', sys.argv[1:])
    data = load_multiple(sys.argv[1:], CONFIG)
    segments = extract_lock_segments(data, CONFIG)
    print(f'共 {len(segments)} 个堵转段 (>= 10 pts):')
    for i, s in enumerate(segments):
        print(f"  [Seg {i+1}] t={s['t'][0]:.0f}~{s['t'][-1]:.0f}s ({len(s['t'])} pts), "
              f"I_rms={np.sqrt(np.mean(s['I2'])):.2f}A, I2_mean={np.mean(s['I2']):.1f}A², "
              f"T1: {s['T1'][0]:.1f}->{s['T1'][-1]:.1f}°C (range={s['T1'].max()-s['T1'].min():.1f}°C)")

    # ── 1) 联合辨识 ──
    print('\n[联合辨识] 所有段共用一组参数:')
    x_joint, c_joint = fit_phase1_segments(segments, CONFIG)
    print(f'  C1={x_joint[0]:.2f} J/K, a_eff={x_joint[1]:.4f} W/A², R13={x_joint[2]:.4f} K/W, '
          f'tau={x_joint[0]*x_joint[2]:.1f}s  [{bound_flags(x_joint)}]')

    # ── 2) 分段辨识 ──
    analyzed = []   # (idx, seg, x_own, rmse_own, rmse_joint)
    for i, seg in enumerate(segments):
        rng = seg['T1'].max() - seg['T1'].min()
        T1_joint = p1_sim_segment(seg, *x_joint, alpha_cu, T_ref)
        rmse_joint = np.sqrt(np.mean((T1_joint - seg['T1'])**2))
        if len(seg['t']) < MIN_PTS or rng < MIN_RANGE:
            print(f'\n[Seg {i+1}] 跳过单段辨识 (pts={len(seg["t"])}, range={rng:.1f}°C 激励不足); '
                  f'联合参数 RMSE={rmse_joint:.2f}°C')
            continue
        print(f'\n[Seg {i+1}] 单段辨识:')
        x_own, _ = fit_phase1_segments([seg], CONFIG, verbose=True)
        T1_own = p1_sim_segment(seg, *x_own, alpha_cu, T_ref)
        rmse_own = np.sqrt(np.mean((T1_own - seg['T1'])**2))
        print(f'  C1={x_own[0]:.2f} J/K, a_eff={x_own[1]:.4f} W/A², R13={x_own[2]:.4f} K/W, '
              f'tau={x_own[0]*x_own[2]:.1f}s  [{bound_flags(x_own)}]')
        print(f'  RMSE: 单段={rmse_own:.2f}°C vs 联合={rmse_joint:.2f}°C')
        analyzed.append((i, seg, x_own, T1_own, T1_joint, rmse_own, rmse_joint))

    # ── 对比表 ──
    print('\n' + '=' * 88)
    print('参数对比 (联合 vs 各段独立)')
    print('=' * 88)
    hdr = f"{'':10s} {'C1[J/K]':>9s} {'a_eff[W/A²]':>12s} {'R13[K/W]':>9s} {'tau[s]':>7s} " \
          f"{'RMSE own':>9s} {'RMSE joint':>10s}  碰界"
    print(hdr)
    print(f"{'联合':10s} {x_joint[0]:9.2f} {x_joint[1]:12.4f} {x_joint[2]:9.4f} "
          f"{x_joint[0]*x_joint[2]:7.1f} {'-':>9s} {'-':>10s}  {bound_flags(x_joint)}")
    for i, seg, x, _, _, ro, rj in analyzed:
        print(f"{'Seg '+str(i+1):10s} {x[0]:9.2f} {x[1]:12.4f} {x[2]:9.4f} "
              f"{x[0]*x[2]:7.1f} {ro:9.2f} {rj:10.2f}  {bound_flags(x)}")

    if len(analyzed) >= 2:
        c1s = np.array([x[0] for _, _, x, _, _, _, _ in analyzed])
        aes = np.array([x[1] for _, _, x, _, _, _, _ in analyzed])
        r13 = np.array([x[2] for _, _, x, _, _, _, _ in analyzed])
        print('\n段间离散度 (max/min):')
        print(f'  C1: {c1s.max()/c1s.min():.2f}x   a_eff: {aes.max()/aes.min():.2f}x   '
              f'R13: {r13.max()/r13.min():.2f}x')

    # ── 诊断图: 每段一行, 实测 vs 联合 vs 单段 ──
    n = len(analyzed)
    if n:
        fig, axes = plt.subplots(n, 1, figsize=(12, 3.2 * n), squeeze=False)
        for row, (i, seg, x, T1_own, T1_joint, ro, rj) in enumerate(analyzed):
            ax = axes[row][0]
            tt = seg['t'] - seg['t'][0]
            ax.plot(tt, seg['T1'], color='#6b7280', lw=1.4, label='T1 measured')
            ax.plot(tt, T1_joint, color='#3477c9', lw=1.8, ls='--',
                    label=f'joint fit (RMSE={rj:.2f}°C)')
            ax.plot(tt, T1_own, color='#d97706', lw=1.8, ls='-.',
                    label=f'own fit (RMSE={ro:.2f}°C)')
            ax.set_title(f"Seg {i+1}: t={seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s, "
                         f"I_rms={np.sqrt(np.mean(seg['I2'])):.1f}A | "
                         f"own: C1={x[0]:.0f}, a_eff={x[1]:.2f}, R13={x[2]:.2f} | "
                         f"joint: C1={x_joint[0]:.0f}, a_eff={x_joint[1]:.2f}, R13={x_joint[2]:.2f}",
                         fontsize=10)
            ax.set_ylabel('T1 (°C)')
            ax.grid(True, alpha=.25)
            ax.legend(fontsize=9, loc='best')
        axes[-1][0].set_xlabel('t within segment (s)')
        fig.suptitle('Phase-1 locked-rotor: joint vs per-segment identification',
                     fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(OUTPUT_PLOT, dpi=150, bbox_inches='tight')
        print(f'\nPlot: {OUTPUT_PLOT}')


if __name__ == '__main__':
    main()
