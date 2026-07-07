"""
MOS(驱动板)节点严格能量平衡 —— 独立验证"额外电流"假设
========================================================
把节点2(MOS)孤立辨识: 邻居 T1(绕组)、T3(壳)、Tamb 全用实测值当时变边界,
扣除绕组经 R12 灌入的热, 残差即真正的驱动损耗 P2。
  dT2/dt = [P2 - (T2-T1)/R12 - (T2-T3)/R23 - (T2-Tamb)/R2a] / C2
  P2 = a2_e·I² + f    (a2_e 每个堵转事件独立, f 共享待机)
判据: 若 MOS 的驱动损耗系数 a2_E1/a2_E2 ≈ 2× (与绕组 ~2.06× 同量级)
      → 逆变器也流过更多电流 → 独立坐实"额外电流(id)"假设
      若 ≈1× → 绕组额外损耗未伴随额外逆变器电流 → 反对 id 假设
用法: python mos_energy_balance.py <CSV>
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import CONFIG, load_multiple

LOCK = (25, 703)     # 堵转段
EVENTS = [(25, 70), (132, 375)]   # E1(冷启), E2(温热再加热)
OUT = 'mos_energy_balance.png'
# C2, R12, R23, R2a, f, a2_E1, a2_E2
BND = [(1, 400), (0.1, 100), (0.1, 100), (0.5, 80), (0, 5), (0.0, 5), (0.0, 5)]


def build(csv):
    d = load_multiple([csv], CONFIG)
    m = (d['t'] >= LOCK[0]) & (d['t'] <= LOCK[1])
    d = d[m].reset_index(drop=True)
    t = d['t'].values
    seg = dict(t=t, T1=d['T1'].values, T2=d['T2'].values, T3=d['T3'].values,
               Tamb=d['Tamb'].values, I2=d['I'].values**2, I=d['I'].values,
               dt=np.diff(t))
    seg['dt_safe'] = np.where(seg['dt'] > 0, seg['dt'], 1.0)
    seg['dT2_meas'] = np.diff(seg['T2']) / seg['dt_safe']
    # 事件标签
    eid = np.full(len(t), -1, int)
    for e, (lo, hi) in enumerate(EVENTS):
        eid[(t >= lo) & (t <= hi)] = e
    seg['eid'] = eid
    return seg


def a2_pt(seg, a2s):
    off = float(np.mean(a2s))
    return np.where(seg['eid'] >= 0, np.asarray(a2s)[seg['eid']], off)


def sim(seg, p):
    C2, R12, R23, R2a, f = p[:5]; a2s = p[5:7]
    ap = a2_pt(seg, a2s)
    T1, T3, Ta, I2 = seg['T1'], seg['T3'], seg['Tamb'], seg['I2']
    dt = seg['dt']; N = len(seg['t'])
    T2 = np.empty(N); T2[0] = seg['T2'][0]
    iC = 1/C2; i12 = 1/R12; i23 = 1/R23; i2a = 1/R2a
    for k in range(N-1):
        P2 = ap[k]*I2[k] + f
        dT2 = (P2 - (T2[k]-T1[k])*i12 - (T2[k]-T3[k])*i23 - (T2[k]-Ta[k])*i2a)*iC
        T2[k+1] = T2[k] + dT2*dt[k]
    return T2


def cost(p, seg):
    try:
        T2 = sim(seg, p)
        if np.any(np.isnan(T2)) or np.any(np.abs(T2) > 300): return 1e6
        dT2 = np.diff(T2)/seg['dt_safe']
        return np.mean((T2-seg['T2'])**2) + 0.3*np.mean((dT2-seg['dT2_meas'])**2)
    except Exception: return 1e6


def fit(seg, seeds=(42, 123, 7)):
    bc, bx = 1e10, None
    for s in seeds:
        r = differential_evolution(cost, BND, args=(seg,), seed=s, maxiter=300,
            popsize=16, tol=1e-10, polish=False, workers=1, mutation=(0.5, 1.5), recombination=0.9)
        r2 = minimize(cost, r.x, args=(seg,), method='L-BFGS-B', bounds=BND,
            options={'maxiter': 600, 'ftol': 1e-15})
        if r2.fun < bc: bc, bx = r2.fun, r2.x
    return bx


def ev_rmse(sim_, seg, e):
    m = seg['eid'] == e; return np.sqrt(np.mean((sim_[m]-seg['T2'][m])**2))


def main():
    if len(sys.argv) < 2: print(__doc__); sys.exit(1)
    seg = build(sys.argv[1])
    print(f"堵转段 t={seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s, {len(seg['t'])}pts")
    for e, (lo, hi) in enumerate(EVENTS):
        m = seg['eid'] == e
        print(f"  E{e+1}: t={lo}~{hi}s, I²={seg['I2'][m].mean():.0f}, "
              f"T2(MOS) {seg['T2'][m][0]:.0f}->{seg['T2'][m][-1]:.0f}°C, "
              f"T1(绕组) {seg['T1'][m][0]:.0f}->{seg['T1'][m][-1]:.0f}°C")

    print("\n[MOS 节点孤立辨识 + 分事件驱动损耗系数]"); t0 = time.time()
    x = fit(seg)
    T2s = sim(seg, x)
    C2, R12, R23, R2a, f = x[:5]; a2 = x[5:7]
    print(f"  C2={C2:.1f} J/K  R12={R12:.2f} R23={R23:.2f} R2a={R2a:.2f} K/W  f={f:.3f}W  ({time.time()-t0:.0f}s)")
    print(f"  拟合 RMSE: E1={ev_rmse(T2s,seg,0):.2f} E2={ev_rmse(T2s,seg,1):.2f} 全段={np.sqrt(np.mean((T2s-seg['T2'])**2)):.2f}°C")
    print(f"  驱动损耗系数 a2: E1={a2[0]:.3f}  E2={a2[1]:.3f}  W/A²")
    print(f"  ★ MOS 比值 a2_E1/a2_E2 = {a2[0]/a2[1]:.2f}×   (绕组对比: ~2.06×)")
    for v, (lo, hi), n in zip(x, BND, ['C2','R12','R23','R2a','f','a2_E1','a2_E2']):
        if abs(v-lo) < 1e-3*(hi-lo): print(f"    ⚠ {n} 碰下界 {lo}")
        if abs(v-hi) < 1e-3*(hi-lo): print(f"    ⚠ {n} 碰上界 {hi}")

    # 顺带算各事件平均"绕组→MOS 灌入热" 占比, 说明能量平衡扣除量
    for e, (lo, hi) in enumerate(EVENTS):
        m = seg['eid'] == e
        q_w2m = np.mean((seg['T1'][m]-T2s[m])/R12)   # 绕组流入MOS(正=流入)
        p2 = a2[e]*seg['I2'][m].mean()+f
        print(f"  E{e+1}: 驱动损耗 P2≈{p2:.2f}W, 绕组灌入≈{q_w2m:.2f}W (占比 {q_w2m/(p2+q_w2m)*100:.0f}%)")

    print("\n  判据: ", end="")
    r = a2[0]/a2[1]
    if r >= 1.6:
        print(f"MOS 也需 ~{r:.1f}× → 逆变器同样流过更多电流 → 独立坐实'额外电流(id/相电流)'假设 ✓")
    elif r <= 1.3:
        print(f"MOS 系数基本不变({r:.2f}×) → 额外绕组损耗未伴随额外逆变器电流 → 反对 id 假设")
    else:
        print(f"MOS 比值 {r:.2f}× 介于中间, 证据不强")

    # 图
    t = seg['t'] - seg['t'][0]
    fig, (a1, a2ax) = plt.subplots(2, 1, figsize=(13, 8), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 2]})
    a1.plot(t, seg['T2'], color='#6b7280', lw=1.5, label='T2 MOS measured')
    a1.plot(t, T2s, color='#d97706', lw=1.8, ls='--', label=f'MOS fit (RMSE={np.sqrt(np.mean((T2s-seg["T2"])**2)):.2f}°C)')
    a1.plot(t, seg['T1'], color='#c0392b', lw=1.0, alpha=.5, label='T1 winding (input)')
    a1.plot(t, seg['T3'], color='#3477c9', lw=1.0, alpha=.5, label='T3 shell (input)')
    for e, (lo, hi) in enumerate(EVENTS):
        a1.axvspan(lo-seg['t'][0], hi-seg['t'][0], alpha=.08, color=['#d97706','#0e7a5f'][e])
    a1.set_ylabel('T (°C)'); a1.legend(fontsize=9, loc='upper right'); a1.grid(True, alpha=.25)
    a1.set_title(f'MOS node isolated energy balance | a2_E1/a2_E2 = {a2[0]/a2[1]:.2f}× (winding ~2.06×)', fontweight='bold')
    # 下: 驱动损耗 P2 与 绕组灌入热
    ap = a2_pt(seg, a2); P2 = ap*seg['I2']+f
    q = (seg['T1']-T2s)/R12
    a2ax.plot(t, P2, color='#d97706', lw=1.4, label='driver loss P2 = a2·I²+f')
    a2ax.plot(t, q, color='#c0392b', lw=1.2, ls=':', label='winding→MOS heat (T1-T2)/R12')
    a2ax.set_ylabel('Power (W)'); a2ax.set_xlabel('t within segment (s)')
    a2ax.legend(fontsize=9, loc='upper right'); a2ax.grid(True, alpha=.25)
    plt.tight_layout(); plt.savefig(OUT, dpi=150, bbox_inches='tight')
    print(f'\nPlot: {OUT}')


if __name__ == '__main__':
    main()
