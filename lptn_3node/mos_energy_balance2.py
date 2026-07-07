"""
MOS 节点能量平衡 v2 —— 每事件"总驱动功率"(不拆 a2/f, 避免堵转内简并)
====================================================================
堵转段内电流近似恒定 → 驱动损耗近似恒定, 故直接把每事件驱动功率当一个
自由常数 P2_e 辨识, 再比 P2/Iq² 的 E1/E2 比值。
  dT2/dt = [P2_e - (T2-T1)/R12 - (T2-T3)/R23 - (T2-Tamb)/R2a]/C2
判据: (P2_E1/I²_E1)/(P2_E2/I²_E2)
  ≈2×  → 逆变器也多耗 → 支持"额外电流(id)"
  ≈1×  → 驱动损耗随 Iq² 正常 → 绕组额外损耗不伴随额外逆变器电流
用法: python mos_energy_balance2.py <CSV>
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lptn_3node_final import CONFIG, load_multiple

LOCK = (25, 703); EVENTS = [(25, 70), (132, 375)]
OUT = 'mos_energy_balance2.png'
# C2, R12, R23, R2a, P2_E1, P2_E2, P2_off
BND = [(1, 400), (0.5, 1000), (0.1, 500), (0.5, 500), (0, 50), (0, 50), (0, 50)]


def build(csv):
    d = load_multiple([csv], CONFIG)
    d = d[(d['t'] >= LOCK[0]) & (d['t'] <= LOCK[1])].reset_index(drop=True)
    t = d['t'].values
    seg = dict(t=t, T1=d['T1'].values, T2=d['T2'].values, T3=d['T3'].values,
               Tamb=d['Tamb'].values, I2=d['I'].values**2, dt=np.diff(t))
    seg['dt_safe'] = np.where(seg['dt'] > 0, seg['dt'], 1.0)
    seg['dT2_meas'] = np.diff(seg['T2'])/seg['dt_safe']
    eid = np.full(len(t), -1, int)
    for e, (lo, hi) in enumerate(EVENTS):
        eid[(t >= lo) & (t <= hi)] = e
    seg['eid'] = eid
    return seg


def P2_pt(seg, P2s, P2off):
    return np.where(seg['eid'] >= 0, np.asarray(P2s)[np.maximum(seg['eid'], 0)], P2off)


def sim(seg, p):
    C2, R12, R23, R2a = p[:4]; P2s = p[4:6]; P2off = p[6]
    pp = P2_pt(seg, P2s, P2off)
    T1, T3, Ta = seg['T1'], seg['T3'], seg['Tamb']; dt = seg['dt']; N = len(seg['t'])
    T2 = np.empty(N); T2[0] = seg['T2'][0]
    iC = 1/C2; i12 = 1/R12; i23 = 1/R23; i2a = 1/R2a
    for k in range(N-1):
        dT2 = (pp[k] - (T2[k]-T1[k])*i12 - (T2[k]-T3[k])*i23 - (T2[k]-Ta[k])*i2a)*iC
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
        r = differential_evolution(cost, BND, args=(seg,), seed=s, maxiter=300, popsize=18,
            tol=1e-10, polish=False, workers=1, mutation=(0.5, 1.5), recombination=0.9)
        r2 = minimize(cost, r.x, args=(seg,), method='L-BFGS-B', bounds=BND,
            options={'maxiter': 800, 'ftol': 1e-15})
        if r2.fun < bc: bc, bx = r2.fun, r2.x
    return bx


def erm(s_, seg, e):
    m = seg['eid'] == e; return np.sqrt(np.mean((s_[m]-seg['T2'][m])**2))


def main():
    if len(sys.argv) < 2: print(__doc__); sys.exit(1)
    seg = build(sys.argv[1])
    I2 = [seg['I2'][seg['eid'] == e].mean() for e in range(2)]
    print(f"堵转段 {len(seg['t'])}pts;  I²: E1={I2[0]:.0f}  E2={I2[1]:.0f}  (Iq²比={I2[0]/I2[1]:.2f})")
    for e, (lo, hi) in enumerate(EVENTS):
        m = seg['eid'] == e
        print(f"  E{e+1}: T2(MOS) {seg['T2'][m][0]:.0f}->{seg['T2'][m][-1]:.0f}°C")

    print("\n[MOS 孤立辨识, 每事件总驱动功率]"); t0 = time.time()
    x = fit(seg); T2s = sim(seg, x)
    C2, R12, R23, R2a = x[:4]; P2 = x[4:6]; P2off = x[6]
    print(f"  C2={C2:.1f} R12={R12:.1f} R23={R23:.1f} R2a={R2a:.1f}  P2off={P2off:.2f}W  ({time.time()-t0:.0f}s)")
    print(f"  τ2≈C2·(R23∥R2a∥R12)={C2/(1/R23+1/R2a+1/R12):.0f}s")
    print(f"  拟合RMSE: E1={erm(T2s,seg,0):.2f} E2={erm(T2s,seg,1):.2f} 全段={np.sqrt(np.mean((T2s-seg['T2'])**2)):.2f}°C")
    print(f"  驱动功率 P2: E1={P2[0]:.2f}W  E2={P2[1]:.2f}W")
    r_p = P2[0]/P2[1]; r_i = I2[0]/I2[1]
    print(f"  P2/Iq²: E1={P2[0]/I2[0]:.4f}  E2={P2[1]/I2[1]:.4f}")
    print(f"  ★ MOS 归一比值 (P2/Iq²)_E1/(P2/Iq²)_E2 = {(P2[0]/I2[0])/(P2[1]/I2[1]):.2f}×   (绕组: ~2.06×)")
    hit = [n for v,(lo,hi),n in zip(x,BND,['C2','R12','R23','R2a','P2_E1','P2_E2','P2off']) if abs(v-lo)<1e-3*(hi-lo) or abs(v-hi)<1e-3*(hi-lo)]
    print(f"  碰界: {hit if hit else '无 ✓'}")

    ratio = (P2[0]/I2[0])/(P2[1]/I2[1])
    print("\n  判据: ", end="")
    if ratio >= 1.6: print(f"MOS 也需 ~{ratio:.1f}× → 逆变器同样多耗电流 → 支持'额外电流(id)' ✓")
    elif ratio <= 1.3: print(f"MOS 归一功率基本不变({ratio:.2f}×) → 额外绕组损耗不伴随额外逆变器电流 → 反对 id, 指向绕组侧/热路")
    else: print(f"MOS 比值 {ratio:.2f}× 居中, 证据不强")

    t = seg['t'] - seg['t'][0]
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(t, seg['T2'], color='#6b7280', lw=1.6, label='T2 MOS measured')
    ax.plot(t, T2s, color='#d97706', lw=1.8, ls='--', label=f'MOS fit (RMSE={np.sqrt(np.mean((T2s-seg["T2"])**2)):.2f}°C)')
    ax.plot(t, seg['T1'], color='#c0392b', lw=1.0, alpha=.4, label='T1 winding (input)')
    ax.plot(t, seg['Tamb'], color='#888', lw=1.0, alpha=.4, label='Tamb (input)')
    for e, (lo, hi) in enumerate(EVENTS):
        ax.axvspan(lo-seg['t'][0], hi-seg['t'][0], alpha=.08, color=['#d97706','#0e7a5f'][e])
    ax.set_xlabel('t within segment (s)'); ax.set_ylabel('T (°C)')
    ax.set_title(f'MOS node isolated | driver-loss/Iq² ratio E1/E2 = {ratio:.2f}× (winding ~2.06×)', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(True, alpha=.25)
    plt.tight_layout(); plt.savefig(OUT, dpi=150, bbox_inches='tight')
    print(f'\nPlot: {OUT}')


if __name__ == '__main__':
    main()
