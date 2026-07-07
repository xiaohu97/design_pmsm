"""
两节点(局部热点+本体) + 分事件损耗系数 —— 闭环确认
====================================================
对照两种结构, 都用分事件 a_e (每个堵转加热事件一个系数):
  SN: 单节点         (C1, R13, a_1..a_E)
  TN: 两节点(热点)   (C_h, C_b, R_h, R13, γ, a_1..a_E)
预期确认:
  单节点分事件 a 比值≈2.07x (结构缺陷+真实差 都压进系数)
  两节点分事件 a 比值应降到 ~1.2x (结构吃掉瞬态伪影, 只剩真实功率差)
  且冷启峰能顶到 ~99°C
用法: python phase1_twonode_perevent.py <CSV>
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import CONFIG, load_multiple, extract_lock_segments, p1_sim_segment

MAX_SUB = 40
I2_THRESH = 5.0
MERGE_GAP = 10.0
COLD = (25, 72)
STEADY = (330, 375)
OUT = 'phase1_twonode_perevent.png'


def tag_events(seg):
    t = seg['t']; on = seg['I2'] > I2_THRESH; N = len(t)
    runs = []; k = 0
    while k < N:
        if on[k]:
            j = k
            while j < N and on[j]: j += 1
            runs.append([k, j]); k = j
        else: k += 1
    merged = []
    for r in runs:
        if merged and t[r[0]] - t[merged[-1][1]-1] < MERGE_GAP:
            merged[-1][1] = r[1]
        else: merged.append(list(r))
    eid = np.full(N, -1, int); ne = 0; info = []
    for r in merged:
        if r[1]-r[0] >= 5:
            eid[r[0]:r[1]] = ne
            info.append((t[r[0]], t[r[1]-1], np.mean(seg['I2'][r[0]:r[1]]),
                         seg['T1'][r[0]], seg['T1'][r[1]-1]))
            ne += 1
    seg['eid'] = eid
    return ne, info


def a_pt(seg, a_es):
    a_off = float(np.mean(a_es))
    return np.where(seg['eid'] >= 0, np.asarray(a_es)[seg['eid']], a_off)


# ── 单节点 + 分事件 a ──
def sn_sim(seg, C1, R13, a_es, alpha, Tref):
    ap = a_pt(seg, a_es); I2 = seg['I2']; T3 = seg['T3']; dt = seg['dt']; N = len(seg['t'])
    T = np.empty(N); T[0] = seg['T1'][0]; iC = 1/C1; iR = 1/R13
    for k in range(N-1):
        P = ap[k]*(1+alpha*(T[k]-Tref))*I2[k]
        T[k+1] = T[k] + (P - (T[k]-T3[k])*iR)*iC*dt[k]
    return T

def sn_cost(p, seg, nE, alpha, Tref):
    C1, R13 = p[0], p[1]; a_es = p[2:2+nE]
    try:
        T = sn_sim(seg, C1, R13, a_es, alpha, Tref)
        if np.any(np.isnan(T)) or np.any(np.abs(T) > 500): return 1e6
        dT = np.diff(T)/seg['dt_safe']
        return np.mean((T-seg['T1'])**2) + 0.5*np.mean((dT-seg['dT1_meas'])**2)
    except Exception: return 1e6


# ── 两节点 + 分事件 a ──
def tn_sim(seg, p, nE, alpha, Tref):
    C_h, C_b, R_h, R13, g = p[:5]; a_es = p[5:5+nE]
    ap = a_pt(seg, a_es); I2 = seg['I2']; T3 = seg['T3']; dt_arr = seg['dt']; N = len(seg['t'])
    Th = np.empty(N); Tb = np.empty(N); Th[0] = Tb[0] = seg['T1'][0]
    iCh=1/C_h; iCb=1/C_b; iRh=1/R_h; iR13=1/R13; tau=C_h*R_h
    for k in range(N-1):
        dt = dt_arr[k]; m = int(dt/(0.4*tau))+1
        if m > MAX_SUB: return None, None
        h = dt/m; th, tb = Th[k], Tb[k]; aI2 = ap[k]*I2[k]
        for _ in range(m):
            Ph = g*aI2*(1+alpha*(th-Tref)); Pb = (1-g)*aI2*(1+alpha*(tb-Tref)); q = (th-tb)*iRh
            th += (Ph-q)*iCh*h; tb += (Pb+q-(tb-T3[k])*iR13)*iCb*h
        Th[k+1] = th; Tb[k+1] = tb
    return Th, Tb

def tn_cost(p, seg, nE, alpha, Tref):
    try:
        Th, _ = tn_sim(seg, p, nE, alpha, Tref)
        if Th is None or np.any(np.isnan(Th)) or np.any(np.abs(Th) > 500): return 1e6
        dT = np.diff(Th)/seg['dt_safe']
        return np.mean((Th-seg['T1'])**2) + 0.5*np.mean((dT-seg['dT1_meas'])**2)
    except Exception: return 1e6


def fit(cost, bounds, seg, nE, alpha, Tref, seeds=(42, 123), mi=220):
    bc, bx = 1e10, None
    for s in seeds:
        r = differential_evolution(cost, bounds, args=(seg, nE, alpha, Tref), seed=s,
            maxiter=mi, popsize=15, tol=1e-9, polish=False, workers=1,
            mutation=(0.5, 1.5), recombination=0.9)
        r2 = minimize(cost, r.x, args=(seg, nE, alpha, Tref), method='L-BFGS-B',
            bounds=bounds, options={'maxiter': 500, 'ftol': 1e-14})
        if r2.fun < bc: bc, bx = r2.fun, r2.x
    return bx


def peak(sim, seg, win):
    m = (seg['t'] >= win[0]) & (seg['t'] < win[1]); return sim[m].max() if m.any() else np.nan
def rmse(sim, seg, win=None):
    if win is None: return np.sqrt(np.mean((sim-seg['T1'])**2))
    m = (seg['t'] >= win[0]) & (seg['t'] < win[1]); return np.sqrt(np.mean((sim[m]-seg['T1'][m])**2))
def ev_rmse(sim, seg, e):
    m = seg['eid'] == e; return np.sqrt(np.mean((sim[m]-seg['T1'][m])**2))


def main():
    if len(sys.argv) < 2: print(__doc__); sys.exit(1)
    alpha = CONFIG.get('alpha_cu', 0.00393); Tref = CONFIG.get('T_ref', 25.0)
    data = load_multiple([sys.argv[1]], CONFIG)
    seg = max(extract_lock_segments(data, CONFIG), key=lambda s: len(s['t']))
    nE, info = tag_events(seg)
    print(f"堵转段 t={seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s, {nE} 个加热事件:")
    for i, (t0, t1, i2, T0, T1) in enumerate(info):
        print(f"  Event {i+1}: t={t0:.0f}~{t1:.0f}s, I²={i2:.0f}, T1 {T0:.0f}->{T1:.0f}°C")
    print(f"实测冷启峰 = {peak(seg['T1'], seg, COLD):.0f}°C\n")

    # 单节点 + 分事件
    print("[单节点 + 分事件 a]"); t0 = time.time()
    bnd_sn = [(5, 500), (0.05, 10)] + [(0.05, 15)]*nE
    xs = fit(sn_cost, bnd_sn, seg, nE, alpha, Tref)
    Ts = sn_sim(seg, xs[0], xs[1], xs[2:2+nE], alpha, Tref)
    a_sn = xs[2:2+nE]
    print(f"  C1={xs[0]:.1f} R13={xs[1]:.3f}  a_e={np.round(a_sn,2)}  比值={max(a_sn)/min(a_sn):.2f}x  ({time.time()-t0:.0f}s)")
    print(f"  冷启峰={peak(Ts,seg,COLD):.0f}°C  RMSE: " + " ".join(f"E{i+1}={ev_rmse(Ts,seg,i):.2f}" for i in range(nE)) +
          f"  准稳态={rmse(Ts,seg,STEADY):.2f} 全段={rmse(Ts,seg):.2f}\n")

    # 两节点 + 分事件
    print("[两节点 + 分事件 a]"); t0 = time.time()
    bnd_tn = [(1, 150), (20, 800), (0.02, 5), (0.05, 10), (0.05, 1.0)] + [(0.05, 15)]*nE
    xt = fit(tn_cost, bnd_tn, seg, nE, alpha, Tref)
    Th, Tb = tn_sim(seg, xt, nE, alpha, Tref)
    a_tn = xt[5:5+nE]
    C_h, C_b, R_h, R13, g = xt[:5]
    print(f"  C_h={C_h:.1f} C_b={C_b:.1f} R_h={R_h:.3f} R13={R13:.3f} γ={g:.2f}  ({time.time()-t0:.0f}s)")
    print(f"  τ_hot={C_h*R_h:.1f}s τ_bulk={C_b*R13:.1f}s")
    print(f"  a_e={np.round(a_tn,2)}  比值={max(a_tn)/min(a_tn):.2f}x")
    print(f"  冷启峰={peak(Th,seg,COLD):.0f}°C  RMSE: " + " ".join(f"E{i+1}={ev_rmse(Th,seg,i):.2f}" for i in range(nE)) +
          f"  准稳态={rmse(Th,seg,STEADY):.2f} 全段={rmse(Th,seg):.2f}")
    for v,(lo,hi),n in zip(xt, bnd_tn, ['C_h','C_b','R_h','R13','γ']+[f'a{i+1}' for i in range(nE)]):
        if abs(v-lo) < 1e-3*(hi-lo): print(f"    ⚠ {n} 碰下界 {lo}")
        if abs(v-hi) < 1e-3*(hi-lo): print(f"    ⚠ {n} 碰上界 {hi}")

    # 汇总
    print("\n" + "="*60)
    print(f"{'':22s}{'a比值':>8s}{'冷启峰':>8s}{'E1 RMSE':>9s}{'全段RMSE':>9s}")
    print(f"{'单节点+分事件':22s}{max(a_sn)/min(a_sn):7.2f}x{peak(Ts,seg,COLD):7.0f}°{ev_rmse(Ts,seg,0):9.2f}{rmse(Ts,seg):9.2f}")
    print(f"{'两节点+分事件':22s}{max(a_tn)/min(a_tn):7.2f}x{peak(Th,seg,COLD):7.0f}°{ev_rmse(Th,seg,0):9.2f}{rmse(Th,seg):9.2f}")

    # 图
    t = seg['t'] - seg['t'][0]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 5.2))
    for ax in (a1, a2):
        ax.plot(t, seg['T1'], color='#6b7280', lw=1.5, label='T1 measured (NTC)')
        ax.plot(t, Ts, color='#3477c9', lw=1.6, ls='--', label=f'single-node+per-event (peak={peak(Ts,seg,COLD):.0f}°C)')
        ax.plot(t, Th, color='#d97706', lw=2.0, ls='-.', label=f'two-node+per-event (peak={peak(Th,seg,COLD):.0f}°C)')
        ax.plot(t, Tb, color='#d97706', lw=1.0, ls=':', alpha=.6, label='two-node Tb (bulk)')
        ax.grid(True, alpha=.25); ax.set_xlabel('t within segment (s)'); ax.set_ylabel('T (°C)')
    a1.legend(fontsize=8.5, loc='best'); a1.set_title('Full locked segment', fontsize=11)
    a2.set_xlim(0, COLD[1]-seg['t'][0]+3); a2.set_title('Cold-start transient (E1) zoom', fontsize=11)
    fig.suptitle('Two-node hotspot + per-event coefficient (confirmation)', fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95]); plt.savefig(OUT, dpi=150, bbox_inches='tight')
    print(f'\nPlot: {OUT}')


if __name__ == '__main__':
    main()
