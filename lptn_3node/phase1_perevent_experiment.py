"""
分事件损耗系数实验: 每个堵转加热事件独立的 a_e
==============================================
假设: 两次堵转的转子位置/占空比不同 → 电流集中的相不同, 传感器看到的
等效损耗系数不同。热容/热阻是几何决定的, 保持恒定。

变体 A: 单节点 + 分事件 a_e            (C1, R13, a_1..a_E)
变体 B: 热点两节点 + 分事件 a_e        (C_h, C_b, R_h, R13, γ, a_1..a_E)

事件 = 堵转段内 I² > 5A² 的连续区间 (间隔<10s 合并)。
非事件点 (I²小) 用 mean(a_e), 功率贡献可忽略。

用法: python phase1_perevent_experiment.py <CSV1> [CSV2] ...
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import (CONFIG, load_multiple, extract_lock_segments,
                              p1_sim_segment, fit_phase1_segments)

OUTPUT_PLOT = 'phase1_perevent_compare.png'
MAX_SUBSTEP = 25
I2_THRESH = 5.0
MERGE_GAP = 10.0  # s


def tag_events(segments):
    """给每段打事件标签 seg['eid'] (非事件=-1), 返回事件总数与摘要"""
    n_event = 0
    summaries = []
    for si, seg in enumerate(segments):
        t = seg['t']; on = seg['I2'] > I2_THRESH
        N = len(t)
        runs = []
        k = 0
        while k < N:
            if on[k]:
                j = k
                while j < N and on[j]: j += 1
                runs.append([k, j])  # [start, end)
                k = j
            else:
                k += 1
        # 合并近邻 run
        merged = []
        for r in runs:
            if merged and t[r[0]] - t[merged[-1][1]-1] < MERGE_GAP:
                merged[-1][1] = r[1]
            else:
                merged.append(r)
        eid = np.full(N, -1, int)
        for r in merged:
            if r[1] - r[0] >= 5:
                eid[r[0]:r[1]] = n_event
                summaries.append(
                    f'  Event {n_event+1}: Seg{si+1} t={t[r[0]]:.0f}~{t[r[1]-1]:.0f}s, '
                    f'I2_mean={np.mean(seg["I2"][r[0]:r[1]]):.1f}A², '
                    f'T1: {seg["T1"][r[0]]:.0f}->{seg["T1"][r[1]-1]:.0f}°C')
                n_event += 1
        seg['eid'] = eid
    return n_event, summaries


def a_points(seg, a_es):
    """每个采样点的损耗系数"""
    a_off = float(np.mean(a_es))
    return np.where(seg['eid'] >= 0, np.asarray(a_es)[seg['eid']], a_off)


# ── 变体 A: 单节点 + 分事件 a ──
def simA(seg, C1, R13, a_es, alpha, Tref):
    a_pt = a_points(seg, a_es)
    I2 = seg['I2']; T3 = seg['T3']; dt = seg['dt']
    N = len(seg['t'])
    T = np.empty(N); T[0] = seg['T1'][0]
    iC = 1.0/C1; iR = 1.0/R13
    for k in range(N-1):
        P = a_pt[k] * (1.0 + alpha*(T[k]-Tref)) * I2[k]
        T[k+1] = T[k] + (P - (T[k]-T3[k])*iR) * iC * dt[k]
    return T


def costA(p, segments, nE, alpha, Tref):
    C1, R13 = p[0], p[1]; a_es = p[2:2+nE]
    sq = dsq = 0.0; n = dn = 0
    try:
        for seg in segments:
            T = simA(seg, C1, R13, a_es, alpha, Tref)
            if np.any(np.isnan(T)) or np.any(np.abs(T) > 500): return 1e6
            sq += np.sum((T - seg['T1'])**2); n += len(T)
            dT = np.diff(T)/seg['dt_safe']
            dsq += np.sum((dT - seg['dT1_meas'])**2); dn += len(dT)
        return sq/n + 0.5*dsq/dn
    except Exception:
        return 1e6


# ── 变体 B: 热点两节点 + 分事件 a ──
def simB(seg, p, nE, alpha, Tref):
    C_h, C_b, R_h, R13, gamma = p[:5]; a_es = p[5:5+nE]
    a_pt = a_points(seg, a_es)
    I2 = seg['I2']; T3 = seg['T3']; dt_arr = seg['dt']
    N = len(seg['t'])
    Th = np.empty(N); Tb = np.empty(N)
    Th[0] = Tb[0] = seg['T1'][0]
    iCh = 1.0/C_h; iCb = 1.0/C_b; iRh = 1.0/R_h; iR13 = 1.0/R13
    tau_fast = C_h * R_h
    for k in range(N-1):
        dt = dt_arr[k]
        m = int(dt / (0.4*tau_fast)) + 1
        if m > MAX_SUBSTEP: return None, None
        h = dt/m
        th, tb = Th[k], Tb[k]
        aI2 = a_pt[k] * I2[k]
        for _ in range(m):
            Ph = gamma * aI2 * (1.0 + alpha*(th-Tref))
            Pb = (1.0-gamma) * aI2 * (1.0 + alpha*(tb-Tref))
            q = (th-tb)*iRh
            th += (Ph - q)*iCh*h
            tb += (Pb + q - (tb-T3[k])*iR13)*iCb*h
        Th[k+1] = th; Tb[k+1] = tb
    return Th, Tb


def costB(p, segments, nE, alpha, Tref):
    sq = dsq = 0.0; n = dn = 0
    try:
        for seg in segments:
            Th, _ = simB(seg, p, nE, alpha, Tref)
            if Th is None or np.any(np.isnan(Th)) or np.any(np.abs(Th) > 500):
                return 1e6
            sq += np.sum((Th - seg['T1'])**2); n += len(Th)
            dT = np.diff(Th)/seg['dt_safe']
            dsq += np.sum((dT - seg['dT1_meas'])**2); dn += len(dT)
        return sq/n + 0.5*dsq/dn
    except Exception:
        return 1e6


def fit(cost, bounds, segments, nE, alpha, Tref, seeds=(42, 123), maxiter=250):
    best_c, best_x = 1e10, None
    for seed in seeds:
        print(f'    DE (seed={seed})...', end=' ', flush=True)
        t0 = time.time()
        r = differential_evolution(cost, bounds, args=(segments, nE, alpha, Tref),
            seed=seed, maxiter=maxiter, popsize=15, tol=1e-9, polish=False,
            workers=1, mutation=(0.5, 1.5), recombination=0.9)
        r2 = minimize(cost, r.x, args=(segments, nE, alpha, Tref),
            method='L-BFGS-B', bounds=bounds, options={'maxiter': 400, 'ftol': 1e-14})
        print(f'cost={r2.fun:.4f} ({time.time()-t0:.0f}s)')
        if r2.fun < best_c:
            best_c, best_x = r2.fun, r2.x
    return best_x, best_c


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    alpha = CONFIG.get('alpha_cu', 0.00393)
    Tref = CONFIG.get('T_ref', 25.0)

    data = load_multiple(sys.argv[1:], CONFIG)
    segments = extract_lock_segments(data, CONFIG)
    nE, summaries = tag_events(segments)
    print(f'{len(segments)} 个堵转段, {nE} 个加热事件:')
    for s in summaries: print(s)

    # 基线: 单节点恒定 a (联合)
    print('\n[基线] 单节点 恒定a:')
    x1, _ = fit_phase1_segments(segments, CONFIG, verbose=True)
    base = [p1_sim_segment(s, *x1, alpha, Tref) for s in segments]
    rb = [np.sqrt(np.mean((sim - s['T1'])**2)) for sim, s in zip(base, segments)]

    # 变体 A
    print('\n[变体A] 单节点 + 分事件 a_e:')
    bndA = [(5, 500), (0.05, 5)] + [(0.05, 20)]*nE
    xA, _ = fit(costA, bndA, segments, nE, alpha, Tref)
    print(f'  C1={xA[0]:.1f} J/K, R13={xA[1]:.3f} K/W, tau={xA[0]*xA[1]:.1f}s')
    print('  a_e = ' + ', '.join(f'{v:.3f}' for v in xA[2:2+nE]) +
          f'  (比值 max/min = {max(xA[2:2+nE])/min(xA[2:2+nE]):.2f}x)')
    simsA = [simA(s, xA[0], xA[1], xA[2:2+nE], alpha, Tref) for s in segments]
    rA = [np.sqrt(np.mean((sim - s['T1'])**2)) for sim, s in zip(simsA, segments)]

    # 变体 B
    print('\n[变体B] 热点两节点 + 分事件 a_e:')
    bndB = [(1, 200), (20, 800), (0.05, 5), (0.05, 10), (0.05, 1.0)] + [(0.05, 20)]*nE
    xB, _ = fit(costB, bndB, segments, nE, alpha, Tref)
    print(f'  C_h={xB[0]:.1f}, C_b={xB[1]:.1f} J/K, R_h={xB[2]:.3f}, R13={xB[3]:.3f} K/W, γ={xB[4]:.2f}')
    print(f'  tau_hot={xB[0]*xB[2]:.1f}s, tau_bulk={xB[1]*xB[3]:.1f}s')
    print('  a_e = ' + ', '.join(f'{v:.3f}' for v in xB[5:5+nE]) +
          f'  (比值 max/min = {max(xB[5:5+nE])/min(xB[5:5+nE]):.2f}x)')
    simsB = [simB(s, xB, nE, alpha, Tref) for s in segments]
    rB = [np.sqrt(np.mean((sim[0] - s['T1'])**2)) for sim, s in zip(simsB, segments)]

    # 汇总
    print('\n' + '='*70)
    print(f"{'模型':24s}" + ''.join(f'{"Seg"+str(i+1)+" RMSE[°C]":>14s}' for i in range(len(segments))))
    print(f"{'单节点 恒定a (基线)':24s}" + ''.join(f'{r:14.2f}' for r in rb))
    print(f"{'A: 单节点+分事件a':24s}" + ''.join(f'{r:14.2f}' for r in rA))
    print(f"{'B: 热点+分事件a':24s}" + ''.join(f'{r:14.2f}' for r in rB))

    # 图
    n = len(segments)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3.6*n), squeeze=False)
    for i, seg in enumerate(segments):
        ax = axes[i][0]
        tt = seg['t'] - seg['t'][0]
        ax.plot(tt, seg['T1'], color='#6b7280', lw=1.4, label='T1 measured')
        ax.plot(tt, base[i], color='#3477c9', lw=1.5, ls='--',
                label=f'baseline const-a (RMSE={rb[i]:.2f}°C)')
        ax.plot(tt, simsA[i], color='#0e7a5f', lw=1.7, ls='-.',
                label=f'A: per-event a (RMSE={rA[i]:.2f}°C)')
        ax.plot(tt, simsB[i][0], color='#d97706', lw=1.9, ls=(0, (4, 1)),
                label=f'B: hotspot + per-event a (RMSE={rB[i]:.2f}°C)')
        ax.plot(tt, simsB[i][1], color='#d97706', lw=0.9, ls=':', alpha=.6,
                label='B: bulk Tb (unmeasured)')
        # 事件区间底纹
        ev = seg['eid'] >= 0
        if ev.any():
            yl = ax.get_ylim()
            ax.fill_between(tt, yl[0], yl[1], where=ev, alpha=.07, color='red',
                            label='heating events')
            ax.set_ylim(yl)
        ax.set_title(f"Seg {i+1}: t={seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s", fontsize=10)
        ax.set_ylabel('T (°C)'); ax.grid(True, alpha=.25)
        ax.legend(fontsize=8.5, loc='best')
    axes[-1][0].set_xlabel('t within segment (s)')
    fig.suptitle('Locked-rotor: per-event loss coefficient hypothesis',
                 fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(OUTPUT_PLOT, dpi=150, bbox_inches='tight')
    print(f'\nPlot: {OUTPUT_PLOT}')


if __name__ == '__main__':
    main()
