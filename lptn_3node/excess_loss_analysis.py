"""
额外损耗诊断: a_eff vs 抖动方差 / a_eff vs 电流幅值
====================================================
a_eff(窗) = P_req / [(1+α·ΔT)·Iq²],  P_req = C1·dT1/dt + (T1-T3)/R13
  - 热平衡用 1s 分箱数据 (load_multiple)
  - 抖动方差用 16Hz 原始 CSV (窗内 std)
判读:
  a_eff 随抖动↑而↑        → 抖动/谐波(B5) 是额外损耗来源
  a_eff 随 Iq²↓ 而↑ (反相) → 恒定电压误差型 (死区 B1)
  a_eff 与两者都无关        → 另有原因(未记录id / 冷热/首次-再次 等)
用法: python excess_loss_analysis.py <RAW_CSV>
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import CONFIG, load_multiple, extract_lock_segments, fit_phase1_segments

WIN = 8.0          # 窗长 [s]
I2_MIN = 20.0      # 只取有明显堵转电流的窗
OUTPUT = 'excess_loss_scatter.png'


def load_raw(path, cfg):
    df = pd.read_csv(path, encoding='utf-8-sig')
    df['ts'] = pd.to_datetime(df[cfg['time_col']], format=cfg['time_format'])
    df['t'] = (df['ts'] - df['ts'].iloc[0]).dt.total_seconds()
    return df


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    path = sys.argv[1]
    alpha = CONFIG.get('alpha_cu', 0.00393)
    Tref = CONFIG.get('T_ref', 25.0)

    data = load_multiple([path], CONFIG)          # 1s 分箱
    raw = load_raw(path, CONFIG)                  # 16Hz 原始
    segs = extract_lock_segments(data, CONFIG)

    # 全段联合辨识 C1, R13 (a_eff 用不到其 a)
    xj, _ = fit_phase1_segments(segs, CONFIG, verbose=False)
    C1, _, R13 = xj
    print(f'联合辨识: C1={C1:.1f} J/K, R13={R13:.4f} K/W (tau={C1*R13:.1f}s)\n')

    # 需要带 I/omega/T 列的分箱段: 用 lock_mask + 时间断点重切分箱 DataFrame
    dt_bin = CONFIG.get('resample_dt', 1)
    lock = np.abs(data['omega'].values) < CONFIG['lock_speed_thresh']
    idx = np.where(lock)[0]
    t_all = data['t'].values
    seg_slices, start = [], 0
    for k in range(1, len(idx)+1):
        brk = (k == len(idx) or idx[k] != idx[k-1]+1 or t_all[idx[k]]-t_all[idx[k-1]] > 3*dt_bin)
        if brk:
            if idx[k-1]-idx[start] >= 9:
                seg_slices.append(data.iloc[idx[start:k]])
            start = k

    rows = []
    for si, d in enumerate(seg_slices):
        t = d['t'].values; T1 = d['T1'].values; T3 = d['T3'].values; I2 = d['I'].values**2
        t0, t1 = t[0], t[-1]
        w = t0
        while w + WIN <= t1 + 1e-6:
            m = (t >= w) & (t < w+WIN)
            if np.sum(m) >= 5:
                tt, y = t[m], T1[m]
                slope = np.polyfit(tt-tt.mean(), y, 1)[0]        # dT1/dt 线性拟合(抗噪)
                P_req = C1*slope + np.mean((T1[m]-T3[m]))/R13
                I2m = np.mean(I2[m])
                dTmean = np.mean(T1[m]) - Tref
                denom = (1+alpha*dTmean)*I2m
                if I2m >= I2_MIN and denom > 0 and P_req > 0:
                    # 抖动: 原始 16Hz 窗内 std
                    r = raw[(raw['t'] >= w) & (raw['t'] < w+WIN)]
                    pos = r['Ch1_Pos_rad'].values
                    pos_dt = pos - np.polyval(np.polyfit(np.arange(len(pos)), pos, 1), np.arange(len(pos))) if len(pos) > 2 else pos*0
                    rows.append(dict(
                        seg=si+1, t=w, a_eff=P_req/denom, P_req=P_req, I2=I2m,
                        std_spd=r['转速(rpm)'].std(), std_iq=r['Ch1_Cur_A'].std(),
                        std_pos=np.std(pos_dt), std_tor=r['扭矩(Nm)'].std(),
                    ))
            w += WIN

    R = pd.DataFrame(rows)
    print(f'有效窗数: {len(R)} (窗长{WIN:.0f}s, I²≥{I2_MIN})\n')
    print(R[['seg','t','a_eff','I2','std_spd','std_iq','std_pos','std_tor']].round(3).to_string(index=False))

    # 相关性
    def corr(x, y):
        if len(x) < 4: return (np.nan, np.nan)
        pr = pearsonr(x, y)[0]; sr = spearmanr(x, y)[0]
        return pr, sr
    print('\n相关性 (a_eff vs 指标)   Pearson  Spearman')
    metrics = [('电流 I²', R['I2']), ('转速抖动 std_spd', R['std_spd']),
               ('Iq抖动 std_iq', R['std_iq']), ('位置抖动 std_pos', R['std_pos']),
               ('扭矩抖动 std_tor', R['std_tor'])]
    for name, col in metrics:
        pr, sr = corr(R['a_eff'].values, col.values)
        print(f'  {name:20s}  {pr:+.3f}   {sr:+.3f}')

    # 分段均值
    print('\n分段均值:')
    g = R.groupby('seg').agg(n=('a_eff','size'), a_eff=('a_eff','mean'), I2=('I2','mean'),
        std_spd=('std_spd','mean'), std_iq=('std_iq','mean'),
        std_pos=('std_pos','mean'), std_tor=('std_tor','mean'))
    print(g.round(3).to_string())

    # 散点图
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    pairs = [('I2', 'current amplitude Iq² (A²)'),
             ('std_spd', 'dither: TR10B speed std (rpm)'),
             ('std_iq', 'dither: Iq std (A)')]
    cmap = {1: '#d97706', 2: '#3477c9', 3: '#0e7a5f'}
    for ax, (key, xl) in zip(axes, pairs):
        for s in sorted(R['seg'].unique()):
            sub = R[R['seg'] == s]
            ax.scatter(sub[key], sub['a_eff'], s=45, alpha=.75,
                       color=cmap.get(s, '#888'), label=f'Seg {s}', edgecolor='white', linewidth=.5)
        pr = corr(R['a_eff'].values, R[key].values)[0]
        ax.set_xlabel(xl); ax.set_ylabel('a_eff = P_req/[(1+αΔT)Iq²]')
        ax.set_title(f'r = {pr:+.2f}', fontsize=11)
        ax.grid(True, alpha=.25); ax.legend(fontsize=8)
    fig.suptitle('Excess-loss coefficient vs current & dither (locked-rotor windows)',
                 fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(OUTPUT, dpi=150, bbox_inches='tight')
    print(f'\nPlot: {OUTPUT}')


if __name__ == '__main__':
    main()
