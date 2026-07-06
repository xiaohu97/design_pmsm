"""
三节点 LPTN 热网络参数辨识 (v3.1 - 多阶段 + 铜电阻温度系数)
==========================================================
Phase-1: 堵转段 → C1, a_eff, R13 (消除铜损/铁耗退化)
Phase-2: 全量数据 → 15参数 (C1 + a_lock 固定)

铜损温度系数: P_cu = a·(1 + α_cu·(ΔT))·I²
用法:  python lptn_3node_final.py <CSV1> [CSV2] ...
"""

import sys, time
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 配置
# ============================================================
CONFIG = {
    'resample_dt': 1,
    'col_map': {
        'T1': 'Ch1_Temp_C', 'T2': 'Ch1_MOS_C',
        'T3': 'P1.最高温', 'Tamb': '图像.最低温',
        'I': 'Ch1_Cur_A', 'omega': 'Ch1_Spd_rads',
    },
    'time_col': '系统时间',
    'time_format': '%Y/%m/%d %H:%M:%S.%f',

    'lock_speed_thresh': 0.5,
    'lock_weight': 5.0,
    'deriv_weight': 0.5,
    'rising_weight': 2.0,

    'bounds': {
        'C1': (5, 500), 'C2': (1, 400), 'C3': (5, 1500),
        'R12': (0.1, 100), 'R13': (0.1, 50), 'R23': (0.1, 100),
        'R2a': (0.5, 80), 'R3a': (0.1, 50),
        'a': (0.05, 15), 'b': (0, 2), 'c': (0, 0.5),
        'd': (0, 5), 'e': (0, 5), 'f': (0, 5),
        'a_lock': (0, 15), 'd_lock': (0, 5), 'omega0': (0.1, 10),
    },

    'target_weight_scale': {'T1': 5.0, 'T2': 1.0, 'T3': 1.0},

    'phase1_de_seeds': [42, 123, 7],
    'phase1_de_maxiter': 150,
    'phase1_de_popsize': 12,

    'de_seeds': [42, 123, 7],
    'de_maxiter': 200,
    'de_popsize': 15,

    # ── 铜电阻温度系数 (固定物理常数) ──
    'alpha_cu': 0.00393,   # 铜电阻温度系数 [1/°C]
    'T_ref': 25.0,         # 参考温度 [°C]

    'output_plot': 'lptn_result.png',
    'output_params': 'lptn_params.csv',
}

# ============================================================
# 数据加载
# ============================================================
def load_and_preprocess(csv_path, config):
    df = pd.read_csv(csv_path, encoding='utf-8-sig')
    df['timestamp'] = pd.to_datetime(df[config['time_col']], format=config['time_format'])
    df['t_sec'] = (df['timestamp'] - df['timestamp'].iloc[0]).dt.total_seconds()
    cm = config['col_map']
    data = pd.DataFrame({
        't': df['t_sec'], 'T1': df[cm['T1']].astype(float),
        'T2': df[cm['T2']].astype(float), 'T3': df[cm['T3']].astype(float),
        'Tamb': df[cm['Tamb']].astype(float), 'I': df[cm['I']].astype(float),
        'omega': df[cm['omega']].astype(float),
    })
    dt = config['resample_dt']
    data['t_bin'] = (data['t'] / dt).astype(int) * dt
    data_rs = data.groupby('t_bin').mean().reset_index()
    data_rs['t'] = data_rs['t_bin']
    data_rs.drop(columns=['t_bin'], inplace=True)
    return data_rs

def load_multiple(csv_paths, config):
    dfs = []
    for path in csv_paths:
        d = load_and_preprocess(path, config)
        if dfs:
            d['t'] = d['t'] + dfs[-1]['t'].max() + 5.0
        dfs.append(d)
    return pd.concat(dfs, ignore_index=True)


# ============================================================
# 高速仿真: 预计算 + 向量化 Euler
# ============================================================
class SimContext:
    """预计算所有不依赖参数的量, 减少每次仿真的开销"""
    def __init__(self, t, I, omega, Tamb, T0, T_meas, config, a_eff=None, C1_fixed=None):
        self.N = len(t)
        self.dt = np.diff(t)  # (N-1,)
        self.I = I.copy()
        self.I2 = I ** 2
        self.I_abs = np.abs(I)
        self.omega_abs = np.abs(omega)
        self.omega2 = omega ** 2
        self.Tamb = Tamb.copy()
        self.T0 = np.array(T0, dtype=np.float64)
        self.T_meas = T_meas.copy()  # (3, N)
        self.a_eff = a_eff      # 硬约束: a + a_lock = a_eff
        self.C1_fixed = C1_fixed  # 硬约束: C1 固定为 Phase-1 值

        # 预计算掩码
        self.lock_mask = self.omega_abs < config.get('lock_speed_thresh', 0.5)
        self.n_lock = np.sum(self.lock_mask)

        # 预计算导数相关
        self.dt_safe = np.where(self.dt > 0, self.dt, 1.0)
        self.dT_meas = np.diff(T_meas, axis=1) / self.dt_safe[np.newaxis, :]  # (3, N-1)

        # 温升掩码
        rising = self.dT_meas[0] > 0.05
        self.rising_mask = rising
        self.n_rising = np.sum(rising)

        # 权重
        ranges = [max(T_meas[i].max() - T_meas[i].min(), 1.0) for i in range(3)]
        weights = [1.0 / r**2 for r in ranges]
        scale = config.get('target_weight_scale', {'T1': 1, 'T2': 1, 'T3': 1})
        self.weights = np.array([weights[0]*scale['T1'], weights[1]*scale['T2'], weights[2]*scale['T3']])
        self.lock_weight = config.get('lock_weight', 5.0)
        self.deriv_weight = config.get('deriv_weight', 0.5)
        self.rising_weight = config.get('rising_weight', 2.0)

        # 有效dt掩码 (跳过异常间隔)
        self.valid = (self.dt > 0) & (self.dt <= 30)

        # 铜电阻温度系数
        self.alpha_cu = config.get('alpha_cu', 0.00393)
        self.T_ref = config.get('T_ref', 25.0)


# Phase-2 优化参数顺序 (15个, 不含 C1 和 a_lock):
# C2, C3, R12, R13, R23, R2a, R3a, a, b, c, d, e, f, d_lock, omega0
P2_PARAM_NAMES = ['C2','C3','R12','R13','R23','R2a','R3a',
                  'a','b','c','d','e','f','d_lock','omega0']

def _expand_params(p15, a_eff, C1_fixed):
    """将15参数 + a_eff + C1 展开为完整17参数"""
    # p15: C2, C3, R12, R13, R23, R2a, R3a, a, b, c, d, e, f, d_lock, omega0
    a = p15[7]  # a is at index 7 in p15
    a_lock = max(a_eff - a, 0.0)
    # 完整17: C1, C2, C3, R12, R13, R23, R2a, R3a, a, b, c, d, e, f, a_lock, d_lock, omega0
    return np.array([C1_fixed, *p15[:13], a_lock, p15[13], p15[14]])


def simulate_fast(params17, ctx):
    """Forward Euler - 接受完整17参数, 铜损含温度系数"""
    C1, C2, C3, R12, R13, R23, R2a, R3a = params17[:8]
    a, b, c, d_, e, f, a_lock, d_lock, omega0 = params17[8:]

    N = ctx.N
    T1 = np.empty(N)
    T2 = np.empty(N)
    T3 = np.empty(N)
    T1[0], T2[0], T3[0] = ctx.T0

    omega0_eff = max(omega0, 1e-6)
    s_lock = np.exp(-(ctx.omega_abs / omega0_eff) ** 2)

    # 铜损系数 (不含温度修正, 每步乘以温度因子)
    a_eff_arr = (a + a_lock * s_lock)  # (N,)
    # 铁耗 (不依赖温度, 可预计算)
    P1_iron = b * ctx.omega_abs + c * ctx.omega2
    # 驱动板损耗 (不依赖温度)
    P2_all = (d_ + d_lock * s_lock) * ctx.I2 + e * ctx.I_abs + f

    inv_C1 = 1.0 / C1
    inv_C2 = 1.0 / C2
    inv_C3 = 1.0 / C3
    inv_R12 = 1.0 / R12
    inv_R13 = 1.0 / R13
    inv_R23 = 1.0 / R23
    inv_R2a = 1.0 / R2a
    inv_R3a = 1.0 / R3a
    alpha_cu = ctx.alpha_cu
    T_ref = ctx.T_ref

    dt_arr = ctx.dt
    valid = ctx.valid
    Ta = ctx.Tamb
    I2 = ctx.I2

    for k in range(N - 1):
        if not valid[k]:
            T1[k+1] = T1[k]; T2[k+1] = T2[k]; T3[k+1] = T3[k]
            continue
        t1, t2, t3 = T1[k], T2[k], T3[k]
        dt = dt_arr[k]
        # 铜损: a_eff · (1 + α · (T1 - T_ref)) · I²
        P1_cu = a_eff_arr[k] * (1.0 + alpha_cu * (t1 - T_ref)) * I2[k]
        P1 = P1_cu + P1_iron[k]
        dT1 = (P1 - (t1-t2)*inv_R12 - (t1-t3)*inv_R13) * inv_C1
        dT2 = (P2_all[k] - (t2-t1)*inv_R12 - (t2-t3)*inv_R23 - (t2-Ta[k])*inv_R2a) * inv_C2
        dT3 = ((t1-t3)*inv_R13 + (t2-t3)*inv_R23 - (t3-Ta[k])*inv_R3a) * inv_C3
        T1[k+1] = t1 + dT1 * dt
        T2[k+1] = t2 + dT2 * dt
        T3[k+1] = t3 + dT3 * dt

    return np.array([T1, T2, T3])


def cost_fast(p15, ctx):
    """高速代价函数 — 15参数 + C1/a_eff硬约束"""
    try:
        params17 = _expand_params(p15, ctx.a_eff, ctx.C1_fixed)
        T_sim = simulate_fast(params17, ctx)
        if np.any(np.isnan(T_sim)) or np.any(np.abs(T_sim) > 500):
            return 1e6

        w = ctx.weights
        err = 0.0

        # 基本 MSE
        for i in range(3):
            err += w[i] * np.mean((T_sim[i] - ctx.T_meas[i]) ** 2)

        # 堵转加权 (全节点)
        if ctx.n_lock > 0:
            lm = ctx.lock_mask
            for i in range(3):
                err += ctx.lock_weight * w[i] * np.mean((T_sim[i, lm] - ctx.T_meas[i, lm]) ** 2)

        # 导数误差
        if ctx.deriv_weight > 0:
            dT_sim = np.diff(T_sim, axis=1) / ctx.dt_safe[np.newaxis, :]
            for i in range(3):
                err += ctx.deriv_weight * w[i] * np.mean((dT_sim[i] - ctx.dT_meas[i]) ** 2)

        # 温升加权
        if ctx.rising_weight > 0 and ctx.n_rising > 0:
            rm = ctx.rising_mask
            err += ctx.rising_weight * w[0] * np.mean((T_sim[0, 1:][rm] - ctx.T_meas[0, 1:][rm]) ** 2)

        return err
    except Exception:
        return 1e6


# ============================================================
# Phase-1: 堵转段单节点辨识 (分段独立仿真, 每段实测初值)
# ============================================================
def extract_lock_segments(data, config, min_pts=10):
    """按时间连续性把堵转数据切成独立段, 避免跨运行段拼接积分"""
    max_gap = 3.0 * config.get('resample_dt', 1)
    omega = data['omega'].values
    lock_mask = np.abs(omega) < config.get('lock_speed_thresh', 0.5)
    idx = np.where(lock_mask)[0]
    t_all = data['t'].values
    segments, start = [], 0
    for k in range(1, len(idx) + 1):
        end_of_run = (k == len(idx) or idx[k] != idx[k-1] + 1
                      or t_all[idx[k]] - t_all[idx[k-1]] > max_gap)
        if not end_of_run:
            continue
        seg_idx = idx[start:k]
        start = k
        if len(seg_idx) < min_pts:
            continue
        d = data.iloc[seg_idx]
        t = d['t'].values
        dt_arr = np.diff(t)
        dt_safe = np.where(dt_arr > 0, dt_arr, 1.0)
        T1 = d['T1'].values
        segments.append({
            't': t, 'I2': d['I'].values ** 2,
            'T1': T1, 'T3': d['T3'].values,
            'dt': dt_arr, 'dt_safe': dt_safe,
            'dT1_meas': np.diff(T1) / dt_safe,
        })
    return segments


def p1_sim_segment(seg, C1, a_eff, R13, alpha_cu, T_ref):
    """单段前向 Euler, 初值取该段实测 T1[0]"""
    I2 = seg['I2']; T3 = seg['T3']; dt_arr = seg['dt']
    N = len(seg['t'])
    T1s = np.empty(N); T1s[0] = seg['T1'][0]
    inv_C1 = 1.0 / C1; inv_R13 = 1.0 / R13
    for k in range(N - 1):
        # 铜损含温度系数: a_eff · (1 + α · (T1 - T_ref)) · I²
        P = a_eff * (1.0 + alpha_cu * (T1s[k] - T_ref)) * I2[k]
        dT = (P - (T1s[k] - T3[k]) * inv_R13) * inv_C1
        T1s[k+1] = T1s[k] + dT * dt_arr[k]
    return T1s


def p1_cost_segments(pvec, segments, alpha_cu, T_ref):
    """所有段残差合并 (点数加权): MSE + 0.5·导数MSE"""
    C1, a_eff, R13 = pvec
    sq = dsq = 0.0
    n = dn = 0
    try:
        for seg in segments:
            T1s = p1_sim_segment(seg, C1, a_eff, R13, alpha_cu, T_ref)
            if np.any(np.isnan(T1s)) or np.any(np.abs(T1s) > 500):
                return 1e6
            sq += np.sum((T1s - seg['T1']) ** 2); n += len(T1s)
            dT1s = np.diff(T1s) / seg['dt_safe']
            dsq += np.sum((dT1s - seg['dT1_meas']) ** 2); dn += len(dT1s)
        return sq / n + 0.5 * dsq / dn
    except Exception:
        return 1e6


def fit_phase1_segments(segments, config, seeds=None, verbose=True):
    """DE + L-BFGS-B 辨识 (C1, a_eff, R13), 支持任意段子集"""
    alpha_cu = config.get('alpha_cu', 0.00393)
    T_ref = config.get('T_ref', 25.0)
    bounds_p1 = [(5, 500), (0.1, 20), (0.1, 50)]
    if seeds is None:
        seeds = config.get('phase1_de_seeds', [42, 123, 7])
    best_c, best_x = 1e10, None
    for seed in seeds:
        if verbose: print(f"    DE-P1 (seed={seed})...", end=' ', flush=True)
        t0s = time.time()
        r = differential_evolution(p1_cost_segments, bounds_p1,
            args=(segments, alpha_cu, T_ref), seed=seed,
            maxiter=config.get('phase1_de_maxiter', 150),
            popsize=config.get('phase1_de_popsize', 12),
            tol=1e-9, polish=False, workers=1)
        r2 = minimize(p1_cost_segments, r.x, args=(segments, alpha_cu, T_ref),
            method='L-BFGS-B', bounds=bounds_p1,
            options={'maxiter': 300, 'ftol': 1e-14})
        if verbose: print(f"cost={r2.fun:.4f} ({time.time()-t0s:.0f}s)")
        if r2.fun < best_c:
            best_c, best_x = r2.fun, r2.x
    return best_x, best_c


def run_phase1(data, config):
    segments = extract_lock_segments(data, config)
    n_pts = sum(len(s['t']) for s in segments)
    if n_pts < 10:
        print("  Phase-1: 堵转数据不足, 跳过")
        return None

    print(f"  Phase-1 data: {len(segments)} 段共 {n_pts} pts (分段独立仿真)")
    for i, seg in enumerate(segments):
        print(f"    [Seg {i+1}] {seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s, "
              f"I_rms={np.sqrt(np.mean(seg['I2'])):.2f}A, "
              f"T1: {seg['T1'][0]:.1f}→{seg['T1'][-1]:.1f}°C (range={seg['T1'].max()-seg['T1'].min():.1f}°C)")

    best_x, _ = fit_phase1_segments(segments, config)
    C1, a_eff, R13 = best_x

    alpha_cu = config.get('alpha_cu', 0.00393)
    T_ref = config.get('T_ref', 25.0)
    sq, n, maxe = 0.0, 0, 0.0
    for seg in segments:
        T1s = p1_sim_segment(seg, C1, a_eff, R13, alpha_cu, T_ref)
        sq += np.sum((T1s - seg['T1'])**2); n += len(T1s)
        maxe = max(maxe, np.max(np.abs(T1s - seg['T1'])))
    rmse = np.sqrt(sq / n)
    print(f"\n  Phase-1 results:")
    print(f"    C1={C1:.2f} J/K, a_eff={a_eff:.4f} W/A²(@{T_ref}°C), R13={R13:.4f} K/W")
    print(f"    α_cu={alpha_cu:.5f} /°C (P_cu∝(1+α·ΔT))")
    print(f"    Fit: RMSE={rmse:.2f}°C, MaxErr={maxe:.2f}°C")
    return {'C1': C1, 'a_eff': a_eff, 'R13': R13}


# ============================================================
# Phase-2: 全量辨识 (16参数, a_lock = a_eff - a 硬约束)
# ============================================================
FULL_PARAM_NAMES = ['C1','C2','C3','R12','R13','R23','R2a','R3a',
                    'a','b','c','d','e','f','a_lock','d_lock','omega0']

def identify_parameters(data, config, p1=None):
    t = data['t'].values; I = data['I'].values
    omega = data['omega'].values; Tamb = data['Tamb'].values
    T_meas = np.array([data['T1'].values, data['T2'].values, data['T3'].values])
    T0 = [T_meas[0,0], T_meas[1,0], T_meas[2,0]]

    a_eff = p1['a_eff'] if p1 else 1.0
    C1_fixed = p1['C1'] if p1 else 50.0
    ctx = SimContext(t, I, omega, Tamb, T0, T_meas, config, a_eff=a_eff, C1_fixed=C1_fixed)

    bd = dict(config['bounds'])
    if p1:
        R13_p1 = p1['R13']
        bd['R13'] = (max(R13_p1*0.5, 0.1), R13_p1*1.5)
        bd['a'] = (0.05, a_eff)
        print(f"\n  Phase-1 硬约束:")
        print(f"    C1 = {C1_fixed:.2f} J/K (FIXED)")
        print(f"    a + a_lock = {a_eff:.4f} W/A² (FIXED)")
        print(f"    R13:[{bd['R13'][0]:.3f},{bd['R13'][1]:.3f}] (Phase-1参考值: {R13_p1:.3f}, 自由辨识)")
        print(f"    a:[{bd['a'][0]:.3f},{bd['a'][1]:.3f}] → a_lock=[{a_eff-bd['a'][1]:.3f},{a_eff-bd['a'][0]:.3f}]")

    # 15参数边界 (无 C1, 无 a_lock)
    bounds_list = [bd[n] for n in P2_PARAM_NAMES]

    best_cost, best_x15 = 1e10, None
    for seed in config['de_seeds']:
        print(f"\n  DE (seed={seed})...", end=' ', flush=True)
        t0s = time.time()
        r = differential_evolution(cost_fast, bounds_list,
            args=(ctx,), seed=seed,
            maxiter=config['de_maxiter'], popsize=config['de_popsize'],
            tol=1e-8, polish=False, workers=1,
            mutation=(0.5,1.5), recombination=0.9)
        print(f"cost={r.fun:.6f} ({time.time()-t0s:.0f}s)", end=' → ', flush=True)

        r2 = minimize(cost_fast, r.x, args=(ctx,),
            method='L-BFGS-B', bounds=bounds_list,
            options={'maxiter': 500, 'ftol': 1e-12})
        print(f"refined={r2.fun:.6f}")

        if r2.fun < best_cost:
            best_cost, best_x15 = r2.fun, r2.x

    print(f"\n  Best cost: {best_cost:.6f}")

    # 展开为17参数
    best_x17 = _expand_params(best_x15, a_eff, C1_fixed)

    # 碰界 (对15参数检查)
    print("  碰界检查:")
    hit = False
    for i,(name,(lo,hi)) in enumerate(zip(P2_PARAM_NAMES, bounds_list)):
        v = best_x15[i]; m = 1e-4*(hi-lo)
        if abs(v-lo)<m:
            print(f"    ⚠ {name}={v:.4f} 碰下界({lo})")
            hit = True
        elif abs(v-hi)<m:
            print(f"    ⚠ {name}={v:.4f} 碰上界({hi})")
            hit = True
    if not hit: print("    ✓ 全部在界内")

    a_val = best_x17[8]; a_lock_val = best_x17[14]
    print(f"  硬约束验证: C1={C1_fixed:.2f}(FIXED)  a={a_val:.4f}+a_lock={a_lock_val:.4f}={a_val+a_lock_val:.4f}(≡{a_eff:.4f})")

    return best_x17, T0, t, I, omega, Tamb, T_meas, ctx


# ============================================================
# 输出
# ============================================================
PARAM_UNITS = ['J/K','J/K','J/K','K/W','K/W','K/W','K/W','K/W',
               'W/A²','W·s/rad','W·s²/rad²','W/A²','W/A','W','W/A²','W/A²','rad/s']
PARAM_DESC = [
    'Winding thermal capacitance','Driver thermal capacitance','Shell thermal capacitance',
    'Winding-Driver resistance','Winding-Shell resistance','Driver-Shell resistance',
    'Driver-Ambient resistance','Shell-Ambient resistance',
    'Copper loss (I², running)','Iron loss (|ω|)','Iron loss (ω²)',
    'Driver conduction (I²)','Driver switching (|I|)','Driver standby',
    'Locked-rotor extra motor (FIXED)','Locked-rotor extra driver','Lock decay speed',
]
NODE_CN = ['节点1(绕组)', '节点2(驱动板)', '节点3(外壳)']


def print_results(params, config=None):
    if config is None: config = CONFIG
    C1,C2,C3,R12,R13,R23,R2a,R3a,a,b,c,d,e,f,a_lock,d_lock,omega0 = params
    print("\n" + "="*60 + "\nIDENTIFIED PARAMETERS\n" + "="*60)
    for n,v,u in zip(FULL_PARAM_NAMES, params, PARAM_UNITS):
        print(f"  {n:8s} = {v:12.6f}  {u}")
    print(f"""
  s(ω) = exp(-(|ω|/{omega0:.4f})²)
  P1 = ({a:.4f}+{a_lock:.4f}·s)·(1+{config.get('alpha_cu',0.00393):.5f}·(T1-{config.get('T_ref',25)}))·I² + {b:.4f}·|ω| + {c:.6f}·ω²
  P2 = ({d:.4f}+{d_lock:.4f}·s)·I² + {e:.4f}·|I| + {f:.4f}
  τ1≈{C1*(1/(1/R12+1/R13)):.1f}s  τ2≈{C2*(1/(1/R12+1/R23+1/R2a)):.1f}s  τ3≈{C3*(1/(1/R13+1/R23+1/R3a)):.1f}s""")


def evaluate_and_plot(params, T0, t, I, omega, Tamb, T_meas, config):
    T_sim = simulate_fast(params, SimContext(t, I, omega, Tamb, T0, T_meas, config))

    print("\nFITTING ACCURACY\n" + "-"*50)
    for i in range(3):
        rmse = np.sqrt(np.mean((T_sim[i]-T_meas[i])**2))
        mae = np.mean(np.abs(T_sim[i]-T_meas[i]))
        maxe = np.max(np.abs(T_sim[i]-T_meas[i]))
        print(f"  {NODE_CN[i]}: RMSE={rmse:.2f}°C, MAE={mae:.2f}°C, MaxErr={maxe:.2f}°C")

    lock_mask = np.abs(omega) < config.get('lock_speed_thresh', 0.5)
    if np.any(lock_mask):
        print("  ── 堵转段 ──")
        for i in range(3):
            r = np.sqrt(np.mean((T_sim[i,lock_mask]-T_meas[i,lock_mask])**2))
            m = np.max(np.abs(T_sim[i,lock_mask]-T_meas[i,lock_mask]))
            print(f"    {NODE_CN[i]}: RMSE={r:.2f}°C, MaxErr={m:.2f}°C")

    a,b,c = params[8:11]; d_,e,f = params[11:14]
    a_lock,d_lock,omega0 = params[14:17]
    alpha_cu = config.get('alpha_cu', 0.00393)
    T_ref = config.get('T_ref', 25.0)
    s = np.exp(-(np.abs(omega)/max(omega0,1e-6))**2)
    a_eff_arr = a + a_lock * s
    # 铜损用T_sim[0] (T1) 计算温度修正
    P1 = a_eff_arr * (1.0 + alpha_cu * (T_sim[0] - T_ref)) * I**2 + b*np.abs(omega) + c*omega**2
    P2 = (d_+d_lock*s)*I**2 + e*np.abs(I) + f

    fig, axes = plt.subplots(5, 1, figsize=(14,20), sharex=True,
                             gridspec_kw={'height_ratios':[3,3,3,2,2]})
    cm = ['#e74c3c','#3498db','#2ecc71']
    cs = ['#c0392b','#2980b9','#27ae60']
    titles = ['Node 1: Winding','Node 2: Driver (MOS)','Node 3: Shell']
    lm = ['T1 meas','T2 meas','T3 meas']

    for i in range(3):
        ax = axes[i]
        rmse = np.sqrt(np.mean((T_sim[i]-T_meas[i])**2))
        ax.plot(t, T_meas[i], color=cm[i], alpha=.6, lw=1, label=lm[i])
        ax.plot(t, T_sim[i], color=cs[i], lw=1.8, ls='--', label=f'Sim (RMSE={rmse:.2f}°C)')
        if np.any(lock_mask):
            yl = (min(T_meas[i].min(),T_sim[i].min())-2, max(T_meas[i].max(),T_sim[i].max())+2)
            ax.fill_between(t, yl[0], yl[1], where=lock_mask, alpha=.06, color='orange', label='Locked')
            ax.set_ylim(yl)
        if i==2: ax.plot(t, Tamb, color='gray', alpha=.4, lw=.8, label='Tamb')
        ax.set_ylabel('Temp (°C)'); ax.legend(fontsize=9, loc='upper left')
        ax.set_title(titles[i], fontweight='bold'); ax.grid(True, alpha=.3)

    ax = axes[3]
    for i in range(3):
        ax.plot(t, T_sim[i]-T_meas[i], color=cm[i], alpha=.7, lw=.8, label=f'{NODE_CN[i]} err')
    ax.axhline(0, color='k', alpha=.3, lw=.5)
    ax.set_ylabel('Error (°C)'); ax.legend(fontsize=9, loc='upper left')
    ax.set_title('Error', fontweight='bold'); ax.grid(True, alpha=.3)

    ax = axes[4]; ax2 = ax.twinx()
    l1=ax.plot(t,I,color='#9b59b6',alpha=.3,lw=.5,label='I(A)')
    l2=ax2.plot(t,omega,color='#e67e22',alpha=.3,lw=.5,label='ω(rad/s)')
    l3=ax.plot(t,P1,color='#e74c3c',alpha=.6,lw=1,label='P1(W)')
    l4=ax.plot(t,P2,color='#3498db',alpha=.6,lw=1,label='P2(W)')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('I(A)/Loss(W)'); ax2.set_ylabel('ω(rad/s)')
    lines = l1+l2+l3+l4
    ax.legend(lines,[l.get_label() for l in lines], fontsize=9, loc='upper left')
    ax.set_title('Inputs & Losses', fontweight='bold'); ax.grid(True, alpha=.3)

    fig.suptitle('3-Node LPTN v3: Multi-Phase Identification', fontsize=14, fontweight='bold', y=.995)
    plt.tight_layout(rect=[0,0,1,.98])
    plt.savefig(config['output_plot'], dpi=150, bbox_inches='tight'); plt.close()
    print(f"\nPlot: {config['output_plot']}")

    pd.DataFrame({'Parameter':FULL_PARAM_NAMES,'Value':params,'Unit':PARAM_UNITS,'Description':PARAM_DESC}
                 ).to_csv(config['output_params'], index=False, encoding='utf-8-sig')
    print(f"Params: {config['output_params']}")


# ============================================================
# 主入口
# ============================================================
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)

    csv_files = sys.argv[1:]
    print("="*60)
    print("3-Node LPTN v3: Multi-Phase Identification")
    print("="*60)

    print(f"\n[1/4] Loading {len(csv_files)} file(s)...")
    data = load_multiple(csv_files, CONFIG)
    print(f"  {len(data)} pts, {data['t'].max()-data['t'].min():.0f}s")
    print(f"  T1:{data['T1'].min():.1f}~{data['T1'].max():.1f} T2:{data['T2'].min():.1f}~{data['T2'].max():.1f} T3:{data['T3'].min():.1f}~{data['T3'].max():.1f}")
    print(f"  I:{data['I'].min():.1f}~{data['I'].max():.1f}A  ω:{data['omega'].min():.1f}~{data['omega'].max():.1f}rad/s")
    lf = np.mean(np.abs(data['omega'].values) < CONFIG['lock_speed_thresh'])
    print(f"  Locked-rotor: {lf*100:.1f}%")

    t_total = time.time()

    print(f"\n[2/4] Phase-1: locked-rotor → C1, a_eff, R13")
    p1 = run_phase1(data, CONFIG)

    print(f"\n[3/4] Phase-2: full data → 15 params (C1 FIXED, a_lock = a_eff - a)")
    params, T0, t, I, omega, Tamb, T_meas, ctx = identify_parameters(data, CONFIG, p1)
    print_results(params, CONFIG)

    print(f"\n[4/4] Validation & plot...")
    evaluate_and_plot(params, T0, t, I, omega, Tamb, T_meas, CONFIG)

    print(f"\nTotal time: {time.time()-t_total:.0f}s")
    print("Done!")
