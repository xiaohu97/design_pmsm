"""
两节点绕组(局部热点+本体) + 恒定损耗系数 + 冷启瞬态加权
=======================================================
检验: 单一恒定 a 下, 两时间尺度结构能否同时拟合
      冷启尖峰(E1, ~98°C) 和温热准稳态(~90°C)?
  能  → "额外损耗"主因是热结构(缺快热容), 非电气
  不能→ 冷启仍需更高功率, 存在真实额外电气损耗

模型(T3实测作边界):
  dTh/dt = [γ·P(Th)        - (Th-Tb)/R_h] / C_h    (Th = NTC 读数)
  dTb/dt = [(1-γ)P(Tb) + (Th-Tb)/R_h - (Tb-T3)/R13] / C_b
  P = a·(1+α(T-Tref))·I²   ← 单一恒定 a
判据: 冷启窗(E1)峰值与RMSE 单列, 对比单节点(同样加权)
用法: python phase1_twonode_weighted.py <CSV>
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import CONFIG, load_multiple, extract_lock_segments, p1_sim_segment

MAX_SUB = 40
K_TREND = 8.0        # 瞬态加权强度: w = 1 + K·relu(dT1/dt)
COLD = (25, 72)      # 冷启窗(E1) 用于单列判据
STEADY = (330, 375)  # 温热准稳态窗
OUT = 'phase1_twonode_weighted.png'

# 两节点参数: C_h, C_b, R_h, R13, a, gamma
TN_BOUNDS = [(1,150),(20,800),(0.02,5),(0.05,10),(0.1,12),(0.05,1.0)]
# 单节点参数: C1, R13, a
SN_BOUNDS = [(5,500),(0.05,10),(0.1,20)]


def weights(seg):
    w = 1.0 + K_TREND*np.maximum(seg['dT1_meas']/seg['dt_safe'], 0.0)  # len N-1
    return np.append(w, w[-1])


def tn_sim(seg, p, alpha, Tref):
    C_h,C_b,R_h,R13,a,g = p
    I2=seg['I2']; T3=seg['T3']; dt_arr=seg['dt']; N=len(seg['t'])
    Th=np.empty(N); Tb=np.empty(N); Th[0]=Tb[0]=seg['T1'][0]
    iCh=1/C_h; iCb=1/C_b; iRh=1/R_h; iR13=1/R13; tau=C_h*R_h
    for k in range(N-1):
        dt=dt_arr[k]; m=int(dt/(0.4*tau))+1
        if m>MAX_SUB: return None,None
        h=dt/m; th,tb=Th[k],Tb[k]; aI2=a*I2[k]
        for _ in range(m):
            Ph=g*aI2*(1+alpha*(th-Tref)); Pb=(1-g)*aI2*(1+alpha*(tb-Tref))
            q=(th-tb)*iRh
            th+=(Ph-q)*iCh*h; tb+=(Pb+q-(tb-T3[k])*iR13)*iCb*h
        Th[k+1]=th; Tb[k+1]=tb
    return Th,Tb


def tn_cost(p, seg, w, alpha, Tref):
    try:
        Th,_=tn_sim(seg,p,alpha,Tref)
        if Th is None or np.any(np.isnan(Th)) or np.any(np.abs(Th)>500): return 1e6
        return np.sum(w*(Th-seg['T1'])**2)/np.sum(w)
    except Exception: return 1e6


def sn_cost(p, seg, w, alpha, Tref):
    try:
        C1,R13,a=p
        T=p1_sim_segment(seg,C1,a,R13,alpha,Tref)
        if np.any(np.isnan(T)): return 1e6
        return np.sum(w*(T-seg['T1'])**2)/np.sum(w)
    except Exception: return 1e6


def fit(cost, bounds, seg, w, alpha, Tref, seeds=(42,123), mi=200):
    bc,bx=1e10,None
    for s in seeds:
        r=differential_evolution(cost,bounds,args=(seg,w,alpha,Tref),seed=s,
            maxiter=mi,popsize=15,tol=1e-9,polish=False,workers=1,
            mutation=(0.5,1.5),recombination=0.9)
        r2=minimize(cost,r.x,args=(seg,w,alpha,Tref),method='L-BFGS-B',
            bounds=bounds,options={'maxiter':400,'ftol':1e-14})
        if r2.fun<bc: bc,bx=r2.fun,r2.x
    return bx


def wrmse(sim, seg, lo, hi):
    t=seg['t']; m=(t>=lo)&(t<hi)
    return np.sqrt(np.mean((sim[m]-seg['T1'][m])**2)) if m.any() else np.nan


def main():
    if len(sys.argv)<2: print(__doc__); sys.exit(1)
    alpha=CONFIG.get('alpha_cu',0.00393); Tref=CONFIG.get('T_ref',25.0)
    data=load_multiple([sys.argv[1]],CONFIG)
    segs=extract_lock_segments(data,CONFIG)
    seg=max(segs,key=lambda s:len(s['t']))   # 取最长堵转段(含E1+再加热+冷却)
    w=weights(seg)
    print(f"堵转段: t={seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s, {len(seg['t'])}pts")
    print(f"实测冷启峰 T1_max(在{COLD}) = {seg['T1'][(seg['t']>=COLD[0])&(seg['t']<COLD[1])].max():.0f}°C\n")

    print("[单节点 恒定a, 瞬态加权]"); t0=time.time()
    xs=fit(sn_cost,SN_BOUNDS,seg,w,alpha,Tref)
    Ts=p1_sim_segment(seg,xs[0],xs[2],xs[1],alpha,Tref)
    print(f"  C1={xs[0]:.1f} R13={xs[1]:.3f} a={xs[2]:.3f}  ({time.time()-t0:.0f}s)")
    print(f"  冷启峰模拟={Ts[(seg['t']>=COLD[0])&(seg['t']<COLD[1])].max():.0f}°C  "
          f"RMSE冷启={wrmse(Ts,seg,*COLD):.2f} 准稳态={wrmse(Ts,seg,*STEADY):.2f} 全段={np.sqrt(np.mean((Ts-seg['T1'])**2)):.2f}\n")

    print("[两节点 恒定a, 瞬态加权]"); t0=time.time()
    xt=fit(tn_cost,TN_BOUNDS,seg,w,alpha,Tref)
    Th,Tb=tn_sim(seg,xt,alpha,Tref)
    C_h,C_b,R_h,R13,a,g=xt
    print(f"  C_h={C_h:.1f} C_b={C_b:.1f} R_h={R_h:.3f} R13={R13:.3f} a={a:.3f} γ={g:.2f}  ({time.time()-t0:.0f}s)")
    print(f"  τ_hot=C_h·R_h={C_h*R_h:.1f}s  τ_bulk=C_b·R13={C_b*R13:.1f}s")
    print(f"  冷启峰模拟={Th[(seg['t']>=COLD[0])&(seg['t']<COLD[1])].max():.0f}°C  "
          f"RMSE冷启={wrmse(Th,seg,*COLD):.2f} 准稳态={wrmse(Th,seg,*STEADY):.2f} 全段={np.sqrt(np.mean((Th-seg['T1'])**2)):.2f}")
    # 碰界
    for v,(lo,hi),n in zip(xt,TN_BOUNDS,['C_h','C_b','R_h','R13','a','γ']):
        if abs(v-lo)<1e-3*(hi-lo): print(f"    ⚠ {n} 碰下界 {lo}")
        if abs(v-hi)<1e-3*(hi-lo): print(f"    ⚠ {n} 碰上界 {hi}")

    # 图: 全段 + E1 放大
    t=seg['t']-seg['t'][0]
    fig,(a1,a2)=plt.subplots(1,2,figsize=(15,5.2))
    for ax in (a1,a2):
        ax.plot(t,seg['T1'],color='#6b7280',lw=1.5,label='T1 measured (NTC)')
        ax.plot(t,Ts,color='#3477c9',lw=1.6,ls='--',label=f'single-node (peak={Ts[(seg["t"]>=COLD[0])&(seg["t"]<COLD[1])].max():.0f}°C)')
        ax.plot(t,Th,color='#d97706',lw=2.0,ls='-.',label=f'two-node Th (peak={Th[(seg["t"]>=COLD[0])&(seg["t"]<COLD[1])].max():.0f}°C)')
        ax.plot(t,Tb,color='#d97706',lw=1.0,ls=':',alpha=.6,label='two-node Tb (bulk)')
        ax.grid(True,alpha=.25); ax.set_xlabel('t within segment (s)'); ax.set_ylabel('T (°C)')
    a1.legend(fontsize=8.5,loc='best'); a1.set_title('Full locked segment',fontsize=11)
    a2.set_xlim(0, COLD[1]-seg['t'][0]+3); a2.set_title('Cold-start transient (E1) zoom',fontsize=11)
    fig.suptitle('Two-node (hotspot) winding, constant coefficient, transient-weighted',fontweight='bold')
    plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig(OUT,dpi=150,bbox_inches='tight')
    print(f'\nPlot: {OUT}')


if __name__=='__main__':
    main()
