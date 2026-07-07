"""
三节点(绕组–铁心–壳) + 热历史, 单一恒定系数 —— 检验慢核能否解释 2×
================================================================
拓扑: 铜损 P 全入绕组 Tw(=NTC实测); 隐藏铁心 Tc(慢热质量);
      壳 T3 实测作边界。 Tw—R_ws—Shell(直接快路),
      Tw—R_wc—Tc—R_cs—Shell(经核慢路)。
整段连续仿真(含 E1→冷却→E2→冷却), 核状态自动携带热历史。
  dTw/dt = [P - (Tw-Tc)/R_wc - (Tw-T3)/R_ws] / C_w
  dTc/dt = [(Tw-Tc)/R_wc - (Tc-T3)/R_cs]      / C_c
  P = a·(1+α(Tw-Tref))·I²     ← 单一恒定 a
判据: 恒定 a 能否同时到 E1 冷启峰(99°C) 和 E2 温热(90°C)?
  能  → 2× 由热历史(慢核)解释, 非电气
  不能→ 峰仍够不到, 排除热历史, 支持真实额外电气损耗
用法: python phase1_3node_history.py <CSV>
"""
import sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lptn_3node_final import CONFIG, load_multiple, extract_lock_segments

MAX_SUB = 60
K_TREND = 8.0
COLD = (25, 72)
STEADY = (330, 375)
OUT = 'phase1_3node_history.png'
# C_w, C_c, R_wc, R_ws, R_cs, a
BND = [(2,200),(20,3000),(0.02,5),(0.05,20),(0.02,10),(0.1,12)]


def weights(seg):
    w = 1.0 + K_TREND*np.maximum(seg['dT1_meas']/seg['dt_safe'], 0.0)
    return np.append(w, w[-1])


def sim(seg, p, alpha, Tref):
    C_w,C_c,R_wc,R_ws,R_cs,a = p
    I2=seg['I2']; T3=seg['T3']; dt_arr=seg['dt']; N=len(seg['t'])
    Tw=np.empty(N); Tc=np.empty(N); Tw[0]=Tc[0]=seg['T1'][0]
    iCw=1/C_w; iCc=1/C_c; iwc=1/R_wc; iws=1/R_ws; ics=1/R_cs
    tau_min=min(C_w*R_wc, C_w*R_ws, C_c*R_cs)
    for k in range(N-1):
        dt=dt_arr[k]; m=int(dt/(0.4*tau_min))+1
        if m>MAX_SUB: return None,None
        h=dt/m; tw,tc=Tw[k],Tc[k]; aI2=a*I2[k]; t3=T3[k]
        for _ in range(m):
            P=aI2*(1+alpha*(tw-Tref))
            dtw=(P-(tw-tc)*iwc-(tw-t3)*iws)*iCw
            dtc=((tw-tc)*iwc-(tc-t3)*ics)*iCc
            tw+=dtw*h; tc+=dtc*h
        Tw[k+1]=tw; Tc[k+1]=tc
    return Tw,Tc


def cost(p, seg, w, alpha, Tref):
    try:
        Tw,_=sim(seg,p,alpha,Tref)
        if Tw is None or np.any(np.isnan(Tw)) or np.any(np.abs(Tw)>500): return 1e6
        return np.sum(w*(Tw-seg['T1'])**2)/np.sum(w)
    except Exception: return 1e6


def fit(seg, w, alpha, Tref, seeds=(42,123,7), mi=250):
    bc,bx=1e10,None
    for s in seeds:
        r=differential_evolution(cost,BND,args=(seg,w,alpha,Tref),seed=s,
            maxiter=mi,popsize=16,tol=1e-9,polish=False,workers=1,
            mutation=(0.5,1.5),recombination=0.9)
        r2=minimize(cost,r.x,args=(seg,w,alpha,Tref),method='L-BFGS-B',
            bounds=BND,options={'maxiter':500,'ftol':1e-14})
        if r2.fun<bc: bc,bx=r2.fun,r2.x
    return bx


def pk(sim_, seg, win):
    m=(seg['t']>=win[0])&(seg['t']<win[1]); return sim_[m].max() if m.any() else np.nan
def rm(sim_, seg, win=None):
    if win is None: return np.sqrt(np.mean((sim_-seg['T1'])**2))
    m=(seg['t']>=win[0])&(seg['t']<win[1]); return np.sqrt(np.mean((sim_[m]-seg['T1'][m])**2))


def main():
    if len(sys.argv)<2: print(__doc__); sys.exit(1)
    alpha=CONFIG.get('alpha_cu',0.00393); Tref=CONFIG.get('T_ref',25.0)
    data=load_multiple([sys.argv[1]],CONFIG)
    seg=max(extract_lock_segments(data,CONFIG),key=lambda s:len(s['t']))
    w=weights(seg)
    T3c=seg['T3'][(seg['t']>=COLD[0])&(seg['t']<COLD[1])].mean()
    T3s=seg['T3'][(seg['t']>=STEADY[0])&(seg['t']<STEADY[1])].mean()
    print(f"堵转段 t={seg['t'][0]:.0f}~{seg['t'][-1]:.0f}s")
    print(f"实测冷启峰={pk(seg['T1'],seg,COLD):.0f}°C (壳T3≈{T3c:.0f}°C)  温热末T1≈{seg['T1'][(seg['t']>=STEADY[0])&(seg['t']<STEADY[1])].max():.0f}°C (壳T3≈{T3s:.0f}°C)")
    print(f"→ 恒定a下 E1渐近线≈T3_cold+P·R, 比E2低约 {T3s-T3c:.0f}°C (壳温差)\n")

    print("[三节点+热历史, 恒定a, 瞬态加权]"); t0=time.time()
    x=fit(seg,w,alpha,Tref)
    Tw,Tc=sim(seg,x,alpha,Tref)
    C_w,C_c,R_wc,R_ws,R_cs,a=x
    print(f"  C_w={C_w:.1f} C_c={C_c:.0f} R_wc={R_wc:.3f} R_ws={R_ws:.3f} R_cs={R_cs:.3f} a={a:.3f}  ({time.time()-t0:.0f}s)")
    print(f"  τ_core=C_c·R_cs={C_c*R_cs:.0f}s  τ_wind≈C_w·(R_wc∥R_ws)={C_w/(1/R_wc+1/R_ws):.1f}s")
    print(f"  冷启峰模拟={pk(Tw,seg,COLD):.0f}°C (实测99)  RMSE冷启={rm(Tw,seg,COLD):.2f} 准稳态={rm(Tw,seg,STEADY):.2f} 全段={rm(Tw,seg):.2f}")
    for v,(lo,hi),n in zip(x,BND,['C_w','C_c','R_wc','R_ws','R_cs','a']):
        if abs(v-lo)<1e-3*(hi-lo): print(f"    ⚠ {n} 碰下界 {lo}")
        if abs(v-hi)<1e-3*(hi-lo): print(f"    ⚠ {n} 碰上界 {hi}")
    print(f"\n  判据: 冷启峰{'达到' if pk(Tw,seg,COLD)>=97 else '够不到'}99°C → "
          f"{'热历史可解释2×' if pk(Tw,seg,COLD)>=97 else '排除热历史, 支持真实额外电气损耗'}")

    t=seg['t']-seg['t'][0]
    fig,(a1,a2)=plt.subplots(1,2,figsize=(15,5.2))
    for ax in (a1,a2):
        ax.plot(t,seg['T1'],color='#6b7280',lw=1.5,label='Tw measured (NTC)')
        ax.plot(t,Tw,color='#d97706',lw=2.0,ls='-.',label=f'3-node const-a (E1 peak={pk(Tw,seg,COLD):.0f}°C)')
        ax.plot(t,Tc,color='#0e7a5f',lw=1.2,ls=':',label='core Tc (hidden, thermal history)')
        ax.plot(t,seg['T3'],color='#3477c9',lw=1.0,alpha=.5,label='shell T3 (boundary)')
        ax.grid(True,alpha=.25); ax.set_xlabel('t within segment (s)'); ax.set_ylabel('T (°C)')
    a1.legend(fontsize=8.5,loc='best'); a1.set_title('Full locked segment',fontsize=11)
    a2.set_xlim(0,COLD[1]-seg['t'][0]+3); a2.set_title('Cold-start (E1) zoom',fontsize=11)
    fig.suptitle('3-node winding–core–shell + thermal history, constant coefficient',fontweight='bold')
    plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig(OUT,dpi=150,bbox_inches='tight')
    print(f'\nPlot: {OUT}')


if __name__=='__main__':
    main()
