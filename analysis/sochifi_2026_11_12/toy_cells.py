# toy_cells.py — SOCHIFI 2026-11-12: comprobación con un toy (ruido independiente por tira, dos tiras ruidosas) de la resta de accidentales por
# celda de impact_map.py: la resta de mezclados celda a celda está sesgada; el estimador condicionado reproduce la verdad. Uso: toy_cells.py
# toy: sesgo por celda de on − mezclados (+ rXrY·den con la forma verdadera) frente a la verdad q·N·s_c
import numpy as np
rng=np.random.default_rng(1); N=400000; NC=31
def cells(H):
    n=H.sum(1); f=np.argmax(H,1); c=np.full(len(H),-1); c[n==1]=2*f[n==1]
    adj=(n==2)&(f<15)&H[np.arange(len(H)),np.minimum(f+1,15)]; c[adj]=2*f[adj]+1; return c
hist=lambda c: np.bincount(np.where(c<0,NC,c),minlength=NC+1).astype(float)
for pn in (0.019,0.005):   # ruido por tira: aX≈0.26 (5 p.e.) y ≈0.08 (10 p.e.)
    noisy=np.full(16,pn); noisy[[2,9]]*=4   # dos tiras ruidosas
    q, rx_only, ry_only = 0.028, 0.034, 0.047
    u=rng.random(N); tx=u<q+rx_only; ty=(u<q)|((u>=q+rx_only)&(u<q+rx_only+ry_only))
    def true_pat(m):
        H=np.zeros((N,16),bool); pos=rng.uniform(-16.5,16.5,N); ctr=-15+2*np.arange(16)
        H[:]=(np.abs(pos[:,None]-ctr[None,:])<=1.5)&m[:,None]; return H
    TX,TY=true_pat(tx),true_pat(ty)
    def noise(): return rng.random((N,16))<noisy[None,:]
    HX=TX|noise(); HY=TY|noise(); HXc=noise(); HYc=noise()
    Xs,Ys,Xc,Yc=HX.any(1),HY.any(1),HXc.any(1),HYc.any(1)
    aX,aY=Xc.mean(),Yc.mean(); rX=(Xs.mean()-aX)/(1-aX); rY=(Ys.mean()-aY)/(1-aY); den=(1-aX)*(1-aY)
    cx,cxc=cells(HX),cells(HXc); on=Xs&Ys; s=N/(N-1)
    m1=np.zeros(N,bool); m1[:-1]=Xs[:-1]&Ys[1:]; m2=np.zeros(N,bool); m2[:-1]=Ys[:-1]&Xs[1:]
    mix=0.5*(hist(cx[:-1][m1[:-1]])+hist(cx[1:][m2[:-1]]))*s
    shape=np.clip(hist(cx[Xs])-hist(cxc[Xc]),0,None); shape/=shape.sum()
    sub=hist(cx[on])-mix+rX*rY*N*den*shape
    truth=hist(cells(TX)[tx&ty])
    nc=hist(cxc)/N   # patrón de ruido solo (ventana fuera de tiempo)
    pred_bias=-N*nc*(q+ (rX-q)*aY)   # sesgo esperado por migración del ruido
    print(f"ruido/tira {pn}: aX {aX:.3f}; total sub {sub.sum():.0f}, truth {truth.sum():.0f}, q·N·den {q*N*den:.0f}")
    for c in (0,2,3,4,5,8,18,19,30):
        print(f"  celda {c:2d}: verdad {truth[c]:6.0f}  sub {sub[c]:7.1f}  sub/(1-aX) {sub[c]/(1-aX):7.1f}  (sub - pred_bias)/(1-aX) {(sub[c]-pred_bias[c])/(1-aX):7.1f}")
print("--- estimador condicionado: raw_c = (q/P_XYs)(on_c − N n_c P(Ys&noX)), P_XYs = q + (rX − q) aY, q del mezclado")
for pn in (0.019,0.005,0.002):
    noisy=np.full(16,pn); noisy[[2,9]]*=4
    q, rx_only, ry_only = 0.028, 0.034, 0.047
    u=rng.random(N); tx=u<q+rx_only; ty=(u<q)|((u>=q+rx_only)&(u<q+rx_only+ry_only))
    pos=rng.uniform(-16.5,16.5,N); ctr=-15+2*np.arange(16)
    TX=(np.abs(pos[:,None]-ctr[None,:])<=1.5)&tx[:,None]; pos2=rng.uniform(-16.5,16.5,N); TY=(np.abs(pos2[:,None]-ctr[None,:])<=1.5)&ty[:,None]
    nz=lambda: rng.random((N,16))<noisy[None,:]
    HX=TX|nz(); HY=TY|nz(); HXc=nz(); HYc=nz()
    Xs,Ys,Xc,Yc=HX.any(1),HY.any(1),HXc.any(1),HYc.any(1)
    aX,aY=Xc.mean(),Yc.mean(); rX=(Xs.mean()-aX)/(1-aX); rY=(Ys.mean()-aY)/(1-aY)
    m1=(Xs[:-1]&Ys[1:]).sum(); m2=(Ys[:-1]&Xs[1:]).sum(); nxy=(Xs&Ys).sum(); acc=0.5*(m1+m2)*N/(N-1)
    qe=(nxy-acc)/(N*(1-aX)*(1-aY))+rX*rY
    PXYs=qe+(rX-qe)*aY; PYnoX=Ys.mean()-PXYs
    cx=cells(HX); on=hist(cx[Xs&Ys]); nc=hist(cells(HXc))/N
    raw=(qe/PXYs)*(on-N*nc*PYnoX)
    truth=hist(cx[tx&ty])
    print(f"ruido {pn}: aX {aX:.3f} q_est {qe:.4f} (q 0.028); suma raw celdas+otra {raw.sum():.0f}, verdad {truth.sum():.0f}")
    for c in (0,2,3,4,8,18,19,30,31):
        print(f"  celda {c:2d}: verdad(obs) {truth[c]:6.0f}  raw {raw[c]:7.1f}")
