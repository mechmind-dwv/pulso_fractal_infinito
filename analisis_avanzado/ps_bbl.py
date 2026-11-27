import sys
import numpy as np
# E26: P(s) de los ceros frente a GUE (sine) y al nucleo CUE con N_eff continuo (Bogomolny et al. 2006).
# Uso: python3 ps_bbl.py ceros_100k.txt [t_min] [t_max]
f = sys.argv[1]
tmin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 1e99
z = np.loadtxt(f); z = z[(z >= tmin) & (z <= tmax)]
def Nsm(t): return t/(2*np.pi)*np.log(t/(2*np.pi)) - t/(2*np.pi) + 7/8 + 1/(48*np.pi*t)
u = Nsm(z); s = np.diff(u)
tm = z.mean(); Lam = 1.57314
Neff = np.log(tm/(2*np.pi))/np.sqrt(12*Lam)
print(f"{len(z)} ceros, t medio {tm:.4g}, N_eff = {Neff:.3f}, <s> = {s.mean():.4f}")
s = s / s.mean()
gx, gw = np.polynomial.legendre.leggauss(40)
def E0(x, N):
    # gap probability det(1-K) en (0,x), nucleo sine (N=inf) o CUE continuo
    t = 0.5*x*(gx+1); w = 0.5*x*gw
    d = t[:, None]-t[None, :]
    if N is None: K = np.sinc(d)
    else:
        with np.errstate(all="ignore"):
            K = np.sin(np.pi*d)/(N*np.sin(np.pi*d/N))
        K[np.abs(d) < 1e-14] = 1.0
    A = np.eye(len(t)) - np.sqrt(w)[:, None]*K*np.sqrt(w)[None, :]
    return np.linalg.det(A)
def Pcurve(N, grid):
    h = 0.02; x = np.arange(h, 4.0+2*h, h)
    e = np.array([E0(v, N) for v in x]); e = np.concatenate([[1.0], e])  # E(0)=1
    xx = np.concatenate([[0.0], x])
    p = (e[2:]-2*e[1:-1]+e[:-2])/h**2
    return np.interp(grid, xx[1:-1], p)
edges = np.arange(0, 3.01, 0.1); mid = 0.5*(edges[1:]+edges[:-1])
cnt, _ = np.histogram(s, bins=edges); n = len(s)
fine = np.linspace(0, 3, 301)
res = {}
for nombre, N in (("GUE (sine)", None), (f"CUE N_eff={Neff:.2f}", Neff), ("CUE N=3", 3.0), ("CUE N=4", 4.0)):
    p = Pcurve(N, fine); cdf = np.cumsum(p)*0.01
    pb = np.interp(edges, fine, cdf); pb = np.diff(pb)
    ex = pb*n
    ok = ex > 5
    chi = np.sum((cnt[ok]-ex[ok])**2/ex[ok])
    res[nombre] = (chi, ok.sum(), pb.sum())
    print(f"{nombre:22s} chi2/bin = {chi/ok.sum():.2f}  (bins {ok.sum()}, masa {pb.sum():.3f})")
