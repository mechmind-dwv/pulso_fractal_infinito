#!/usr/bin/env python3
"""Extrae fase desenrollada, velocidad, curvatura y compara con theta(t) de
Riemann-Siegel. Entrada: salida_pfi01_N{N}/trayectoria.csv
Uso: python3 lab_phase.py 1000
"""
import sys, csv, math
import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
base = f"salida_pfi01_N{N}"

# 1) cargar trayectoria
ts, xs, ys, Zs = [], [], [], []
with open(f"{base}/trayectoria.csv") as f:
    r = csv.DictReader(f)
    # nombres tolerantes
    fields = r.fieldnames or []
    def pick(*names):
        for n in names:
            if n in fields: return n
        return fields[0]
    kt, kx, ky, kZ = pick("t"), pick("Re(Gamma)", "Re", "X"), pick("Im(Gamma)", "Im", "Y"), pick("Z(t)", "Z")
    for row in r:
        ts.append(float(row[kt])); xs.append(float(row[kx]))
        ys.append(float(row[ky])); Zs.append(float(row[kZ]))

t = np.asarray(ts); x = np.asarray(xs); y = np.asarray(ys); Z = np.asarray(Zs)
print(f"N={N}  puntos={len(t)}  t∈[{t[0]:.2f},{t[-1]:.2f}]")

# 2) fase desenrollada
phi = np.unwrap(np.arctan2(y, x))

# 3) theta(t) Riemann-Siegel
from mpmath import mp, gamma, log, pi, arg, mpf
mp.dps = 15
def theta_rs(tv):
    return float(arg(gamma(mpf(0.25) + 0.5j*tv)) - (tv/2.0)*log(pi))
# vectorizar solo puntos muestreados (mpmath es lento)
sample_idx = np.linspace(0, len(t)-1, 200).astype(int)
thetas = np.array([theta_rs(float(t[i])) for i in sample_idx])
phi_s = phi[sample_idx]
ts_s  = t[sample_idx]

# 4) velocidad y curvatura (diferencias finitas centradas)
def deriv(arr, t):
    out = np.zeros_like(arr)
    out[1:-1] = (arr[2:] - arr[:-2]) / (t[2:] - t[:-2])
    return out
xp, yp = deriv(x, t), deriv(y, t)
xpp, ypp = deriv(xp, t), deriv(yp, t)
speed = np.sqrt(xp**2 + yp**2)
k = np.abs(xp*ypp - yp*xpp) / (speed**3 + 1e-30)

# 5) guardar
with open(f"{base}/fase_curvatura.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t", "phi", "dphi_dt", "speed", "curvature", "R"])
    dphi = np.zeros_like(phi)
    dphi[1:-1] = (phi[2:] - phi[:-2]) / (t[2:] - t[:-2])
    for i in range(len(t)):
        w.writerow([f"{t[i]:.10f}", f"{phi[i]:.10f}", f"{dphi[i]:.10f}",
                    f"{speed[i]:.10f}", f"{k[i]:.10f}",
                    f"{math.hypot(x[i],y[i]):.10f}"])

with open(f"{base}/theta_comparacion.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t", "phi", "theta_rs", "phi_minus_theta", "ratio_phi_theta"])
    for tt, pp, th in zip(ts_s, phi_s, thetas):
        w.writerow([f"{tt:.10f}", f"{pp:.10f}", f"{th:.10f}",
                    f"{pp-th:.10f}", f"{(pp-th)/tt if tt>0 else 0:.6e}"])

# 6) resumen
print(f"fase final={phi[-1]:.3f}  theta_rs final={thetas[-1]:.3f}")
print(f"|phi-theta|/t  último={abs(phi_s[-1]-thetas[-1])/ts_s[-1]:.3e}")
print(f"velocidad media={speed[1:-1].mean():.4f}  mediana curvatura={np.median(k[1:-1]):.4f}")
print(f"OK: {base}/fase_curvatura.csv, {base}/theta_comparacion.csv")
