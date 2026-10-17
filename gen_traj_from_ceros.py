#!/usr/bin/env python3
"""Genera trayectoria.csv para un N dado leyendo ceros.csv existente.
Uso: python3 gen_traj_from_ceros.py 5000
"""
import sys, csv
from mpmath import mp, zeta, gamma, sin, cos, pi, re, im, log, exp, mpc, mpf

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
base = f"salida_pfi01_N{N}"
DENSIDAD = 15.0
DPS = 15  # bajamos dps solo para la viz; el CSV de ceros sigue siendo el bueno

# 1) leer t_max de ceros.csv
ts = []
with open(f"{base}/ceros.csv") as f:
    for row in csv.DictReader(f):
        ts.append(float(row["t"]))
t_max = max(ts)
n_pts = int(t_max * DENSIDAD)
print(f"N={N}  t_max={t_max:.3f}  puntos={n_pts}  dps={DPS}")

mp.dps = DPS
with open(f"{base}/trayectoria.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t", "Re(Gamma)", "Im(Gamma)", "Z(t)"])
    for i in range(n_pts + 1):
        t = t_max * i / n_pts
        s = mpc(0.5, t)
        g = zeta(s)
        # Z(t) = Re(e^{i theta(t)} zeta(1/2+it)), theta = Riemann-Siegel
        theta = mp.arg(gamma(mpf(0.25) + 1j * t / 2)) - (t / 2) * log(pi)
        Z = re(exp(1j * theta) * g)
        w.writerow([f"{t:.10f}", f"{re(g):.15g}", f"{im(g):.15g}", f"{Z:.15g}"])
        if i % 5000 == 0:
            print(f"  {i}/{n_pts}")
print("listo")
