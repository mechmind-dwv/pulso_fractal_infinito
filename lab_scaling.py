#!/usr/bin/env python3
"""Autosimilitud: compara P(s) entre N=1000 y N=5000 en la ventana común.
Uso: python3 lab_scaling.py
"""
import csv, os
import numpy as np

def load_gaps(N):
    path = f"salida_pfi01_N{N}/separaciones.csv"
    if not os.path.exists(path): return None
    g = []
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r: g.append(float(row[r.fieldnames[0]]))
    return np.asarray(g)

def P_s(g, bins=40, rng=(0, 4)):
    g = g / g.mean()
    c, e = np.histogram(g, bins=bins, range=rng)
    return (e[:-1] + e[1:]) / 2, c / c.sum()

g1 = load_gaps(1000); g5 = load_gaps(5000)
if g1 is None or g5 is None:
    print("faltan N=1000 o N=5000"); raise SystemExit(1)

s, P1 = P_s(g1)
_, P5 = P_s(g5)
sup_diff = float(np.abs(P1 - P5).max())
l1_diff  = float(np.abs(P1 - P5).sum())
print(f"N=1000 (n={len(g1)})  N=5000 (n={len(g5)})")
print(f"sup|P1-P5|={sup_diff:.4f}   L1={l1_diff:.4f}")
print(f"criterio §8: si NO decrece respecto a N=500 vs N=1000, autosimilitud falla")

with open("salida_pfi01_N1000/scaling_1000_5000.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["s", "P1000", "P5000", "diff"])
    for i in range(len(s)):
        w.writerow([f"{s[i]:.4f}", f"{P1[i]:.6f}", f"{P5[i]:.6f}", f"{abs(P1[i]-P5[i]):.6f}"])
print("OK: salida_pfi01_N1000/scaling_1000_5000.csv")
