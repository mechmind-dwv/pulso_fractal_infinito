#!/usr/bin/env python3
"""Controles de §13. Genera tres secuencias sintéticas con la misma N que
salida_pfi01_N{N}/separaciones.csv:
  1) GUE muestreado (eigenvalores de matriz hermítica aleatoria)
  2) Poisson (gaps exponenciales)
  3) Gaps de ceros barajados (destruye correlación, preserva marginales)
Luego compara P(s), <r>, ACF(1) con los ceros reales.

Uso: python3 lab_controls.py 1000
"""
import sys, csv, os
import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
base = f"salida_pfi01_N{N}"

# 1) cargar gaps reales
g_real = []
with open(f"{base}/separaciones.csv") as f:
    r = csv.DictReader(f)
    key = r.fieldnames[0] if r.fieldnames else None
    for row in r:
        g_real.append(float(row[key]))
g_real = np.asarray(g_real)

# 2) GUE: eigenvalores de una matriz hermítica aleatoria
def gue_gaps(n):
    m = n + 1
    A = np.random.randn(m, m) + 1j*np.random.randn(m, m)
    A = (A + A.conj().T) / 2.0  # hermítica
    ev = np.sort(np.linalg.eigvalsh(A).real)
    gaps = np.diff(ev)
    # normalizar por densidad local media
    return gaps / gaps.mean()

# 3) Poisson
def poisson_gaps(n):
    return np.random.exponential(1.0, n)

# 4) gaps de ceros barajados
def shuffled(g):
    out = g.copy(); np.random.shuffle(out); return out

# 5) estadísticos
def P_s(g, bins=40, rng=(0, 4)):
    g = g / g.mean()
    c, e = np.histogram(g, bins=bins, range=rng)
    return (e[:-1] + e[1:]) / 2, c / c.sum()

def r_ratio(g):
    g = g / g.mean()
    r = np.minimum(g[:-1], g[1:]) / np.maximum(g[:-1], g[1:])
    return float(r.mean())

def acf1(g, lag=1):
    g = (g - g.mean()) / g.std()
    return float((g[:-lag] * g[lag:]).mean())

def P_lt_03(g):
    g = g / g.mean()
    return float((g < 0.3).mean())

np.random.seed(42)
n = len(g_real)
controls = {
    "zeros_reales":  g_real,
    "gue":           gue_gaps(n),
    "poisson":       poisson_gaps(n),
    "shuffled":      shuffled(g_real),
}

print(f"N={N}  n_gaps={n}")
print(f"{'secuencia':<15} {'<r>':>8} {'ACF(1)':>9} {'P(s<0.3)':>10}")
print("-" * 46)
rows = []
for name, g in controls.items():
    r_m = r_ratio(g); a1 = acf1(g); p03 = P_lt_03(g)
    print(f"{name:<15} {r_m:>8.4f} {a1:>9.4f} {p03:>10.4f}")
    rows.append([name, f"{r_m:.6f}", f"{a1:.6f}", f"{p03:.6f}", str(len(g))])

os.makedirs(f"{base}", exist_ok=True)
with open(f"{base}/controles.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["secuencia", "r_ratio", "acf1", "P_lt_03", "n"])
    w.writerows(rows)

# histogramas superpuestos
with open(f"{base}/controles_Ps.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["s", "zeros", "gue", "poisson", "shuffled"])
    all_P = {name: P_s(g) for name, g in controls.items()}
    s_grid = all_P["zeros_reales"][0]
    for i in range(len(s_grid)):
        w.writerow([f"{s_grid[i]:.4f}"] +
                   [f"{all_P[k][1][i]:.6f}" for k in ("zeros_reales","gue","poisson","shuffled")])

print(f"OK: {base}/controles.csv, {base}/controles_Ps.csv")
