import sys
import numpy as np

# Uso: python3 rs_completa.py ARCHIVO_CEROS [T0] [SIGMA]
# E17: ceros de Z(t) con la formula de Riemann-Siegel completa (suma hasta N=floor(sqrt(t/2pi)) + termino C0),
# frente a los ceros reales: error de posicion y silencio |S(x)| en x genericos.
ceros = np.loadtxt(sys.argv[1])
T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 40000.0
sig = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
a, b = T0 - 3 * sig, T0 + 3 * sig
vera = ceros[(ceros >= a) & (ceros <= b)]
esp = 2 * np.pi / np.log(T0 / (2 * np.pi))
gen = np.linspace(0.3, 2.6, 41)
gen = gen[np.min(abs(gen[:, None] - np.log(np.arange(2, 14))[None, :]), axis=1) > 0.05]


def S(z, xs):
    w = np.exp(-(z - T0) ** 2 / (2 * sig ** 2))
    return np.array([np.sum(w * np.cos(x * z)) for x in xs])


def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def Z_rs(t, corr=True):
    a_ = np.sqrt(t / (2 * np.pi))
    N = np.floor(a_).astype(int)
    p = a_ - N
    th = theta(t)
    s = np.zeros_like(t)
    for n in range(1, N.max() + 1):
        s += np.where(n <= N, np.cos(th - t * np.log(n)) / np.sqrt(n), 0.0)
    s *= 2
    if not corr:
        return s
    p2 = np.where(np.abs(np.cos(2 * np.pi * p)) < 1e-6, p + 1e-5, p)
    C0 = np.cos(2 * np.pi * (p2 ** 2 - p2 - 1 / 16)) / np.cos(2 * np.pi * p2)
    return s + (-1.0) ** (N - 1) * (2 * np.pi / t) ** 0.25 * C0


def raices(t, Z):
    i = np.nonzero(Z[:-1] * Z[1:] < 0)[0]
    return t[i] - Z[i] * (t[i + 1] - t[i]) / (Z[i + 1] - Z[i])


t = np.arange(a - 5, b + 5, 0.002)
print(f"Ventana T0={T0:.0f}, sigma={sig:.0f}: {len(vera)} ceros reales, espaciado medio {esp:.3f}, N_RS = {int(np.sqrt(T0/(2*np.pi)))}")
print(f"{'variante':>26s} {'ceros':>6s} {'|S| gener.':>11s} {'err mediano':>12s} {'err max':>9s} {'S(ln2)':>8s} {'S(ln5)':>8s}")
for nombre, corr in (("solo suma principal", False), ("suma + C0 (completa)", True)):
    z = raices(t, Z_rs(t, corr))
    z = z[(z >= a) & (z <= b)]
    k = np.searchsorted(z, vera)
    d = np.minimum(abs(z[np.clip(k, 0, len(z) - 1)] - vera), abs(z[np.clip(k - 1, 0, len(z) - 1)] - vera)) / esp
    print(f"{nombre:>26s} {len(z):6d} {np.median(abs(S(z, gen))):11.3f} {np.median(d):12.2e} {d.max():9.2e} {S(z, [np.log(2)])[0]:8.2f} {S(z, [np.log(5)])[0]:8.2f}")
print(f"{'ceros reales':>26s} {len(vera):6d} {np.median(abs(S(vera, gen))):11.3f} {0.0:12.2e} {0.0:9.2e} {S(vera, [np.log(2)])[0]:8.2f} {S(vera, [np.log(5)])[0]:8.2f}")
