import sys
import numpy as np

# Uso: python3 ceros_desde_enteros.py ARCHIVO_CEROS [T0] [ANCHO]
# Ceros aproximados de Z(t) usando solo los primeros M enteros:
#   Z_M(t) = 2 * sum_{n<=M} cos(theta(t) - t ln n) / sqrt(n)
ceros = np.loadtxt(sys.argv[1])
T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 10000.0
W = float(sys.argv[3]) if len(sys.argv) > 3 else 1500.0


def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def razon(z):
    s = np.diff(z)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r.mean()


def ceros_Z(M, t, th):
    Z = np.zeros_like(t)
    for n in range(1, M + 1):
        Z += 2 * np.cos(th - t * np.log(n)) / np.sqrt(n)
    i = np.nonzero(Z[:-1] * Z[1:] < 0)[0]
    return t[i] - Z[i] * (t[i + 1] - t[i]) / (Z[i + 1] - Z[i])


t = np.arange(T0, T0 + W, 0.005)
th = theta(t)
vera = ceros[(ceros >= T0) & (ceros < T0 + W)]
dens = np.log(T0 / (2 * np.pi)) / (2 * np.pi)
N_RS = int(np.sqrt((T0 + W / 2) / (2 * np.pi)))
print(f"Ventana [{T0:.0f}, {T0+W:.0f}]: {len(vera)} ceros reales, <r> real = {razon(vera):.4f}")
print(f"Entero de corte de Riemann-Siegel: floor(sqrt(t/2pi)) = {N_RS}")
print(f"{'M':>5s} {'M/N_RS':>7s} {'ceros':>6s} {'encontrados':>12s} {'error rms/espaciado':>20s} {'<r>':>7s}")
for M in sorted({1, 2, 3, 5, 8, 13, 20, N_RS // 2, N_RS, int(1.5 * N_RS), 2 * N_RS}):
    z = ceros_Z(M, t, th)
    if len(z) < 10:
        print(f"{M:5d} {M/N_RS:7.2f} {len(z):6d}")
        continue
    k = np.searchsorted(z, vera)
    d = np.minimum(abs(z[np.clip(k, 0, len(z) - 1)] - vera), abs(z[np.clip(k - 1, 0, len(z) - 1)] - vera))
    enc = np.mean(d < 0.25 / dens)
    print(f"{M:5d} {M/N_RS:7.2f} {len(z):6d} {enc:12.3f} {np.sqrt(np.mean(np.minimum(d, 1/dens) ** 2)) * dens:20.3f} {razon(z):7.4f}")
