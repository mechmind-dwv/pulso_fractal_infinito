import sys
import numpy as np

# Uso: python3 silencio.py ARCHIVO_CEROS [T0] [SIGMA]
# E16: el "silencio" |S(x)| en x genericos como medida de la precision de posicion de los niveles.
#  1) ceros reales desplazados con ruido de eps espaciados
#  2) ceros aproximados de Z_M (M enteros de Riemann-Siegel) en la misma ventana
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


def ceros_Z(M, t, th):
    Z = np.zeros_like(t)
    for n in range(1, M + 1):
        Z += 2 * np.cos(th - t * np.log(n)) / np.sqrt(n)
    i = np.nonzero(Z[:-1] * Z[1:] < 0)[0]
    return t[i] - Z[i] * (t[i + 1] - t[i]) / (Z[i + 1] - Z[i])


rng = np.random.default_rng(0)
print(f"Ventana T0={T0:.0f}, sigma={sig:.0f}: {len(vera)} ceros, espaciado medio {esp:.3f}")
print("\n1) Ceros reales con ruido de posicion (5 repeticiones)")
print(f"{'eps (espaciados)':>17s} {'|S| generico':>13s}")
for eps in (0, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 0.4):
    r = [np.median(abs(S(vera + eps * esp * rng.normal(size=len(vera)), gen))) for _ in range(5)]
    print(f"{eps:17.0e} {np.mean(r):13.3f}")
t = np.arange(a - 5, b + 5, 0.005)
th = theta(t)
N_RS = int(np.sqrt(T0 / (2 * np.pi)))
print(f"\n2) Ceros de Z_M con M enteros (N_RS = {N_RS}): silencio, error de posicion y picos")
print(f"{'M':>5s} {'|S| generico':>13s} {'err rms (esp.)':>15s} {'S(ln2)':>8s} {'S(ln5)':>8s}")
for M in (5, 20, N_RS // 2, N_RS, int(1.5 * N_RS)):
    z = ceros_Z(M, t, th)
    z = z[(z >= a) & (z <= b)]
    k = np.searchsorted(z, vera)
    d = np.minimum(abs(z[np.clip(k, 0, len(z) - 1)] - vera), abs(z[np.clip(k - 1, 0, len(z) - 1)] - vera))
    err = np.sqrt(np.mean(np.minimum(d, esp) ** 2)) / esp
    print(f"{M:5d} {np.median(abs(S(z, gen))):13.3f} {err:15.3f} {S(z, [np.log(2)])[0]:8.2f} {S(z, [np.log(5)])[0]:8.2f}")
