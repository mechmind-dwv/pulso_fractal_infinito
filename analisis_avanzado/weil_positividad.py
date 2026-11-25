import sys
import numpy as np

# Uso: python3 weil_positividad.py ceros_100k.txt [s] [separacion]   (s=3 sep=6 | s=1.5 sep=5 | s=0.8 sep=3)
Z = np.loadtxt(sys.argv[1])
s = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
sep = float(sys.argv[3]) if len(sys.argv) > 3 else 6.0
rj = np.arange(12, 54, sep)
N = len(rj)

trapecio = getattr(np, "trapezoid", None) or np.trapz
B2K = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6]


def digamma(z):
    z = np.asarray(z, complex)
    out = np.zeros_like(z)
    while np.any(z.real < 10):
        m = z.real < 10
        out[m] -= 1 / z[m]
        z = np.where(m, z + 1, z)
    r = np.log(z) - 1 / (2 * z)
    for k, b in enumerate(B2K, start=1):
        r -= b / (2 * k * z ** (2 * k))
    return out + r


umax = np.sqrt(80.0) / (s / 2)
nmax = int(np.exp(min(umax, 16))) + 1
sv = np.ones(nmax + 1, bool)
sv[:2] = False
for i in range(2, int(nmax ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = False
lam = np.zeros(nmax + 1)
for p in np.nonzero(sv)[0]:
    q = p
    while q <= nmax:
        lam[q] = np.log(p)
        q *= p
n = np.nonzero(lam)[0]
L = lam[n]
u = np.log(n)


def h(r, m, C):
    return C * np.exp(-(r - m) ** 2 / s ** 2)


Wp = np.zeros((N, N))
Wz = np.zeros((N, N))
gz = np.concatenate([Z, -Z])
for j in range(N):
    for k in range(j, N):
        m = (rj[j] + rj[k]) / 2
        C = np.exp(-(rj[j] - rj[k]) ** 2 / (4 * s ** 2))
        g0 = C * s / (2 * np.sqrt(np.pi))
        polo = 2 * np.real(h(0.5j, m, C))
        r = np.linspace(m - 12 * s, m + 12 * s, 40001)
        arch = trapecio(h(r, m, C) * digamma(0.25 + 0.5j * r).real, r) / (2 * np.pi)
        pr = -np.sum(L / np.sqrt(n) * 2 * g0 * np.cos(m * u) * np.exp(-s ** 2 * u ** 2 / 4))
        Wp[j, k] = Wp[k, j] = polo - g0 * np.log(np.pi) + arch + pr
        Wz[j, k] = Wz[k, j] = np.sum(h(gz, m, C))
print(f"s={s}, {N} funciones de prueba r_j={rj[0]:.0f}..{rj[-1]:.0f}; primos hasta n={nmax}")
print(f"max|W_primos - W_ceros| = {np.abs(Wp - Wz).max():.2e}   (escala W ~ {np.abs(Wz).max():.2f})")
print(f"autovalores minimos: primos {np.linalg.eigvalsh(Wp).min():.3e}   ceros {np.linalg.eigvalsh(Wz).min():.3e}")
i0 = np.argmin(abs(Z - 31.0))
g0_ = Z[i0]
print(f"\nControl sintetico: el cero gamma={g0_:.3f} se sustituye por un cuarteto fuera de la recta (peso 1/2 c/u)")
print(" delta   min autovalor")
zs = np.concatenate([np.delete(Z, i0), -np.delete(Z, i0)])
for d in (0.0, 0.05, 0.1, 0.2, 0.4):
    W = np.zeros((N, N))
    for j in range(N):
        for k in range(j, N):
            m = (rj[j] + rj[k]) / 2
            C = np.exp(-(rj[j] - rj[k]) ** 2 / (4 * s ** 2))
            v = np.sum(h(zs, m, C))
            for gg in (g0_ + 1j * d, g0_ - 1j * d, -g0_ + 1j * d, -g0_ - 1j * d):
                v = v + 0.5 * h(gg, m, C)
            W[j, k] = W[k, j] = np.real(v)
    print(f" {d:5.2f}  {np.linalg.eigvalsh(W).min():+.3e}")
