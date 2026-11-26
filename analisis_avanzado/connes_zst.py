import sys
import numpy as np

# Prueba numerica de "Zeta Spectral Triples" (Connes, Consani, Moscovici, arXiv:2511.22755, nov 2025).
# Matriz de QW_lambda en la base V_n, |n| <= N, con las formulas (4.2)-(4.4) copiadas del articulo:
#   QW = W02 - WR - sum_p Wp ;  L = 2 ln(lambda) = ln(X), X = lambda^2 (primos y potencias <= X)
# Vector xi = autovector del menor autovalor; ceros previstos = raices reales de sum_n c_n/(z - 2 pi n/L)
# (transformada de Fourier de xi: xi^(z) ~ (e^{-izL} - 1) * sum_n c_n /(i(2 pi n/L - z)) ).
# Uso: python3 connes_zst.py X N ARCHIVO_CEROS [control]
#   X = lambda^2 (primos y potencias de primo <= X), N = tamano de la base (|n| <= N), por ejemplo 13 40 ceros_100k.txt
#   'control' baraja los pesos Lambda(k) entre los k (debe fallar). Subespacio par; doble precision (suelo ~1e-10).
X = float(sys.argv[1]) if len(sys.argv) > 1 else 13.0
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
L = np.log(X)
EG = 0.5772156649015329
ns = np.arange(-N, N + 1)

def q(m, n, y):
    """q(U_m,U_n)(y), y >= 0 (formula 2.9-2.10)"""
    m = np.asarray(m)[:, None, None]; n = np.asarray(n)[None, :, None]; y = np.asarray(y)[None, None, :]
    dif = (n - m).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        gen = (np.sin(2 * np.pi * m * y / L) - np.sin(2 * np.pi * n * y / L)) / (np.pi * dif)
    diag = 2 * (1 - y / L) * np.cos(2 * np.pi * n * y / L)
    return np.where(dif == 0, diag, gen)

# W_{0,2}
m_, n_ = np.meshgrid(ns, ns, indexing="ij")
W02 = 32 * L * np.sinh(L / 4) ** 2 * (L ** 2 - 16 * np.pi ** 2 * m_ * n_) / ((L ** 2 + 16 * np.pi ** 2 * m_ ** 2) * (L ** 2 + 16 * np.pi ** 2 * n_ ** 2))

# W_p : sobre potencias de primo k <= X
cr = np.ones(int(X) + 1, bool); cr[:2] = False
for i in range(2, int(X ** 0.5) + 1):
    if cr[i]: cr[i * i::i] = False
ks, lam = [], []
for p in np.nonzero(cr)[0]:
    k = int(p)
    while k <= X:
        ks.append(k); lam.append(np.log(p)); k *= int(p)
ks = np.array(ks, float); lam = np.array(lam)
if len(sys.argv) > 4 and sys.argv[4] == "control":
    lam = np.random.default_rng(3).permutation(lam)        # control: pesos de Lambda desordenados entre los k
Wp = np.einsum("k,mnk->mn", lam / np.sqrt(ks), q(ns, ns, np.log(ks)))

# W_R
xg, wg = np.polynomial.legendre.leggauss(600)
x = (xg + 1) * L / 2; w = wg * L / 2
om = q(ns, ns, x); om0 = q(ns, ns, np.array([0.0]))[:, :, 0]
integrando = (np.exp(x / 2)[None, None, :] * om - om0[:, :, None]) / (np.exp(x) - np.exp(-x))[None, None, :]
WR = om0 / 2 * (EG + np.log(4 * np.pi * (np.exp(L) - 1) / (np.exp(L) + 1))) + np.einsum("mnk,k->mn", integrando, w)

QW = W02 - WR - Wp
QW = (QW + QW.T) / 2
# subespacio par: c_n = c_-n
Pm = np.zeros((len(ns), N + 1)); Pm[N, 0] = 1
for j in range(1, N + 1):
    Pm[N + j, j] = Pm[N - j, j] = 1 / np.sqrt(2)
Qp = Pm.T @ QW @ Pm
evp, Vp = np.linalg.eigh(Qp)
ev, V = evp, Pm @ Vp
print(f"X = {X}, L = {L:.4f}, N = {N}: 5 menores autovalores de QW: {ev[:5]}")
c = V[:, 0]
print(f"simetria del vector: |c_n - c_-n| max = {np.abs(c - c[::-1]).max():.2e}  (paridad)")
a = 2 * np.pi * ns / L
# raices reales de g(z) = sum c_n/(z - a_n) en (0, zmax)
zmax = min(160.0, a[-1] - 1)
zz = np.arange(0.01, zmax, 0.002)
g = lambda z: np.sum(c[None, :] / (z[:, None] - a[None, :]), axis=1)
gv = g(zz)
cand = []
for i in range(len(zz) - 1):
    if gv[i] * gv[i + 1] < 0 and not np.any((a > zz[i]) & (a < zz[i + 1])):
        lo, hi = zz[i], zz[i + 1]
        for _ in range(60):
            mid = (lo + hi) / 2
            if g(np.array([lo]))[0] * g(np.array([mid]))[0] <= 0: hi = mid
            else: lo = mid
        cand.append((lo + hi) / 2)
cand = np.array(cand)
ref = np.loadtxt(sys.argv[3])[:200]
print(f"{len(cand)} raices reales positivas en (0, {zmax:.0f})")
zr = ref[ref < zmax - 1]
d = np.array([np.min(np.abs(cand - z)) for z in zr])
print(f"ceros de Riemann en el rango: {len(zr)}; coinciden a < 1e-6: {int(np.sum(d < 1e-6))}; a < 1e-3: {int(np.sum(d < 1e-3))}; a < 1e-1: {int(np.sum(d < 1e-1))}")
ok = np.nonzero(d < 1e-6)[0]
print("indices de los ceros reproducidos (a < 1e-6):", (ok + 1).tolist())
print("primeras diferencias (zero k, dif):", [(k + 1, f"{d[k]:.1e}") for k in range(min(10, len(d)))])
