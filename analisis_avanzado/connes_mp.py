import sys, time
import numpy as np
import mpmath as mp

# E23 con precision multiple: misma construccion que connes_zst.py (QW_lambda en la base V_n, subespacio par),
# pero con mpmath (dps digitos). Las integrales se reducen a integrales unidimensionales por indice:
#   m != n:  WR = [J(m) - J(n)] / (pi (n - m)),   J(n) = int_0^L e^{x/2} sin(2 pi n x/L) / (e^x - e^-x) dx
#   m == n:  WR = gamma + ln(4pi (e^L-1)/(e^L+1)) + D(n),  D(n) = int_0^L (e^{x/2} 2(1-x/L) cos(2 pi n x/L) - 2)/(e^x - e^-x) dx
# Uso: python3 connes_mp.py X N [dps]     (por ejemplo 13 40 60)
X = int(sys.argv[1]) if len(sys.argv) > 1 else 13
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
t0 = time.time()
L = mp.log(X); pi = mp.pi

def den(x): return mp.exp(x) - mp.exp(-x)
J = {n: mp.quad(lambda x: mp.exp(x / 2) * mp.sin(2 * pi * n * x / L) / den(x), [0, L / 2, L]) for n in range(0, 2 * N + 1)}
Jf = lambda n: J[n] if n >= 0 else -J[-n]
Dd = {n: mp.quad(lambda x: (mp.exp(x / 2) * 2 * (1 - x / L) * mp.cos(2 * pi * n * x / L) - 2) / den(x), [0, L / 2, L]) for n in range(0, N + 1)}
cte = mp.euler + mp.log(4 * pi * (mp.exp(L) - 1) / (mp.exp(L) + 1))

# potencias de primo <= X
ks, lam = [], []
for p in [q for q in range(2, X + 1) if all(q % r for r in range(2, int(q ** 0.5) + 1))]:
    k = p
    while k <= X:
        ks.append(k); lam.append(mp.log(p)); k *= p

def qf(m, n, y):
    if m == n: return 2 * (1 - y / L) * mp.cos(2 * pi * n * y / L)
    return (mp.sin(2 * pi * m * y / L) - mp.sin(2 * pi * n * y / L)) / (pi * (n - m))

def QW(m, n):
    W02 = 32 * L * mp.sinh(L / 4) ** 2 * (L ** 2 - 16 * pi ** 2 * m * n) / ((L ** 2 + 16 * pi ** 2 * m ** 2) * (L ** 2 + 16 * pi ** 2 * n ** 2))
    Wp = sum(l / mp.sqrt(k) * qf(m, n, mp.log(k)) for k, l in zip(ks, lam))
    WR = (cte + Dd[abs(n)]) if m == n else (Jf(m) - Jf(n)) / (pi * (n - m))
    return W02 - WR - Wp

# bloque par
Q = mp.matrix(N + 1, N + 1)
for j in range(N + 1):
    for k in range(j, N + 1):
        if j == 0 and k == 0: v = QW(0, 0)
        elif j == 0: v = mp.sqrt(2) * QW(0, k)
        else: v = QW(j, k) + QW(j, -k)
        Q[j, k] = Q[k, j] = v
print(f"matriz construida en {time.time() - t0:.0f} s", flush=True)
E, V = mp.eigsy(Q)
ev = [E[i] for i in range(N + 1)]
order = sorted(range(N + 1), key=lambda i: ev[i])
print("menores autovalores de QW (par):", [mp.nstr(ev[i], 6) for i in order[:6]])
i0 = order[0]
v = [V[j, i0] for j in range(N + 1)]
c = {0: v[0]}
for j in range(1, N + 1): c[j] = c[-j] = v[j] / mp.sqrt(2)
a = {n: 2 * pi * n / L for n in range(-N, N + 1)}
def g(z): return mp.fsum(c[n] / (z - a[n]) for n in range(-N, N + 1))

# raices reales positivas: barrido en doble precision, refinado en mp
zmax = float(a[N]) - 1
zz = [mp.mpf(k) / 100 for k in range(1, int(zmax * 100))]
raices = []
gprev = g(zz[0])
for i in range(len(zz) - 1):
    gnext = g(zz[i + 1])
    if gprev * gnext < 0 and not any(zz[i] < a[n] < zz[i + 1] for n in range(0, N + 1)):
        try:
            r = mp.findroot(g, (zz[i], zz[i + 1]), solver="anderson", tol=mp.mpf(10) ** (-mp.mp.dps // 2), maxsteps=200, verify=False)
            if zz[i] <= r <= zz[i + 1]: raices.append(r)
        except Exception: pass
    gprev = gnext
raices = sorted(raices)
uniq = []
for r in raices:
    if not uniq or abs(r - uniq[-1]) > mp.mpf(10) ** -8: uniq.append(r)
raices = uniq
print(f"{len(raices)} raices reales positivas en (0, {zmax:.0f})  (t = {time.time() - t0:.0f} s)")
nz = int(sum(1 for k in range(1, 200) if float(mp.zetazero(k).imag) < zmax - 1)) if False else None
print("  k   cero de Riemann (mp.zetazero)        raiz mas cercana - cero     coincide?")
zeros = []
k = 1
while True:
    z = mp.zetazero(k).imag
    if z > zmax - 1: break
    zeros.append(z); k += 1
usadas = set()
for k, z in enumerate(zeros, 1):
    d, j = min(((abs(r - z), j) for j, r in enumerate(raices)), key=lambda t: t[0])
    usadas.add(j)
    if k <= 12 or k % 5 == 0 or d > 1e-6:
        print(f"{k:3d}   {mp.nstr(z, 15):>22s}   {mp.nstr(raices[j] - z, 3):>22s}")
print(f"ceros en el rango: {len(zeros)}; reproducidos a < 1e-30: {sum(1 for z in zeros if min(abs(r - z) for r in raices) < mp.mpf(10) ** -30)}; a < 1e-6: {sum(1 for z in zeros if min(abs(r - z) for r in raices) < mp.mpf(10) ** -6)}")
esp = [r for j, r in enumerate(raices) if j not in usadas]
print(f"raices que no son ceros de Riemann (espurias): {len(raices) - len(usadas)}", [mp.nstr(r, 6) for r in esp][:12])
