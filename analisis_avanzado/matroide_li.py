import sys, time, random
import mpmath as mp
# Matroide de dependencia lineal racional de los gamma_n (hipotesis LI).
# Uso: python3 matroide_li.py [Nceros=40] [k=8] [dps=80] [maxcoef=1e5]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
k = int(sys.argv[2]) if len(sys.argv) > 2 else 8
dps = int(sys.argv[3]) if len(sys.argv) > 3 else 80
mc = int(float(sys.argv[4])) if len(sys.argv) > 4 else 10**5
mp.mp.dps = dps
t0 = time.time()
g = [mp.zetazero(n).imag for n in range(1, N + 1)]
print(f"{N} ceros a {dps} digitos ({time.time()-t0:.0f} s); ventanas de {k}, coef <= {mc}")
def busca(v):
    return mp.pslq(v, tol=mp.mpf(10) ** (-(dps - 15)), maxcoeff=mc, maxsteps=100000)
# 1) ventanas consecutivas y subconjuntos al azar
rel = 0; tot = 0
for i in range(0, N - k + 1):
    tot += 1
    if busca(g[i:i + k]) is not None: rel += 1
random.seed(1)
for _ in range(60):
    S = sorted(random.sample(range(N), k)); tot += 1
    if busca([g[j] for j in S]) is not None: rel += 1
print(f"zeros: {tot} subconjuntos de {k}, relaciones enteras halladas: {rel}")
# 2) control: reales al azar (no deben dar relacion) y relacion plantada
r = [mp.rand() * 50 for _ in range(k)]
print("control azar, relaciones:", busca(r) is not None)
v = list(g[:k - 1]); v.append(3 * g[0] - 2 * g[1] + 5 * g[2])   # plantada
rp = busca(v)
print("control plantada (3,-2,5,...), detectada:", rp)
# 3) hay relaciones triviales entre n y 2n? comprobacion de rango sobre Q del conjunto completo con coef pequenos
for kk in (4, 6, 8, 10, 12):
    if kk <= N:
        print(f"primeros {kk}: relacion?", busca(g[:kk]) is not None)
