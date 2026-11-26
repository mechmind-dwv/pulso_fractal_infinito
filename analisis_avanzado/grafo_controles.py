import sys, time
import numpy as np
# Grafo cuantico completo K_V. Longitudes de arista: ln(p) (primos) frente a controles.
# Niveles = ceros EXACTOS de f(k)=prod sin(theta_j/2) (real), por cambios de signo; control con ley de Weyl.
# Uso: python3 grafo_controles.py [V=6] [K=40] [realizaciones=12]
V = int(sys.argv[1]) if len(sys.argv) > 1 else 6
K = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
R = int(sys.argv[3]) if len(sys.argv) > 3 else 12
edges = [(i, j) for i in range(V) for j in range(i+1, V)]
B = len(edges)
def primos(n):
    out, p = [], 2
    while len(out) < n:
        if all(p % q for q in out): out.append(p)
        p += 1
    return np.array(out, float)
def haar(n, rng):
    A = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    Q, Rr = np.linalg.qr(A)
    return Q * (np.diag(Rr)/np.abs(np.diag(Rr)))
def coe(n, rng):
    W = haar(n, rng)
    return W.T @ W
def matriz_S(rng, trs):
    S = np.zeros((2*B, 2*B), complex)
    for v in range(V):
        inc = []
        for e, (u, w) in enumerate(edges):
            if u == v: inc.append((2*e+1, 2*e))
            elif w == v: inc.append((2*e, 2*e+1))
        U = coe(len(inc), rng) if trs else haar(len(inc), rng)
        for i, (_, oi) in enumerate(inc):
            for j, (ij, _) in enumerate(inc):
                S[oi, ij] = U[i, j]
    return S
def niveles(S, L, K, dk):
    n = 2*B
    Ld = np.repeat(L, 2)
    ks = np.arange(dk, K, dk)
    phi0 = np.angle(np.linalg.det(S))
    f = np.empty(len(ks))
    for i, k in enumerate(ks):
        Sig = S * np.exp(1j*k*Ld)[None, :]
        d = np.linalg.det(np.eye(n) - Sig)
        z = d*np.exp(-1j*(phi0 + k*Ld.sum())/2)/((-2j)**n)
        f[i] = z.real
    idx = np.where(f[:-1]*f[1:] < 0)[0]
    return ks[idx] - f[idx]*dk/(f[idx+1]-f[idx])
def razon(E):
    s = np.diff(E); r = np.minimum(s[:-1], s[1:])/np.maximum(s[:-1], s[1:])
    return r.mean()
rng = np.random.default_rng(7)
P = np.log(primos(B))
casos = [("ln p (primos)", lambda: P, False),
         ("longitudes al azar (misma media)", lambda: rng.uniform(0.3, 2*P.mean()-0.3, B), False),
         ("ln p, S simetrica (TRS)", lambda: P, True),
         ("ln n, n=2..B+1 (dependientes)", lambda: np.log(np.arange(2, B+2, dtype=float)), False)]
print(f"K_{V}, B={B} aristas, k en (0,{K}), {R} realizaciones por caso")
print("GUE 0.5996 | GOE 0.5307 | Poisson 0.3863")
for nombre, gen, trs in casos:
    rs, wey = [], []
    for r in range(R):
        L = gen()
        S = matriz_S(rng, trs)
        E = niveles(S, L, K, 0.003)
        rs.append(razon(E)); wey.append(len(E)/(K*L.sum()/np.pi))
    rs = np.array(rs)
    print(f"{nombre:34s} <r> = {rs.mean():.4f} +- {rs.std()/np.sqrt(R):.4f}   niveles/Weyl = {np.mean(wey):.3f}")
