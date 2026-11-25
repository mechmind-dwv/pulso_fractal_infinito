import re, sys
import numpy as np

# E19: los ceros frente a CUE de tamano efectivo N_eff = ln(E/2pi)/sqrt(12*Lambda), Lambda = 1.57314
# (Bogomolny, Bohigas, Leboeuf, Monastra 2006, arXiv:math/0602270). Uso: python3 neff_cue.py DIR_CEROS
# DIR_CEROS contiene ceros_100k.txt, ceros_1millon.txt, zeros3.txt, zeros4.txt, zeros5.txt
D = sys.argv[1] if len(sys.argv) > 1 else "."
LAM = 1.57314
rng = np.random.default_rng(1)


def neff(E):
    return np.log(E / (2 * np.pi)) / np.sqrt(12 * LAM)


def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t)


def estad(s):
    s1, s2 = s[:-1], s[1:]
    r = np.minimum(s1, s2) / np.maximum(s1, s2)
    return r.mean(), (s < 0.3).mean(), np.corrcoef(s1, s2)[0, 1]


def con_error(s, nb=10):
    """estadistica y error estandar entre nb bloques consecutivos"""
    B = len(s) // nb
    v = np.array([estad(s[k * B:(k + 1) * B]) for k in range(nb)])
    return estad(s), v.std(axis=0, ddof=1) / np.sqrt(nb)


def cue(N, M):
    A = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(A)
    d = np.diagonal(R, axis1=1, axis2=2)
    Q = Q * (d / np.abs(d))[:, None, :]
    th = np.sort(np.angle(np.linalg.eigvals(Q)), axis=1)
    s = np.diff(np.concatenate([th, th[:, :1] + 2 * np.pi], axis=1), axis=1) * N / (2 * np.pi)
    s1, s2 = s, np.roll(s, -1, axis=1)
    r = np.minimum(s1, s2) / np.maximum(s1, s2)
    return r.mean(), (s < 0.3).mean(), np.corrcoef(s1.ravel(), s2.ravel())[0, 1]


print("CUE(N): <r>, P(s<0.3), corr(s_n, s_n+1)")
tabla = {}
for N, M in [(2, 400000), (3, 400000), (4, 300000), (5, 200000), (6, 150000), (8, 80000), (11, 40000), (16, 25000), (24, 15000), (50, 5000)]:
    tabla[N] = cue(N, M)
    r, p, c = tabla[N]
    print(f"  N={N:3d}   {r:.4f}   {100 * p:5.2f}%   {c:+.3f}")
Ns = sorted(tabla)
print("  (N=inf, GUE asintotico: <r> = 0.5996; corr ~ -0.308)\n")

print("Ceros: espaciados normalizados con la densidad local")
print(f"{'origen':>16s} {'t':>9s} {'N_eff':>6s} {'<r>':>14s} {'P(s<0.3)':>13s} {'corr':>16s}  corr CUE(N_eff) interpolada")
fuentes = []
z = np.loadtxt(f"{D}/ceros_100k.txt")
fuentes.append(("primeros 1e5", z, z.mean()))
z = np.loadtxt(f"{D}/ceros_1millon.txt")
fuentes.append(("ceros 1e6..2e6", z[1000000:], z[1000000:].mean()))
for nom, f in [("zeros 1e12", "zeros3.txt"), ("zeros 1e21", "zeros4.txt"), ("zeros 1e22", "zeros5.txt")]:
    L = open(f"{D}/{f}").read().splitlines()
    off = float(re.search(r"gamma\s*-\s*(\d+)", " ".join(L[:3])).group(1))
    v = np.array([float(l) for l in L if re.fullmatch(r"\s*[-+]?\d+(\.\d+)?\s*", l)])
    fuentes.append((nom, (off, v), off))
for nom, dat, t0 in fuentes:
    if isinstance(dat, tuple):
        off, v = dat
        s = np.diff(v) * np.log((off + v[:-1]) / (2 * np.pi)) / (2 * np.pi)
    else:
        s = np.diff(theta(dat)) / np.pi
    (r, p, c), (er, ep, ec) = con_error(s)
    ne = neff(t0)
    cc = np.interp(ne, Ns, [tabla[n][2] for n in Ns])
    print(f"{nom:>16s} {t0:9.2e} {ne:6.2f} {r:.4f}+-{er:.4f} {100 * p:5.2f}+-{100 * ep:.2f}% {c:+.3f}+-{ec:.3f}  {cc:+.3f}")
