import re, sys
import numpy as np

# E20: R2(s) de los ceros frente a GUE y frente a la primera correccion de Bogomolny-Keating
# (Bogomolny, Bohigas, Leboeuf, Monastra 2006, arXiv:math/0602270, ecs. 16-17):
#   R2(s) = 1 - sin^2(pi s)/(pi s)^2 - Lambda/(pi^2 rho^2) * sin^2(pi*alpha*s),
#   rho = ln(E/2pi)/(2pi), alpha = 1 + C/ln(E/2pi), Lambda = 1.57314, C = 1.4720.
# Uso: python3 r2_bk.py DIR_CEROS   (ceros_100k.txt, ceros_1millon.txt, zeros3/4/5.txt)
D = sys.argv[1] if len(sys.argv) > 1 else "."
LAM, C = 1.57314, 1.4720
SMAX, DS = 3.0, 0.1


def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t)


def r2_gue(s):
    return 1 - (np.sin(np.pi * s) / (np.pi * s)) ** 2


def r2_bk(s, E):
    L = np.log(E / (2 * np.pi)); rho = L / (2 * np.pi); alpha = 1 + C / L
    return r2_gue(s) - LAM / (np.pi * rho) ** 2 * np.sin(np.pi * alpha * s) ** 2


def hist_r2(u):
    """histograma de distancias de pares (u con espaciado medio 1), normalizado a R2"""
    n = len(u); ed = np.arange(0, SMAX + 1e-9, DS); cnt = np.zeros(len(ed) - 1)
    hi = np.searchsorted(u, u + SMAX, side="right")
    for i in range(n):
        d = u[i + 1:hi[i]] - u[i]
        cnt += np.histogram(d, ed)[0]
    return cnt / (n * DS)


def analiza(nom, u, E):
    nb = 10; B = len(u) // nb
    bl = np.array([hist_r2(u[k * B:(k + 1) * B]) for k in range(nb)])
    r2 = hist_r2(u); err = bl.std(axis=0, ddof=1) / np.sqrt(nb)
    s = (np.arange(len(r2)) + 0.5) * DS
    g, b = r2_gue(s), r2_bk(s, E)
    chi = lambda m: np.sum(((r2 - m) / err) ** 2) / len(s)
    print(f"{nom:>14s} E={E:9.2e} rho={np.log(E / (2 * np.pi)) / (2 * np.pi):5.2f}  RMS(GUE)={np.sqrt(np.mean((r2 - g) ** 2)):.4f}  RMS(BK)={np.sqrt(np.mean((r2 - b) ** 2)):.4f}"
          f"   chi2/bin  GUE: {chi(g):6.2f}  BK: {chi(b):6.2f}")
    return s, r2, g, b, err


z = np.loadtxt(f"{D}/ceros_100k.txt")
res = []
u = theta(z[30000:70000]) / np.pi; res.append(analiza("1e5: 3e4-7e4", u, z[30000:70000].mean()))
z2 = np.loadtxt(f"{D}/ceros_1millon.txt")
u = theta(z2[1000000:1200000]) / np.pi; res.append(analiza("1e6..1.2e6", u, z2[1000000:1200000].mean()))
for nom, f in [("zeros 1e12", "zeros3.txt"), ("zeros 1e21", "zeros4.txt"), ("zeros 1e22", "zeros5.txt")]:
    L = open(f"{D}/{f}").read().splitlines()
    off = float(re.search(r"gamma\s*-\s*(\d+)", " ".join(L[:3])).group(1))
    v = np.array([float(l) for l in L if re.fullmatch(r"\s*[-+]?\d+(\.\d+)?\s*", l)])
    u = np.concatenate([[0], np.cumsum(np.diff(v) * np.log((off + v[:-1]) / (2 * np.pi)) / (2 * np.pi))])
    res.append(analiza(nom, u, off))
s, r2, g, b, err = res[0]
print("\nDetalle primera ventana: s, R2 medido, GUE, BK, error")
for i in range(0, len(s), 3):
    print(f"  {s[i]:.2f}  {r2[i]:.4f}  {g[i]:.4f}  {b[i]:.4f}  +-{err[i]:.4f}")


# --- Amplitud de la correccion frente a la altura ---
print("ajuste de kappa en DeltaR2 = -kappa*sin^2(pi*alpha*s), s<=1.5  (BK: kappa = Lambda/(pi rho)^2)")
print(f"{'indices':>16s} {'E':>9s} {'rho':>5s} {'kappa_med':>10s} {'+-':>6s} {'kappa_BK':>9s} {'cociente':>9s}")
for a, b in [(5000, 25000), (25000, 65000), (100000, 200000), (400000, 600000), (1000000, 1300000), (1700000, 2000000)]:
    zz = z2[a:b]; u = theta(zz) / np.pi
    E = zz.mean(); L = np.log(E / (2 * np.pi)); rho = L / (2 * np.pi); al = 1 + C / L
    nb = 8; B = len(u) // nb
    s = (np.arange(int(SMAX / DS)) + 0.5) * DS; m = s <= 1.5
    def kap(r2):
        x = np.sin(np.pi * al * s[m]) ** 2; y = r2_gue(s[m]) - r2[m]
        return np.sum(x * y) / np.sum(x * x)
    k = kap(hist_r2(u)); ks = np.array([kap(hist_r2(u[i * B:(i + 1) * B])) for i in range(nb)])
    kb = LAM / (np.pi * rho) ** 2
    print(f"{a:7d}-{b:<8d} {E:9.2e} {rho:5.2f} {k:10.4f} {ks.std(ddof=1)/np.sqrt(nb):6.4f} {kb:9.4f} {k/kb:9.2f}")
