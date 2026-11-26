import sys
import numpy as np
import mpmath as mp

# E21: formula completa de Bogomolny-Keating para la correlacion de pares (arXiv:math/0602270, ecs. 9-13), sin expandir:
#   r2(e) = rho^2 + diag(e) + off(e);  R2(s) = r2(s/rho)/rho^2;  rho = ln(E/2pi)/(2pi)
#   diag(e) = -(1/4pi^2) d^2/de^2 [ ln|zeta(1+ie)|^2 + 2 sum_p sum_r (1-r)/(r^2 p^r) cos(e r ln p) ]
#   off(e)  = (1/4pi^2) |zeta(1+ie)|^2 exp(2 pi i rho e) prod_p [1 - (1-p^(ie))^2/(p-1)^2] + c.c.
# Uso: python3 bk_completa.py DIR_CEROS
D = sys.argv[1] if len(sys.argv) > 1 else "."
mp.mp.dps = 25
PMAX = 200000
cr = np.ones(PMAX + 1, bool); cr[:2] = False
for i in range(2, int(PMAX ** 0.5) + 1):
    if cr[i]: cr[i * i::i] = False
P = np.nonzero(cr)[0].astype(float); LP = np.log(P)
RMAX = 40


def r2_completa(e, rho):
    z = mp.zeta(1 + 1j * e); z1 = mp.zeta(1 + 1j * e, 1, 1); z2 = mp.zeta(1 + 1j * e, 1, 2)
    d2lnz2 = float(-2 * mp.re(z2 / z - (z1 / z) ** 2))          # d^2/de^2 ln|zeta(1+ie)|^2
    s = 0.0
    for r in range(2, RMAX):                                     # el termino r=1 se anula: (1-r)=0
        s += np.sum((1 - r) * LP ** 2 / P ** r * np.cos(e * r * LP))
    diag = -(1 / (4 * np.pi ** 2)) * (d2lnz2 - 2 * s)            # d^2/de^2 cos(e r lnp) = -(r lnp)^2 cos
    pi_ = np.prod(1 - (1 - np.exp(1j * e * LP)) ** 2 / (P - 1) ** 2)
    off = 2 * np.real((1 / (4 * np.pi ** 2)) * float(abs(z) ** 2) * np.exp(2j * np.pi * rho * e) * pi_)
    return rho ** 2 + diag + off


def R2_completa(s, E):
    rho = np.log(E / (2 * np.pi)) / (2 * np.pi)
    return np.array([r2_completa(x / rho, rho) for x in s]) / rho ** 2


LAM, C = 1.57314, 1.4720
def r2_gue(s): return 1 - (np.sin(np.pi * s) / (np.pi * s)) ** 2
def r2_bk1(s, E):
    L = np.log(E / (2 * np.pi)); rho = L / (2 * np.pi); a = 1 + C / L
    return r2_gue(s) - LAM / (np.pi * rho) ** 2 * np.sin(np.pi * a * s) ** 2
def theta(t): return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t)

SMAX, DS = 3.0, 0.1
def hist_r2(u):
    ed = np.arange(0, SMAX + 1e-9, DS); cnt = np.zeros(len(ed) - 1)
    hi = np.searchsorted(u, u + SMAX, side="right")
    for i in range(len(u)):
        cnt += np.histogram(u[i + 1:hi[i]] - u[i], ed)[0]
    return cnt / (len(u) * DS)

s = (np.arange(int(SMAX / DS)) + 0.5) * DS
print("Validacion: a E enorme la formula completa debe tender a GUE + correccion de orden rho^-2")
Ev = 1e40
print("  s    completa    GUE+corr(16-17)")
for x in [0.35, 1.05, 2.05]:
    print(f"  {x:.2f}  {R2_completa(np.array([x]), Ev)[0]:.6f}  {r2_bk1(np.array([x]), Ev)[0]:.6f}")
print()
z = np.loadtxt(f"{D}/ceros_100k.txt"); z2 = np.loadtxt(f"{D}/ceros_1millon.txt")
ventanas = [("3e4-7e4 (1e5)", z[30000:70000]), ("1.0-1.2e6", z2[1000000:1200000])]
print(f"{'ventana':>14s} {'rho':>5s} {'RMS GUE':>9s} {'RMS orden1':>11s} {'RMS completa':>13s}   chi2/bin  GUE  orden1  completa")
for nom, zz in ventanas:
    E = zz.mean(); u = theta(zz) / np.pi
    nb = 10; B = len(u) // nb
    bl = np.array([hist_r2(u[k * B:(k + 1) * B]) for k in range(nb)])
    r2 = hist_r2(u); err = bl.std(axis=0, ddof=1) / np.sqrt(nb)
    g, b1, bc = r2_gue(s), r2_bk1(s, E), R2_completa(s, E)
    rms = lambda m: np.sqrt(np.mean((r2 - m) ** 2)); chi = lambda m: np.sum(((r2 - m) / err) ** 2) / len(s)
    print(f"{nom:>14s} {np.log(E / (2 * np.pi)) / (2 * np.pi):5.2f} {rms(g):9.4f} {rms(b1):11.4f} {rms(bc):13.4f}   {chi(g):13.2f} {chi(b1):6.2f} {chi(bc):8.2f}")
    if nom.startswith("3e4"):
        print("   s    medido   GUE     orden1  completa  error")
        for i in range(0, len(s), 3): print(f"  {s[i]:.2f}  {r2[i]:.4f}  {g[i]:.4f}  {b1[i]:.4f}  {bc[i]:.4f}   +-{err[i]:.4f}")
