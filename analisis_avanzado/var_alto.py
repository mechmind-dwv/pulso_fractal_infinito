import re, sys
import numpy as np
D = sys.argv[1] if len(sys.argv) > 1 else "."
GUE_V, GUE_R = 0.1800, 0.5996
def carga(f):
    L = open(f"{D}/{f}").read().splitlines()
    off = float(re.search(r"gamma\s*-\s*(\d+)", " ".join(L[:3])).group(1))
    v = np.array([float(l) for l in L if re.fullmatch(r"\s*[-+]?\d+(\.\d+)?\s*", l)])
    return off, v
def razon(s):
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:]); return r.mean()
def bloques(s, nb=10):
    B = len(s)//nb
    a = np.array([[s[k*B:(k+1)*B].var(), razon(s[k*B:(k+1)*B])] for k in range(nb)])
    return a.std(axis=0, ddof=1)/np.sqrt(nb)
modelos = {
 "A: 0.1748-1.069/L^2": lambda L: 0.1748 - 1.069/L**2,
 "B: 0.1846-0.206/L":   lambda L: 0.1846 - 0.206/L,
 "C: 0.1715-7.292/L^3": lambda L: 0.1715 - 7.292/L**3,
 "D: GUE-0.0193(8.75/L)^1.29": lambda L: GUE_V - 0.0193*(8.75/L)**1.29,
 "E: GUE (constante)":  lambda L: GUE_V + 0*L}
r_mod = lambda L: GUE_R + 0.841/L**2 - 0.0002
print(f"{'origen':>11s} {'ln(t/2pi)':>9s} {'N':>6s} {'var(s)':>16s} {'<r>':>16s}")
filas = []
for nom, f in (("1e12", "zeros3.txt"), ("1e21", "zeros4.txt"), ("1e22", "zeros5.txt")):
    off, v = carga(f)
    s = np.diff(v) * np.log((off + v[:-1]) / (2*np.pi)) / (2*np.pi)
    s = s / s.mean()
    ev, er = bloques(s)
    L = np.log((off + v.mean())/(2*np.pi))
    print(f"{nom:>11s} {L:9.2f} {len(s):6d} {s.var():.4f}+-{ev:.4f} {razon(s):.4f}+-{er:.4f}")
    filas.append((nom, L, s.var(), ev, razon(s), er))
print("\nprediccion de cada modelo de E26 (ajustados con t <= 1e6) y desviacion en sigmas:")
for m, fn in modelos.items():
    cad = []
    for nom, L, vv, ev, rr, er in filas:
        p = fn(L); cad.append(f"{nom}: {p:.4f} ({(vv-p)/ev:+.1f}s)")
    print(f"  {m:30s} " + "  ".join(cad))
print("  <r> (0.5998 + 0.841/L^2): " + "  ".join(f"{nom}: {r_mod(L):.4f} ({(rr-r_mod(L))/er:+.1f}s)" for nom, L, vv, ev, rr, er in filas))
print("\nAviso: una ventana de 1e4 ceros por altura; error de 10 bloques de 1000 espaciados.")
