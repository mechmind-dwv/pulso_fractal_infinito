import sys
import numpy as np
z = np.loadtxt(sys.argv[1]); B = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
VG, RG = 0.1800, 0.5996
def Nsm(t): return t/(2*np.pi)*np.log(t/(2*np.pi)) - t/(2*np.pi) + 7/8 + 1/(48*np.pi*t)
def razon(s):
    r = np.minimum(s[:-1], s[1:])/np.maximum(s[:-1], s[1:]); return r.mean()
T, V, R = [], [], []
for i in range(0, len(z) - B + 1, B):
    w = z[i:i+B]; s = np.diff(Nsm(w)); s = s/s.mean()
    T.append(w.mean()); V.append(s.var()); R.append(razon(s))
T, V, R = map(np.array, (T, V, R)); L = np.log(T/(2*np.pi))
print(f"{len(T)} bloques de {B} ceros")
print("   t medio   ln(t/2pi)   var(s)   GUE-var   <r>    <r>-GUE")
for t, l, v, r in zip(T, L, V, R):
    print(f"{t:10.4g}   {l:7.2f}   {v:.4f}   {VG-v:+.4f}   {r:.4f}  {r-RG:+.4f}")
for nombre, y in (("exceso de var (GUE - var)", VG - V), ("exceso de <r> (r - GUE)", R - RG)):
    ok = y > 0
    if ok.sum() >= 4:
        a, b = np.polyfit(np.log(L[ok]), np.log(y[ok]), 1)
        print(f"ajuste {nombre}: y ~ ln^({a:.2f})  (con {ok.sum()} bloques; heuristico)")
