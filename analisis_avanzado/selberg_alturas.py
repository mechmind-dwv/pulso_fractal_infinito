import re
import sys
import numpy as np
import mpmath as mp

mp.mp.dps = 60
rng = np.random.default_rng(5)

# Ancla: ultimo bloque de los 2.001.052 ceros (t ~ 4e5..1.1e6)
X0, V0 = 2.603, 0.1125
PEND_TEO, PEND_OBS = 1 / (2 * np.pi ** 2), 0.0608


def theta(t):
    return t / 2 * mp.log(t / (2 * mp.pi)) - t / 2 - mp.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def leer(ruta):
    lineas = open(ruta, encoding="utf-8", errors="ignore").read().splitlines()
    cab = " ".join(l.strip() for l in lineas[:12])
    m_off = re.search(r"gamma\s*-\s*(\d+)", cab)
    m_exp = re.search(r"10\^(\d+)\s*\+\s*1\b", cab)
    if not (m_off and m_exp):
        sys.exit(f"{ruta}: no pude leer offset/indice de la cabecera:\n{cab}")
    offset = mp.mpf(m_off.group(1))
    n0 = mp.mpf(10) ** int(m_exp.group(1)) + 1
    vals = [l.strip() for l in lineas if re.fullmatch(r"\s*[-+]?\d+(\.\d+)?\s*", l)]
    return offset, n0, vals, int(m_exp.group(1))


def analizar(ruta):
    offset, n0, vals, ex = leer(ruta)
    N = len(vals)
    d = np.empty(N)
    tf = np.empty(N)
    for i, v in enumerate(vals):
        t = offset + mp.mpf(v)
        n = n0 + i
        d[i] = float((n - mp.mpf(1) / 2) - (theta(t) / mp.pi + 1))
        tf[i] = float(t)
    s = d - d.mean()
    var = s.var()
    kur = np.mean(s ** 4) / var ** 2
    # error estandar: bootstrap de bloques moviles (los S(t_n) vecinos estan muy correlacionados)
    B = 500
    nb = N // B
    bloques = [s[k * B:(k + 1) * B] for k in range(nb)]
    boots = []
    for _ in range(300):
        idx = rng.integers(0, nb, size=nb)
        x = np.concatenate([bloques[j] for j in idx])
        boots.append(x.var())
    se = np.std(boots)
    T = float(np.exp(np.mean(np.log(tf))))
    x = np.log(np.log(T))
    return dict(f=ruta, ex=ex, N=N, T=T, x=x, mean=d.mean(), mx=np.abs(d).max(),
                var=var, se=se, kur=kur)


def main():
    rutas = sys.argv[1:]
    if not rutas:
        sys.exit("uso: python3 selberg_alturas.py zeros3.txt [zeros4.txt zeros5.txt]")
    print("Var S(T) a gran altura (Odlyzko). Prediccion anclada en t~10^6 (var=%.4f, lnlnT=%.3f)\n" % (V0, X0))
    print(f"{'archivo':>12s} {'idx':>6s} {'N':>6s} {'lnln T':>7s} {'media d':>8s} {'max|d|':>7s} "
          f"{'var S':>7s} {'+-':>6s} {'curt':>5s} {'pred b=.0507':>13s} {'pred b=.0608':>13s}")
    for r in rutas:
        a = analizar(r)
        p1 = V0 + PEND_TEO * (a['x'] - X0)
        p2 = V0 + PEND_OBS * (a['x'] - X0)
        print(f"{r:>12s} 1e{a['ex']:<4d} {a['N']:6d} {a['x']:7.3f} {a['mean']:+8.3f} {a['mx']:7.2f} "
              f"{a['var']:7.4f} {a['se']:6.4f} {a['kur']:5.2f} {p1:13.4f} {p2:13.4f}")
    print("\nSi |media d| o max|d| son grandes (>~1.5) revisa el indice/offset: un error de indice desplaza d en un entero.")
    print("Cada ventana es UNA realizacion de 10^4 ceros: la incertidumbre (+-) es grande; no sobreinterpretar.")


if __name__ == "__main__":
    main()
