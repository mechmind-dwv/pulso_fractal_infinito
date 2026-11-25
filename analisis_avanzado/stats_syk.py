import sys
import itertools
import numpy as np
import syk                                   # usa N=argv[1], NREAL=argv[2], semilla=argv[3]

G = np.euler_gamma


def delta3(u, Ls, m=300, rng=None):
    """Delta3(L): minimo cuadratico de la escalera N(E) frente a una recta en [a,a+L], promediado sobre a."""
    rng = rng or np.random.default_rng(0)
    out = []
    for L in Ls:
        v = []
        xs = np.linspace(0.0, 1.0, 200)
        for a in rng.uniform(u[0], u[-1] - L, m):
            e = u[(u >= a) & (u <= a + L)] - a
            N = np.searchsorted(e, xs * L).astype(float)
            p = np.polyfit(xs, N, 1)
            v.append(np.mean((N - np.polyval(p, xs)) ** 2))
        out.append(np.mean(v))
    return np.array(out)


def d3_teo(L, clase):
    if clase == "GUE":
        return (np.log(2 * np.pi * L) + G - 5 / 4) / (2 * np.pi ** 2)
    if clase == "GOE":
        return (np.log(2 * np.pi * L) + G - 5 / 4 - np.pi ** 2 / 8) / np.pi ** 2
    return np.nan                              # GSE: no incluida (no verificada)


def wigner(s, clase):
    if clase == "GOE":
        return np.pi / 2 * s * np.exp(-np.pi * s ** 2 / 4)
    if clase == "GUE":
        return 32 / np.pi ** 2 * s ** 2 * np.exp(-4 * s ** 2 / np.pi)
    return 2 ** 18 / (3 ** 6 * np.pi ** 3) * s ** 4 * np.exp(-64 * s ** 2 / (9 * np.pi))


def main():
    N, NREAL = syk.NF, syk.NREAL
    clase = {2: "GUE", 6: "GUE", 0: "GOE", 4: "GSE"}[N % 8]
    D, s0, par, zsign = syk.preparar(N)
    quads = list(itertools.combinations(range(N), 4))
    Ls = [2, 3, 5, 8, 12, 18]
    S, d3 = [], []
    for r in range(NREAL):
        H = syk.hamiltoniano(N, quads, D, s0, par, zsign)
        E = syk.niveles(np.linalg.eigvalsh(H))
        u = syk.desplegar(E)
        S.append(np.diff(u) / np.mean(np.diff(u)))
        d3.append(delta3(u, [L for L in Ls if L < (u[-1] - u[0]) / 2]))
        print(f"  realizacion {r+1}/{NREAL}", flush=True)
    S = np.concatenate(S)
    print(f"\nN={N}, clase esperada {clase}, {len(S)} espaciados")
    print("\nP(s): histograma vs surmise de Wigner (aproximado, +-2% intrinseco)")
    print(f"{'s':>5s} {'medido':>8s} {'GOE':>7s} {'GUE':>7s} {'GSE':>7s}")
    h, ed = np.histogram(S, bins=np.arange(0, 3.01, 0.25), density=True)
    for i, hv in enumerate(h):
        c = (ed[i] + ed[i + 1]) / 2
        print(f"{c:5.2f} {hv:8.3f} {wigner(c,'GOE'):7.3f} {wigner(c,'GUE'):7.3f} {wigner(c,'GSE'):7.3f}")
    print("\nDelta3(L), media de realizaciones")
    n = min(len(x) for x in d3)
    M = np.mean([x[:n] for x in d3], axis=0)
    print(f"{'L':>4s} {'medido':>8s} {'GUE':>7s} {'GOE':>7s}")
    for L, mv in zip(Ls[:n], M):
        print(f"{L:4d} {mv:8.4f} {d3_teo(L,'GUE'):7.4f} {d3_teo(L,'GOE'):7.4f}")


if __name__ == "__main__":
    main()
