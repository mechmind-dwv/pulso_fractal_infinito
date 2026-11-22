import sys
import numpy as np

RUTA = sys.argv[1] if len(sys.argv) > 1 else "ceros_100k.txt"   # admite .txt y .gz


def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def main():
    z = np.loadtxt(RUTA)
    t = z[:, -1] if z.ndim > 1 else z
    n = np.arange(1, len(t) + 1)
    d = (n - 0.5) - (theta(t) / np.pi + 1.0)          # ~ S(t_n)
    print(f"{len(t)} ceros, t = {t[0]:.1f} .. {t[-1]:.1f}")
    print("Selberg: Var S(T) ~ (1/(2 pi^2)) ln ln T + cte  -> pendiente teorica = %.4f\n"
          % (1 / (2 * np.pi ** 2)))

    # bloques logaritmicos en t (desde t=200 para evitar el regimen de pocos ceros)
    edges = np.geomspace(200, t[-1], 9)
    x, v, k = [], [], []
    print(f"{'t rango':>22s} {'n':>8s} {'lnln T':>7s} {'var S':>7s} {'curtosis':>9s}")
    for a, b in zip(edges[:-1], edges[1:]):
        m = (t >= a) & (t < b)
        if m.sum() < 500:
            continue
        s = d[m]
        s = s - s.mean()
        var = s.var()
        kur = np.mean(s ** 4) / var ** 2          # Gaussiana = 3
        tm = np.exp(np.mean(np.log(t[m])))
        x.append(np.log(np.log(tm)))
        v.append(var)
        k.append(kur)
        print(f"{a:>10.0f}..{b:<10.0f} {m.sum():8d} {x[-1]:7.3f} {var:7.4f} {kur:9.3f}")

    x, v = np.array(x), np.array(v)
    pend, ordenada = np.polyfit(x, v, 1)
    print(f"\nAjuste var S = a + b*lnln T:  b = {pend:.4f}  (teoria {1/(2*np.pi**2):.4f}),  a = {ordenada:+.3f}")
    print("Curtosis ~ 3 => S(T) casi gaussiana (teorema central de Selberg, convergencia muy lenta).")
    print("OJO: ln ln T cambia poco en este rango; la pendiente tiene incertidumbre grande.")


if __name__ == "__main__":
    main()
