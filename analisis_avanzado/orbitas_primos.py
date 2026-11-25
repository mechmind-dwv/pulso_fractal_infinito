import sys
import numpy as np

RUTA = sys.argv[1] if len(sys.argv) > 1 else "ceros_2M.txt.gz"   # admite .txt y .gz
rng = np.random.default_rng(3)


def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def von_mangoldt(n):
    # Lambda(n) = ln p si n = p^m, 0 en otro caso
    if n < 2:
        return 0.0
    m, p = n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            return float(np.log(p)) if m == 1 else 0.0
        p += 1
    return float(np.log(n))          # n primo


def S(t, w, x, bloque=200000):
    # S(x) = sum_n w(t_n) cos(x t_n)
    tot = 0.0
    for i in range(0, len(t), bloque):
        tt = t[i:i + bloque]
        tot += np.sum(w[i:i + bloque] * np.cos(x * tt))
    return tot


def main():
    t = np.loadtxt(RUTA)
    t = t[:, -1] if t.ndim > 1 else t
    N = len(t)
    print(f"{N} ceros, t = {t[0]:.1f} .. {t[-1]:.1f}\n")

    # ===== 1. Formula explicita: los ceros contienen los primos =====
    Tc = t[-1] / 4.0                                   # peso gaussiano ~ 3e-4 en el ultimo cero
    w = np.exp(-0.5 * (t / Tc) ** 2)
    I = Tc * np.sqrt(np.pi / 2)                        # int_0^inf w dt
    print("=== 1. Formula explicita (Weil): S(x) = sum w(t_n) cos(x t_n) ===")
    print(f"peso gaussiano Tc = {Tc:.0f};  prediccion en x = ln n :  S = -Lambda(n)/(2 pi sqrt(n)) * {I:.0f}\n")
    print(f"{'n':>4s} {'tipo':>10s} {'S medido':>11s} {'S predicho':>11s} {'cociente':>9s}")
    for n in range(2, 41):
        lam = von_mangoldt(n)
        pr = -lam / (2 * np.pi * np.sqrt(n)) * I
        s = S(t, w, np.log(n))
        tipo = "pot. primo" if lam > 0 else "otro"
        cq = f"{s / pr:9.4f}" if lam > 0 else "        -"
        print(f"{n:4d} {tipo:>10s} {s:11.1f} {pr:11.1f} {cq}")
    # control: x genericos (no son ln de entero)
    xs = rng.uniform(0.6, 3.7, size=40)
    xs = np.array([x for x in xs if np.min(np.abs(np.exp(x) - np.round(np.exp(x)))) > 1e-3 * np.exp(x)])
    ctl = np.array([S(t, w, x) for x in xs])
    print(f"\nControl, {len(xs)} valores de x genericos: |S| tipico = {np.abs(ctl).mean():.1f}, max = {np.abs(ctl).max():.1f}  "
          f"(frente a picos de decenas de miles)")

    # ===== 2. Factor de forma espectral K(tau): rampa =====
    u = theta(t) / np.pi + 1.0
    Nw = 5000
    nwin = N // Nw
    taus = np.arange(0.05, 2.0001, 0.05)
    K = np.zeros(len(taus))
    Ksh = np.zeros(len(taus))
    Kpo = np.zeros(len(taus))
    for k in range(nwin):
        uu = u[k * Nw:(k + 1) * Nw]
        uu = uu - uu[0]
        d = np.diff(uu)
        sh = np.concatenate([[0.0], np.cumsum(rng.permutation(d))])
        po = np.concatenate([[0.0], np.cumsum(rng.exponential(size=len(d)) * d.mean())])
        for arr, acc in ((uu, K), (sh, Ksh), (po, Kpo)):
            ph = np.exp(2j * np.pi * np.outer(taus, arr))
            acc += np.abs(ph.sum(axis=1)) ** 2 / len(arr)
    K /= nwin; Ksh /= nwin; Kpo /= nwin
    print(f"\n=== 2. Factor de forma espectral K(tau), {nwin} ventanas de {Nw} ceros ===")
    print("GUE (rampa y meseta): K = tau si tau<1, K = 1 si tau>1. Poisson: K = 1.\n")
    print(f"{'tau':>5s} {'ceros':>8s} {'GUE':>6s} {'barajado':>9s} {'Poisson':>8s}")
    for tau, a, b, c in zip(taus, K, Ksh, Kpo):
        print(f"{tau:5.2f} {a:8.3f} {min(tau, 1.0):6.2f} {b:9.3f} {c:8.3f}")


if __name__ == "__main__":
    main()
