import sys
import numpy as np

RUTA = sys.argv[1] if len(sys.argv) > 1 else "ceros_100k.txt"   # admite .txt y .gz
BLOQUE = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
rng = np.random.default_rng(11)


def theta(t):
    # theta de Riemann-Siegel (serie asintotica, muy precisa para t > 10)
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def n_suave(t):
    return theta(t) / np.pi + 1.0


def sigma2(u, L, ventanas=4000):
    x = rng.uniform(u[0], u[-1] - L, size=ventanas)
    c = np.searchsorted(u, x + L) - np.searchsorted(u, x)
    return c.var()


def stats_bloque(t):
    u = n_suave(t)
    s = np.diff(u)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    c = np.corrcoef(s[:-1], s[1:])[0, 1]
    return r.mean(), 100 * np.mean(s < 0.3), c, sigma2(u, 20.0)


def gue_ref(n=600, mats=20):
    R = np.sqrt(2.0 * n)
    F = lambda x: 0.5 + (x * np.sqrt(1 - x * x) + np.arcsin(x)) / np.pi
    S, U = [], []
    for _ in range(mats):
        A = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
        H = (A + A.conj().T) / 2
        x = np.linalg.eigvalsh(H) / R
        x = x[np.abs(x) < 0.5]
        u = n * F(x)
        U.append(u)
        S.append(np.diff(u))
    s = np.concatenate(S)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    c = np.corrcoef(s[:-1], s[1:])[0, 1]
    s2 = np.mean([sigma2(u, 20.0, 800) for u in U])
    return r.mean(), 100 * np.mean(s < 0.3), c, s2


def main():
    z = np.loadtxt(RUTA)
    t = z[:, -1] if z.ndim > 1 else z
    N = len(t)
    n = np.arange(1, N + 1)
    print(f"{N} ceros, t = {t[0]:.2f} .. {t[-1]:.2f}\n")

    # ---------------- 1. Chequeo tipo Turing (conteo) ----------------
    # N(T) = theta(T)/pi + 1 + S(T). En el cero n-esimo, N = n - 1/2 (a mitad del salto).
    d = (n - 0.5) - n_suave(t)             # ~ S(t_n); debe oscilar alrededor de 0, acotado
    sd_teo = np.sqrt(np.log(np.log(t[-1])) / (2 * np.pi ** 2))
    print("=== 1. Conteo: d_n = (n - 1/2) - N_suave(t_n)  [~ S(T)] ===")
    print(f"media={d.mean():+.3f}  desv={d.std():.3f}  max|d|={np.abs(d).max():.2f}  "
          f"(desv. teorica ~ {sd_teo:.2f} a T final)")
    print("media de d por bloque:")
    for k in range(0, N, BLOQUE):
        b = d[k:k + BLOQUE]
        print(f"  ceros {k+1:>9d}..{k+len(b):>9d}   t~{t[k]:>10.1f}   media d = {b.mean():+.4f}")

    # control de sensibilidad: borrar UN cero y ver si el detector lo nota
    k = N // 2
    t2 = np.delete(t, k)
    d2 = (np.arange(1, N) - 0.5) - n_suave(t2)
    antes = d[max(0, k - 5000):k].mean()
    despues = d2[k:k + 5000].mean()
    print(f"\nControl: borro el cero n={k+1}. Salto de la media local de d: "
          f"{despues - antes:+.2f}  (esperado ~ -1)")
    print("Un cero perdido o 'escapado' se veria como un escalon de -1 en d. "
          "Sin escalones: el conteo es consistente.")
    print("OJO: esto es una comprobacion numerica, no la version rigurosa de Turing "
          "(Platt-Trudgian usan aritmetica de bolas y cotas de S(T)).")

    # ---------------- 2. Estadistica por altura ----------------
    print("\n=== 2. Estadistica por bloques de altura (vs GUE) ===")
    print("GUE referencia (matrices simuladas)...")
    g = gue_ref()
    print(f"{'bloque':>20s} {'t medio':>11s} {'<r>':>7s} {'P(s<.3)%':>9s} {'corr':>7s} {'Sigma2(20)':>11s}")
    for k in range(0, N - BLOQUE + 1, BLOQUE):
        tb = t[k:k + BLOQUE]
        r, p, c, s2 = stats_bloque(tb)
        print(f"{k+1:>9d}..{k+BLOQUE:<9d} {tb.mean():>11.1f} {r:7.4f} {p:9.2f} {c:+7.3f} {s2:11.3f}")
    print(f"{'GUE sim':>20s} {'':>11s} {g[0]:7.4f} {g[1]:9.2f} {g[2]:+7.3f} {g[3]:11.3f}")
    print("GUE asintotico: <r>=0.5996, Sigma2(20)=0.650. "
          "Se espera que los bloques se acerquen a GUE al subir t (Bogomolny-Keating).")


if __name__ == "__main__":
    main()
