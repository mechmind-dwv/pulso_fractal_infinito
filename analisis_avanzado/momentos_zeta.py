import sys
import numpy as np
import mpmath as mp

mp.mp.dps = 30
rng = np.random.default_rng(1)


# ---------------- Z(t) de Riemann-Siegel, vectorizada (suma principal + C0) ----------------
def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def Z_rs(t):
    t = np.asarray(t, float)
    a = np.sqrt(t / (2 * np.pi))
    N = np.floor(a).astype(int)
    p = a - N
    th = theta(t)
    s = np.zeros_like(t)
    for n in range(1, N.max() + 1):
        m = n <= N
        s += np.where(m, np.cos(th - t * np.log(n)) / np.sqrt(n), 0.0)
    s *= 2
    c = np.cos(2 * np.pi * p)
    p2 = np.where(np.abs(c) < 1e-6, p + 1e-5, p)     # singularidad evitable de C0
    C0 = np.cos(2 * np.pi * (p2 ** 2 - p2 - 1 / 16)) / np.cos(2 * np.pi * p2)
    return s + (-1.0) ** (N - 1) * (2 * np.pi / t) ** 0.25 * C0


# ---------------- factor aritmetico a(k) de Keating-Snaith ----------------
def primos(n):
    c = np.ones(n + 1, bool)
    c[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if c[i]:
            c[i * i::i] = False
    return np.nonzero(c)[0]


def a_k(k, P=10 ** 6, M=60):
    x = 1.0 / primos(P).astype(float)
    S = np.ones_like(x)
    c = 1.0
    xm = np.ones_like(x)
    for m in range(1, M + 1):
        c *= ((m - 1 + k) / m) ** 2
        xm = xm * x
        S += c * xm
    return float(np.exp(np.sum(k ** 2 * np.log1p(-x) + np.log(S))))


def Mcue(N, k):
    # G(N+1) G(N+2k+1) G(k+1)^2 / ( G(2k+1) G(N+k+1)^2 ),  N real (continuacion con Barnes G)
    G = mp.barnesg
    return float(G(N + 1) * G(N + 2 * k + 1) * G(k + 1) ** 2 / (G(2 * k + 1) * G(N + k + 1) ** 2))


def lead(k):
    G = mp.barnesg
    return float(G(k + 1) ** 2 / G(2 * k + 1))


def pred(t, k, ak):
    L = np.log(t / (2 * np.pi))
    return ak * Mcue(L, k)


# ---------------- main ----------------
def validar():
    print("== Validacion de Z_rs contra mpmath.siegelz ==")
    for t in (1000.3, 54321.7, 250000.123, 1.0e7 + 0.37):
        zz = float(mp.siegelz(t))
        print(f"  t={t:>11.3f}  RS={Z_rs(np.array([t]))[0]:+.8f}  mpmath={zz:+.8f}  dif={abs(Z_rs(np.array([t]))[0]-zz):.1e}")
    print("== a(k): a(1)=1, a(2)=6/pi^2=%.6f ==" % (6 / np.pi ** 2))
    print("  a(1)=%.6f  a(2)=%.6f" % (a_k(1.0), a_k(2.0)))


def main(Tmax=1e6):
    Ts = [x for x in (1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8) if x <= Tmax]
    ks = [0.5, 1.0, 1.5, 2.0]
    AK = {k: a_k(k) for k in ks}
    print("a(k):", {k: round(v, 5) for k, v in AK.items()})
    print("\nVentanas [T, T+W], W=min(20000,T), paso 0.05. Momento M_k = media de |zeta|^{2k}.\n")
    res = {k: [] for k in ks}
    Ls = []
    varlog = []
    for T in Ts:
        W = min(20000.0, T)
        t = np.arange(T, T + W, 0.05)
        z = np.abs(Z_rs(t))
        L = float(np.mean(np.log(t / (2 * np.pi))))
        Ls.append(L)
        lz = np.log(np.maximum(z, 1e-300))
        varlog.append((L, lz.mean(), lz.var()))
        nb = 20
        for k in ks:
            f = z ** (2 * k)
            M = f.mean()
            blk = np.array([b.mean() for b in np.array_split(f, nb)])
            se = blk.std(ddof=1) / np.sqrt(nb)
            # prediccion promediada en la ventana
            tt = np.linspace(T, T + W, 12)
            P = np.mean([pred(x, k, AK[k]) for x in tt])
            Lead = AK[k] * lead(k) * np.mean(np.log(t / (2 * np.pi))) ** (k * k)
            res[k].append((L, M, se, P, Lead))
    for k in ks:
        print(f"--- k={k}  (exponente teorico k^2 = {k*k:g}; coef. principal a*G = {AK[k]*lead(k):.5f}) ---")
        print(f"{'T':>8s} {'L=ln(t/2pi)':>11s} {'M_k medido':>11s} {'+-':>8s} {'KS(N=L)':>10s} {'medido/KS':>10s} {'a*G*L^k2':>10s} {'medido/lead':>11s}")
        for T, (L, M, se, P, Lead) in zip(Ts, res[k]):
            print(f"{T:8.0f} {L:11.3f} {M:11.4f} {se:8.4f} {P:10.4f} {M/P:10.3f} {Lead:10.4f} {M/Lead:11.3f}")
        Lv = np.array([r[0] for r in res[k]])
        Mv = np.array([r[1] for r in res[k]])
        Pv = np.array([r[3] for r in res[k]])
        sd, _ = np.polyfit(np.log(Lv), np.log(Mv), 1)
        sp, _ = np.polyfit(np.log(Lv), np.log(Pv), 1)
        print(f"  exponente efectivo d ln M / d ln L:  datos = {sd:.2f}   prediccion KS en el mismo rango = {sp:.2f}   asintotico = {k*k:g}\n")
    print("Chequeo k=1 (Hardy-Littlewood): la media LOCAL es M_1 = L + 2*gamma  (la acumulada es L + 2*gamma - 1)")
    for T, (L, M, se, P, Lead) in zip(Ts, res[1.0]):
        print(f"  T={T:10.0f}  medido={M:9.4f} +- {se:.4f}   exacto={L + 2*np.euler_gamma:9.4f}   dif={M - L - 2*np.euler_gamma:+.4f}")
    print("\nln|zeta| (Selberg): media ~ 0 y Var ~ (1/2) ln ln T + cte")
    print(f"{'T':>10s} {'ln ln T':>8s} {'media':>8s} {'Var ln|z|':>10s}")
    xs, vs = [], []
    for T, (L, m, v) in zip(Ts, varlog):
        x = np.log(np.log(T))
        xs.append(x); vs.append(v)
        print(f"{T:10.0f} {x:8.3f} {m:+8.4f} {v:10.4f}")
    b, a0 = np.polyfit(xs, vs, 1)
    print(f"  ajuste Var = a + b*lnlnT:  b = {b:.3f}  (teoria 0.500),  a = {a0:+.3f}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "validar":
        validar()
    else:
        main(float(sys.argv[1]) if len(sys.argv) > 1 else 1e6)
