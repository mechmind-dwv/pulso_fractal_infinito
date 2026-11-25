import sys
import time
import itertools
import numpy as np

NF = int(sys.argv[1]) if len(sys.argv) > 1 else 22          # fermiones de Majorana (par); N mod 8 = 2,6 -> GUE
NREAL = int(sys.argv[2]) if len(sys.argv) > 2 else 8        # realizaciones del desorden
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 12)


def preparar(N):
    D = 1 << (N // 2)
    s0 = np.arange(D, dtype=np.int64)
    par = np.zeros(D, np.int64)                     # paridad del numero de unos de cada entero < D
    x = s0.copy()
    while x.any():
        par ^= x & 1
        x >>= 1
    zsign = np.where(par == 1, -1.0, 1.0)
    return D, s0, par, zsign


def termino(quad, s0, zsign):
    """chi_i chi_j chi_k chi_l |s> = fase(s) |s xor mask>  (Jordan-Wigner: chi_{2a}=Z..Z X_a/sqrt2, chi_{2a+1}=Z..Z Y_a/sqrt2)"""
    cur = np.zeros_like(s0)                         # mascara acumulada (los giros no dependen del estado)
    ph = np.ones(s0.size, np.complex128)
    for m in reversed(quad):
        a, tipo = divmod(m, 2)
        s = s0 ^ cur
        ph *= zsign[s & ((1 << a) - 1)]
        if tipo == 1:                               # Y|b> = i(-1)^b |b xor 1>
            ph *= 1j * (1.0 - 2.0 * ((s >> a) & 1))
        cur ^= (1 << a)
    return int(cur[0]), ph / 4.0                    # (1/sqrt2)^4


def hamiltoniano(N, quads, D, s0, par, zsign):
    sigma = np.sqrt(6.0 / N ** 3)                   # <J^2> = 3! J^2 / N^3, J = 1
    J = rng.normal(0.0, sigma, size=len(quads))
    H = np.zeros((D, D), np.complex128)
    for q, quad in enumerate(quads):
        mask, ph = termino(quad, s0, zsign)
        H[s0 ^ mask, s0] += J[q] * ph
    idx = np.nonzero(par == 0)[0]                   # sector de paridad par: [H, (-1)^F] = 0
    return H[np.ix_(idx, idx)]


def niveles(E):
    d = np.diff(E)
    if np.mean(d < 1e-9) > 0.3:                     # degeneracion de Kramers (clase GSE)
        E = E[::2]
    return E


def desplegar(E, frac=0.6, grado=9):
    n = len(E)
    lo, hi = int(n * (1 - frac) / 2), int(n * (1 + frac) / 2)
    e = E[lo:hi]
    c = np.arange(lo, hi, dtype=float)
    x = (e - e.mean()) / e.std()
    return np.polyval(np.polyfit(x, c, grado), x)


def razon(E):
    s = np.diff(E)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r.mean()


def K_teo(tau, clase):
    tau = np.asarray(tau, float)
    if clase == "GUE":
        return np.minimum(tau, 1.0)
    if clase == "GOE":
        a = 2 * tau - tau * np.log(1 + 2 * tau)
        b = 2 - tau * np.log((2 * tau + 1) / np.abs(2 * tau - 1) + 1e-300)
        return np.where(tau < 1, a, b)
    return np.where(tau < 2, tau / 2 - tau / 4 * np.log(np.abs(1 - tau) + 1e-300), 1.0)


def main():
    N = NF
    clase = {2: "GUE", 6: "GUE", 0: "GOE", 4: "GSE"}[N % 8]
    t0 = time.time()
    print(f"SYK q=4, N={N} Majorana (N mod 8 = {N % 8} -> clase esperada {clase}), dim H = 2^{N//2}, "
          f"sector de paridad par = {1 << (N // 2 - 1)} estados")
    D, s0, par, zsign = preparar(N)
    quads = list(itertools.combinations(range(N), 4))
    print(f"terminos de 4 Majorana: {len(quads)}\n")
    taus = np.arange(0.05, 2.0001, 0.1)             # evita tau=1 (divergencia logaritmica de GSE)
    Kacc = np.zeros(len(taus))
    rs, nl = [], 0
    for r in range(NREAL):
        H = hamiltoniano(N, quads, D, s0, par, zsign)
        if r == 0:
            print(f"  comprobacion: max|H - H^dagger| = {np.abs(H - H.conj().T).max():.1e}")
        E = niveles(np.linalg.eigvalsh(H))
        rs.append(razon(E[len(E) // 5: 4 * len(E) // 5]))
        u = desplegar(E)
        Kacc += np.abs(np.exp(2j * np.pi * np.outer(taus, u)).sum(axis=1)) ** 2 / len(u)
        nl = len(u)
        print(f"  realizacion {r+1}/{NREAL}: {len(E)} niveles, <r> = {rs[-1]:.4f}   ({time.time()-t0:.0f} s)", flush=True)
    K = Kacc / NREAL
    rm = np.mean(rs)
    rse = np.std(rs, ddof=1) / np.sqrt(len(rs)) if len(rs) > 1 else float("nan")
    print(f"\n<r> medio = {rm:.4f} +- {rse:.4f}   (GOE 0.5307, GUE 0.5996, GSE 0.6744;  ceros de Riemann a t~1e6: 0.6056)")
    print(f"\nFactor de forma K(tau), {NREAL} realizaciones, {nl} niveles centrales cada una (ruidoso: +-0.05..0.1 por punto)")
    print(f"{'tau':>5s} {'SYK':>8s} {'GUE':>6s} {'GOE':>6s} {'GSE':>6s}")
    for tau, k, a, b, c in zip(taus, K, K_teo(taus, "GUE"), K_teo(taus, "GOE"), K_teo(taus, "GSE")):
        print(f"{tau:5.2f} {k:8.3f} {a:6.2f} {b:6.2f} {c:6.2f}")


if __name__ == "__main__":
    main()
