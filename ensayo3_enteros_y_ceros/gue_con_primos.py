import sys
import numpy as np

# Uso: python3 gue_con_primos.py ARCHIVO_CEROS [T0] [ANCHO] [NREAL]
# GUE desplegado -> niveles t = N_mod^{-1}(u), con N_mod(t) = N_suave(t) + S_P(t) y
# S_P(t) = -(1/pi) sum_{n<=P} w_n Lambda(n)/(sqrt(n) ln n) sin(t ln n)   (serie de primos truncada, peso de Fejer w_n)
ceros = np.loadtxt(sys.argv[1])
T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 40000.0
W = float(sys.argv[3]) if len(sys.argv) > 3 else 400.0
NREAL = int(sys.argv[4]) if len(sys.argv) > 4 else 8
sig = W / 4


def Ns(t):
    th = t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t)
    return th / np.pi + 1


def von_mangoldt(nmax):
    out = np.zeros(nmax + 1)
    for m in range(2, nmax + 1):
        for q in range(2, m + 1):
            if m % q == 0:
                r = m
                while r % q == 0:
                    r //= q
                out[m] = np.log(q) if r == 1 else 0.0
                break
    return out


def razon(z):
    s = np.diff(z)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r.mean()


def S(z, xs):
    w = np.exp(-(z - T0) ** 2 / (2 * sig ** 2))
    return np.array([np.sum(w * np.cos(xx * z)) for xx in xs])


vera = ceros[(ceros >= T0 - 3 * sig) & (ceros <= T0 + 3 * sig)]
nv = len(vera)
tt = np.arange(T0 - 4 * sig, T0 + 4 * sig, 0.002)
Nbase = Ns(tt)
lam = von_mangoldt(300)
I = sig * np.sqrt(2 * np.pi)
rng = np.random.default_rng(1)
US = []                                       # mismas realizaciones GUE para todas las filas (comparacion pareada)
for _ in range(NREAL):
    K = 2 * nv
    A = rng.normal(size=(K, K)) + 1j * rng.normal(size=(K, K))
    E = np.linalg.eigvalsh((A + A.conj().T) / 2)
    x = (E - E.mean()) / E.std()
    u = np.polyval(np.polyfit(x, np.arange(K), 11), x)
    US.append(u[K // 2 - nv // 2:K // 2 + nv // 2])
gen = np.linspace(0.3, 2.6, 41)
gen = gen[np.min(abs(gen[:, None] - np.log(np.arange(2, 14))[None, :]), axis=1) > 0.05]
print(f"Ventana T0={T0:.0f}, sigma={sig:.0f}: {nv} ceros reales. {NREAL} realizaciones GUE por fila.")
print(f"Ceros reales: <r> = {razon(vera):.4f}; picos S(ln n) n=2,3,5,7: " + " ".join(f"{S(vera,[np.log(n)])[0]:7.2f}" for n in (2, 3, 5, 7)) + f";  |S| generico {np.median(abs(S(vera, gen))):.2f}")
print(f"{'P':>4s} {'N_mod mono':>10s} {'<r>':>13s} {'S(2)':>7s} {'S(3)':>7s} {'S(5)':>7s} {'S(7)':>7s} {'|S| gener.':>10s}")
for P in (1, 3, 7, 13, 50, 200):
    nn = np.arange(2, P + 1)
    wn = 1 - np.log(nn) / np.log(P + 1)
    Sp = np.zeros_like(tt)
    for n in nn:
        if lam[n]:
            Sp += -(1 / np.pi) * wn[n - 2] * lam[n] / (np.sqrt(n) * np.log(n)) * np.sin(tt * np.log(n))
    Nmod = Nbase + Sp
    mono = np.mean(np.diff(Nmod) > 0)
    Nmon = np.maximum.accumulate(Nmod)
    rs, peaks, gens = [], [], []
    for u in US:
        z = np.interp(Ns(vera[0]) + (u - u[0]), Nmon, tt)
        rs.append(razon(z))
        peaks.append([S(z, [np.log(n)])[0] for n in (2, 3, 5, 7)])
        gens.append(np.median(abs(S(z, gen))))
    pk = np.mean(peaks, axis=0)
    print(f"{P:4d} {mono:10.3f} {np.mean(rs):7.4f}+-{np.std(rs)/np.sqrt(NREAL):.4f} " + " ".join(f"{v:7.2f}" for v in pk) + f" {np.mean(gens):10.2f}")
print("\nPrediccion de la formula explicita con peso de Fejer, S(ln n) para n=2,3,5,7:")
for P in (3, 7, 13, 50, 200):
    print(f"P={P:4d}: " + " ".join(f"{-lam[n]/(2*np.pi*np.sqrt(n))*I*(1-np.log(n)/np.log(P+1)):7.2f}" if n <= P else "      -" for n in (2, 3, 5, 7)))
