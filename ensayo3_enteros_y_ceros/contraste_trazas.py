import sys
import numpy as np

# Uso: python3 contraste_trazas.py ARCHIVO_CEROS [T0] [ANCHO]
# Los primos (picos en x = ln n) en la suma de cosenos de los ceros, frente a un espectro GUE de igual densidad.
ceros = np.loadtxt(sys.argv[1])
T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 40000.0
W = float(sys.argv[3]) if len(sys.argv) > 3 else 400.0
sig = W / 4


def Ns(t):
    th = t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t)
    return th / np.pi + 1


def von_mangoldt(nmax):
    out = {}
    for m in range(2, nmax + 1):
        for q in range(2, m + 1):
            if m % q == 0:
                r = m
                while r % q == 0:
                    r //= q
                out[m] = np.log(q) if r == 1 else 0.0
                break
    return out


lam = von_mangoldt(60)
vera = ceros[(ceros >= T0 - 3 * sig) & (ceros <= T0 + 3 * sig)]
rng = np.random.default_rng(0)
K = 2 * len(vera)
A = rng.normal(size=(K, K)) + 1j * rng.normal(size=(K, K))
E = np.linalg.eigvalsh((A + A.conj().T) / 2)
x = (E - E.mean()) / E.std()
u = np.polyval(np.polyfit(x, np.arange(K), 11), x)
u = u[K // 2 - len(vera) // 2:K // 2 + len(vera) // 2]
tt = np.linspace(T0 - 4 * sig, T0 + 4 * sig, 200001)
sur = np.interp(Ns(vera[0]) + (u - u[0]), Ns(tt), tt)


def S(z, xs):
    w = np.exp(-(z - T0) ** 2 / (2 * sig ** 2))
    return np.array([np.sum(w * np.cos(xx * z)) for xx in xs])


I = sig * np.sqrt(2 * np.pi)
print(f"Ventana centrada en T0={T0:.0f}, sigma={sig:.0f}: {len(vera)} ceros reales y {len(sur)} niveles GUE de igual densidad")
print(f"{'n':>3s} {'x=ln n':>7s} {'pred. Weil':>11s} {'ceros':>9s} {'GUE':>9s}")
for n in (2, 3, 4, 5, 7, 8, 9, 11, 13):
    xx = np.log(n)
    pred = -lam[n] / (2 * np.pi * np.sqrt(n)) * I if lam[n] else 0.0
    print(f"{n:3d} {xx:7.4f} {pred:11.2f} {S(vera, [xx])[0]:9.2f} {S(sur, [xx])[0]:9.2f}")
gen = np.linspace(0.3, 2.6, 41)
gen = gen[np.min(abs(gen[:, None] - np.log(np.arange(2, 14))[None, :]), axis=1) > 0.05]
print(f"\nControl, {len(gen)} x genericos: |S| tipico ceros {np.median(abs(S(vera, gen))):.2f} max {abs(S(vera, gen)).max():.2f}"
      f"  |  GUE {np.median(abs(S(sur, gen))):.2f} max {abs(S(sur, gen)).max():.2f}")
