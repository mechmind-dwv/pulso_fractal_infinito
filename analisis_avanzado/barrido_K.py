import sys, time
import numpy as np

# Uso: python3 barrido_K.py ARCHIVO_CEROS [XMAX] [KMAX]
# RMS(D - reconstruccion con K ceros) segun K, con D = (psi(x)-x)/sqrt(x), x de 20 a XMAX.
# KMAX por defecto 100000 (unos 30 s en un PC; 10^6 son ~5 min).
XMAX = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10_000_000
KMAX = int(float(sys.argv[3])) if len(sys.argv) > 3 else 100_000
ceros = np.loadtxt(sys.argv[1])

lam = np.zeros(XMAX + 1)
cr = np.ones(XMAX + 1, bool)
cr[:2] = False
for i in range(2, int(XMAX ** 0.5) + 1):
    if cr[i]:
        cr[i * i::i] = False
for p in np.nonzero(cr)[0]:
    q = p
    while q <= XMAX:
        lam[q] = np.log(p)
        q *= p
psi = np.cumsum(lam)

x = np.exp(np.linspace(np.log(20), np.log(XMAX), 6000))
lx = np.log(x)
D = (psi[np.floor(x).astype(int)] - x) / np.sqrt(x)
acum = -(np.log(2 * np.pi) + 0.5 * np.log(1 - x ** -2.0)) / np.sqrt(x)
S = np.zeros_like(x)
k0 = 0
print(f"x hasta {XMAX:.0e}; ceros disponibles {len(ceros)}; RMS(D) = {np.sqrt(np.mean(D ** 2)):.3f}")
print("      K    gamma_max   RMS(D - recon)   tiempo(s)")
for K in [100, 1000, 10000, 100000, 1000000, 2000000]:
    if K > len(ceros) or K > KMAX:
        break
    t0 = time.time()
    g = ceros[k0:K]
    for i in range(0, len(g), 1000):
        gg = g[i:i + 1000][None, :]
        a = (0.5 * np.cos(gg * lx[:, None]) + gg * np.sin(gg * lx[:, None])) / (0.25 + gg ** 2)
        S += -2 * a.sum(axis=1)
    k0 = K
    print(f"{K:8d} {g[-1]:11.0f} {np.sqrt(np.mean((D - S - acum) ** 2)):14.3f} {time.time() - t0:11.0f}", flush=True)
