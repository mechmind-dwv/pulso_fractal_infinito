import sys
import numpy as np

# Uso: python3 psi_desde_ceros.py [archivo_ceros]   (por defecto ceros_100k.txt)
t = np.loadtxt(sys.argv[1] if len(sys.argv) > 1 else "ceros_100k.txt")


def psi_exacta(x):
    n = int(x)
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    tot = 0.0
    for p in np.nonzero(s)[0]:
        q = p
        while q <= x:
            tot += np.log(p) * (0.5 if q == x else 1.0)
            q *= p
    return tot


def psi_ceros(x, g):
    lx = np.log(x)
    a = 0.5 * np.cos(g * lx) + g * np.sin(g * lx)
    s = 2 * np.sqrt(x) * np.sum(a / (0.25 + g ** 2))
    return x - s - np.log(2 * np.pi) - 0.5 * np.log(1 - x ** -2.0)


xs = np.linspace(2.0, 60.0, 1500) + 1e-9
ex = np.array([psi_exacta(x) for x in xs])
print("Reconstruccion de psi(x) = suma de ln p sobre p^k <= x, con K ceros (formula explicita de von Mangoldt)")
print(f"{'K ceros':>9s} {'RMS en [2,60]':>14s}")
for K in (10, 100, 1000, 10000, 100000):
    if K > len(t):
        break
    ap = np.array([psi_ceros(x, t[:K]) for x in xs])
    print(f"{K:9d} {np.sqrt(np.mean((ap - ex) ** 2)):14.3f}")
rng = np.random.default_rng(0)
g = np.sort(rng.uniform(t[0], t[999], 1000))
ap = np.array([psi_ceros(x, g) for x in xs])
print(f"control: 1000 gammas uniformes al azar  RMS = {np.sqrt(np.mean((ap - ex) ** 2)):.3f}")
