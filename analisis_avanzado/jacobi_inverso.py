import sys
import numpy as np
import mpmath as mp

# Uso: python3 jacobi_inverso.py [K] [archivo_ceros]   (K=80 tarda ~1 min)
mp.mp.dps = 250
K = int(sys.argv[1]) if len(sys.argv) > 1 else 80
t = np.loadtxt(sys.argv[2] if len(sys.argv) > 2 else "ceros_100k.txt")[:K]


def lanczos(ev):
    n = len(ev)
    ev = [mp.mpf(float(x)) for x in ev]
    q = [mp.mpf(1) / mp.sqrt(n)] * n
    qprev = [mp.mpf(0)] * n
    a, b = [], []
    bb = mp.mpf(0)
    for k in range(n):
        Aq = [ev[i] * q[i] for i in range(n)]
        ak = sum(q[i] * Aq[i] for i in range(n))
        a.append(ak)
        r = [Aq[i] - ak * q[i] - bb * qprev[i] for i in range(n)]
        bb = mp.sqrt(sum(x * x for x in r))
        if k < n - 1:
            b.append(bb)
            qprev = q
            q = [x / bb for x in r]
    return np.array([float(x) for x in a]), np.array([float(x) for x in b])


def suave(K):
    """espectro sin fluctuaciones: N_smooth(E) = theta(E)/pi + 1 = n - 1/2 (biseccion)"""
    out = []
    for n in range(1, K + 1):
        lo, hi = 10.0, float(t[-1]) + 50
        for _ in range(60):
            mid = (lo + hi) / 2
            if float(mp.siegeltheta(mid)) / np.pi + 1 < n - 0.5:
                lo = mid
            else:
                hi = mid
        out.append((lo + hi) / 2)
    return np.array(out)


a, b = lanczos(t)
s = suave(K)
a0, b0 = lanczos(s)
print("n   a_n(ceros)  a_n(suave)   b_n(ceros)  b_n(suave)")
for n in range(0, 10):
    print(f"{n:3d} {a[n]:10.4f} {a0[n]:10.4f} {b[n]:10.4f} {b0[n]:10.4f}")
m = slice(K // 4, 3 * K // 4)
print(f"\nParte central: rms(a-a_suave) = {(a[m]-a0[m]).std():.4f}, rms(b-b_suave) = {(b[m]-b0[m]).std():.4f}")
rng = np.random.default_rng(0)
res = []
for _ in range(20):
    A = rng.normal(size=(2 * K, 2 * K)) + 1j * rng.normal(size=(2 * K, 2 * K))
    E = np.linalg.eigvalsh((A + A.conj().T) / 2)
    x = (E - E.mean()) / E.std()
    u = np.polyval(np.polyfit(x, np.arange(len(E)), 9), x)
    u = u[len(u) // 2 - K // 2:len(u) // 2 + K // 2]
    u = u - u[0]
    gE = np.interp(np.clip(u, 0, K - 1), np.arange(K), s)
    a1, _ = lanczos(gE[:K])
    res.append((a1[m] - a0[m]).std())
print(f"GUE desplegado a la misma densidad: rms(a-a_suave) = {np.mean(res):.4f} +- {np.std(res):.4f}")
