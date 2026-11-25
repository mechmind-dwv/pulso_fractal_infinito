import sys
import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
a = 1.0
b = float(sys.argv[2]) if len(sys.argv) > 2 else 200.0


def razon(E):
    s = np.diff(E)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r.mean()


# A) xp exacto en u = ln x: H = -i d/du, caja de longitud ln(b/a)
Lu = np.log(b / a)
n = np.arange(1, 400)
EA = 2 * np.pi * n / Lu
print(f"A) xp en caja log: niveles equiespaciados E_n = 2*pi*n/ln(b/a); <r> = {razon(EA):.3f} (peine rigido: 1.0)")
print(f"   densidad constante = {Lu/(2*np.pi):.3f}; ceros: a t=100 {np.log(100/2/np.pi)/(2*np.pi):.3f}, a t=1e6 {np.log(1e6/2/np.pi)/(2*np.pi):.3f}")

# B) rejilla en x, diferencias centrales, H = (xp + px)/2, contorno Dirichlet
x = np.linspace(a, b, N)
h = x[1] - x[0]
P = np.zeros((N, N), complex)
for i in range(N - 1):
    P[i, i + 1] = -1j / (2 * h)
    P[i + 1, i] = 1j / (2 * h)
X = np.diag(x)
H = (X @ P + P @ X) / 2
H = (H + H.conj().T) / 2
E = np.linalg.eigvalsh(H)
E = E[E > 1]
print(f"B) rejilla: {len(E)} niveles positivos entre {E.min():.2f} y {E.max():.2f}")
m = slice(len(E) // 4, 3 * len(E) // 4)
print(f"   <r> (centro) = {razon(E[m]):.3f}   (GUE 0.600, Poisson 0.386, peine 1.0)")
for q in (0.2, 0.5, 0.8):
    i = int(q * len(E))
    d = np.mean(np.diff(E[i - 20:i + 20]))
    print(f"   E={E[i]:8.2f}: densidad medida {1/d:.4f}   (1/2pi)ln(E/2pi) = {np.log(E[i]/2/np.pi)/(2*np.pi):.4f}")
