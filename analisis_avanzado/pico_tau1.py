import sys
import numpy as np

# E05: de donde sale el pico K(tau=1) ~ 41 en 400 ventanas de 5000 ceros.
# Con u_n = N_suave(t_n) y d_n = (n - 1/2) - u_n (el conteo de Turing, ~S(t_n)):
#   exp(2 pi i u_n) = -exp(-2 pi i d_n)   =>   K(1) = |sum exp(2 pi i d_n)|^2 / N   (identidad exacta)
# Uso: python3 pico_tau1.py ceros_1millon.txt
def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t)

z = np.loadtxt(sys.argv[1] if len(sys.argv) > 1 else "ceros_1millon.txt")
Nw = 5000
u = theta(z) / np.pi + 1               # N_suave = theta/pi + 1
n = np.arange(1, len(z) + 1)
d = (n - 0.5) - u                      # d_n ~ S(t_n)
nv = len(z) // Nw
def K(tau, v):
    return np.abs(np.sum(np.exp(2j * np.pi * tau * v))) ** 2 / len(v)
for tau in [0.95, 1.00, 1.05, 2.00]:
    k = np.mean([K(tau, u[i * Nw:(i + 1) * Nw]) for i in range(nv)])
    print(f"K({tau:.2f}) medio en {nv} ventanas = {k:8.3f}")
k1 = np.array([K(1.0, u[i * Nw:(i + 1) * Nw]) for i in range(nv)])
k1d = np.array([np.abs(np.sum(np.exp(2j * np.pi * d[i * Nw:(i + 1) * Nw]))) ** 2 / Nw for i in range(nv)])
print(f"\nidentidad exacta  K(1) = |sum e^(2 pi i d_n)|^2 / N :  max |diferencia| entre ventanas = {np.abs(k1 - k1d).max():.2e}")
m = np.array([np.abs(np.mean(np.exp(2j * np.pi * d[i * Nw:(i + 1) * Nw]))) ** 2 for i in range(nv)])
print(f"|<exp(2 pi i d)>|^2 medio = {m.mean():.5f}  ->  N * eso = {Nw * m.mean():.2f}")
s2 = d.var()
print(f"var(d) = {s2:.4f}; gaussiano: exp(-4 pi^2 var) = {np.exp(-4 * np.pi ** 2 * s2):.5f} -> N * eso = {Nw * np.exp(-4 * np.pi ** 2 * s2):.1f}")
# Nota: el control 'barajado' de orbitas_primos.py (espaciados desordenados) da K(1.00) = 0.974; en CUE, E|Tr U^N|^2 = N, es decir K(1) = 1.
