import sys
import numpy as np

# Uso: python3 psi_escalada.py [XMAX] [ARCHIVO_CEROS] [K]   (XMAX por defecto 2e6; en el movil 2e6 va sobrado)
XMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2_000_000
ceros = np.loadtxt(sys.argv[2] if len(sys.argv) > 2 else "ceros_100k.txt")
K = int(sys.argv[3]) if len(sys.argv) > 3 else 100
g = ceros[:K]

# psi(x) exacta: suma de ln p sobre p^k <= x
lam = np.zeros(XMAX + 1)
crib = np.ones(XMAX + 1, bool)
crib[:2] = False
for i in range(2, int(XMAX ** 0.5) + 1):
    if crib[i]:
        crib[i * i::i] = False
for p in np.nonzero(crib)[0]:
    q = p
    while q <= XMAX:
        lam[q] = np.log(p)
        q *= p
psi = np.cumsum(lam)

u = np.linspace(np.log(20), np.log(XMAX), 6000)
x = np.exp(u)
D = (psi[np.floor(x).astype(int)] - x) / np.sqrt(x)           # (psi(x) - x)/sqrt(x)


def D_ceros(gam):
    out = np.empty_like(x)
    for i in range(0, len(x), 500):
        lx = np.log(x[i:i + 500])[:, None]
        a = (0.5 * np.cos(gam * lx) + gam * np.sin(gam * lx)) / (0.25 + gam ** 2)
        out[i:i + 500] = -2 * a.sum(axis=1)
    return out - (np.log(2 * np.pi) + 0.5 * np.log(1 - x ** -2.0)) / np.sqrt(x)


DK = D_ceros(g)
rng = np.random.default_rng(0)
Dr = D_ceros(np.sort(rng.uniform(g[0], g[-1], K)))
print(f"x hasta {XMAX:.0e}, {K} ceros. Cociente (psi(x)-x)/sqrt(x):")
print(f"  rango medido [{D.min():.3f}, {D.max():.3f}],  max |D| = {np.abs(D).max():.3f}")
print(f"  RMS(D - reconstruccion con {K} ceros) = {np.sqrt(np.mean((D - DK) ** 2)):.3f}")
print(f"  RMS(D - control con {K} gammas al azar) = {np.sqrt(np.mean((D - Dr) ** 2)):.3f}")
print(f"  RMS(D) = {np.sqrt(np.mean(D ** 2)):.3f}")
m = x >= 101
cota = np.log(x) * (np.log(x) - np.log(np.log(x))) / (8 * np.pi)      # cota condicional a RH, x >= 101 (arXiv:2312.05628, Teor. 1.2)
print(f"  max |D| = {np.abs(D[m]).max():.3f} en x >= 101;  cota bajo RH en x = 1e7: {cota[-1]:.2f};  max |D| / cota = {np.max(np.abs(D[m]) / cota[m]):.3f}")
np.savez("psi_escalada.npz", u=u, x=x, D=D, DK=DK, Dr=Dr)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A, B, C = "#257de0", "#f47230", "#8a8a8a"
    fig, ax = plt.subplots(3, 1, figsize=(9, 10), dpi=150, sharex=True, gridspec_kw={"height_ratios": [3, 3, 2.4]})
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
        a.grid(axis="y", color="#e6e6e6", lw=0.6)
        a.tick_params(colors="#666666")
    ax[0].plot(u, D, color=A, lw=0.8, label=r"exacta: $(\psi(x)-x)/\sqrt{x}$")
    ax[0].plot(u, DK, color=B, lw=0.8, label=f"fórmula explícita con {K} ceros")
    ax[0].set_ylabel(r"$(\psi(x)-x)/\sqrt{x}$")
    ax[0].legend(frameon=False, loc="upper left", ncol=2)
    ax[0].set_title("Prueba y control: la oscilación de los primos viene de las posiciones de los ceros", loc="left", fontsize=11, color="#191919")
    ax[1].plot(u, D, color=A, lw=0.8, label="exacta")
    ax[1].plot(u, Dr, color=C, lw=0.8, label=f"control: {K} \u03b3 al azar en el mismo rango")
    ax[1].set_ylabel(r"$(\psi(x)-x)/\sqrt{x}$")
    ax[1].legend(frameon=False, loc="upper left", ncol=2)
    mm = x >= 101
    ax[2].plot(u[mm], np.maximum.accumulate(np.abs(D[mm])), color=A, lw=1.2, label=r"máx. acumulado de $|(\psi(x)-x)/\sqrt{x}|$")
    ax[2].plot(u[mm], cota[mm], color="#191919", lw=1.0, ls="--", label=r"cota bajo RH: $\ln x(\ln x-\ln\ln x)/8\pi$ (x$\geq$101)")
    ax[2].set_ylabel("cociente")
    ax[2].set_xlabel(r"$\ln x$")
    ax[2].legend(frameon=False, loc="upper left")
    fig.text(0.01, 0.005, f"Datos: primos exactos hasta {XMAX:.0e}; ceros de Odlyzko. La cota condicional es de arXiv:2312.05628; que se cumpla no prueba RH.", fontsize=7.5, color="#666666")
    fig.tight_layout(rect=(0, 0.015, 1, 1))
    fig.savefig("psi_escalada.png")
    print("figura guardada en psi_escalada.png")
except ImportError:
    print("(sin matplotlib: no se dibuja la figura; los numeros de arriba son el resultado)")
