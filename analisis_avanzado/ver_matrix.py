import sys
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

RUTA = sys.argv[1] if len(sys.argv) > 1 else "ceros_100k.txt"
rng = np.random.default_rng(7)


def cargar_ceros(ruta):
    z = np.loadtxt(ruta)
    return z[:, -1] if z.ndim > 1 else z


def unfold_riemann(t):
    # N_suave(t) = (t/2pi) ln(t/2pi e) + 7/8  -> densidad media 1
    return t / (2 * np.pi) * np.log(t / (2 * np.pi * np.e)) + 7.0 / 8.0


def gue_unfolded(n=800, mats=40, ventana=0.5):
    # H = (A + A^dagger)/2, A=(randn+i randn)/sqrt2  ->  radio del semicirculo R = sqrt(2n)
    R = np.sqrt(2.0 * n)
    F = lambda x: 0.5 + (x * np.sqrt(1 - x * x) + np.arcsin(x)) / np.pi
    trozos = []
    for _ in range(mats):
        A = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
        H = (A + A.conj().T) / 2
        x = np.linalg.eigvalsh(H) / R
        x = x[np.abs(x) < ventana]          # solo el centro del espectro
        trozos.append(n * F(x))             # unfolding con la ley del semicirculo
    return trozos


def espaciados(trozos):
    return np.concatenate([np.diff(u) for u in trozos])


def stats(nombre, s):
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    c = np.corrcoef(s[:-1], s[1:])[0, 1]
    print(f"{nombre:10s} n={len(s):6d}  <s>={s.mean():.3f}  <r>={r.mean():.4f}  "
          f"P(s<0.3)={100*np.mean(s<0.3):5.2f}%  corr(s_n,s_n+1)={c:+.3f}")


def sigma2(u_list, Ls, ventanas=20000):
    out = []
    for L in Ls:
        v = []
        for u in u_list:
            lo, hi = u[0], u[-1] - L
            if hi <= lo:
                continue
            x = rng.uniform(lo, hi, size=ventanas // len(u_list) + 1)
            cnt = np.searchsorted(u, x + L) - np.searchsorted(u, x)
            v.append(cnt)
        v = np.concatenate(v)
        out.append(v.var())
    return np.array(out)


def correlacion_pares(u, K=12, xmax=3.0, bins=60):
    h = np.zeros(bins)
    for m in range(1, K + 1):
        d = u[m:] - u[:-m]
        d = d[d < xmax]
        h += np.histogram(d, bins=bins, range=(0, xmax))[0]
    dx = xmax / bins
    return (np.arange(bins) + 0.5) * dx, h / (len(u) * dx)


def main():
    print("Cargando ceros...")
    t = cargar_ceros(RUTA)
    u = unfold_riemann(t)
    s_z = np.diff(u)
    print(f"{len(t)} ceros, altura t = {t[0]:.1f} .. {t[-1]:.1f}")

    print("Generando GUE (40 matrices de 800x800)...")
    trozos = gue_unfolded()
    s_g = espaciados(trozos)
    s_g = s_g / s_g.mean()

    print("\n=== Estadistica de espaciados (unfolding correcto en ambos) ===")
    stats("ceros", s_z)
    stats("GUE sim", s_g)
    stats("barajado", rng.permutation(s_z))
    stats("Poisson", rng.exponential(size=len(s_z)))
    print("(GUE asintotico: <r>=0.5996; Wigner P(s<0.3)=2.73%)")

    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        "P(s): ceros vs GUE vs Poisson", "Rigidez: varianza del conteo Sigma2(L)",
        "Correlacion de pares R2(x)"))

    bins = np.linspace(0, 3.5, 70)
    cz, bz = np.histogram(s_z, bins=bins, density=True)
    cg, _ = np.histogram(s_g, bins=bins, density=True)
    c = 0.5 * (bz[:-1] + bz[1:])
    x = np.linspace(0.001, 3.5, 300)
    fig.add_trace(go.Scatter(x=c, y=cz, mode="markers", name="ceros de Riemann",
                             marker=dict(color="magenta", size=5)), row=1, col=1)
    fig.add_trace(go.Scatter(x=c, y=cg, mode="lines", name="GUE (matrices simuladas)",
                             line=dict(color="cyan")), row=1, col=1)
    fig.add_trace(go.Scatter(x=x, y=np.exp(-x), mode="lines", name="Poisson",
                             line=dict(color="lime", dash="dot")), row=1, col=1)

    Ls = np.array([0.5, 1, 2, 3, 5, 8, 12, 20, 30, 50])
    s2z = sigma2([u], Ls)
    s2g = sigma2(trozos, Ls)
    g = np.euler_gamma
    teo = (np.log(2 * np.pi * Ls) + g + 1) / np.pi ** 2
    fig.add_trace(go.Scatter(x=Ls, y=s2z, mode="markers+lines", name="ceros",
                             marker=dict(color="magenta"), showlegend=False), row=1, col=2)
    fig.add_trace(go.Scatter(x=Ls, y=s2g, mode="lines", name="GUE sim",
                             line=dict(color="cyan"), showlegend=False), row=1, col=2)
    fig.add_trace(go.Scatter(x=Ls, y=teo, mode="lines", name="teoria GUE (L grande)",
                             line=dict(color="white", dash="dash"), showlegend=False), row=1, col=2)
    fig.add_trace(go.Scatter(x=Ls, y=Ls, mode="lines", name="Poisson (=L)",
                             line=dict(color="lime", dash="dot"), showlegend=False), row=1, col=2)
    fig.update_xaxes(type="log", row=1, col=2)
    fig.update_yaxes(range=[0, 1.5], row=1, col=2)

    xc, r2 = correlacion_pares(u)
    xt = np.linspace(0.001, 3, 300)
    fig.add_trace(go.Scatter(x=xc, y=r2, mode="markers", name="ceros",
                             marker=dict(color="magenta", size=5), showlegend=False), row=1, col=3)
    fig.add_trace(go.Scatter(x=xt, y=1 - np.sinc(xt) ** 2, mode="lines",
                             name="Montgomery 1 - sinc^2", line=dict(color="white", dash="dash"),
                             showlegend=False), row=1, col=3)

    fig.update_layout(template="plotly_dark", height=520,
                      title="Ceros de Riemann vs GUE: datos reales, sin formulas pintadas",
                      legend=dict(orientation="h", y=-0.15))
    fig.write_html("matrix.html")
    print("\nListo: matrix.html")


if __name__ == "__main__":
    main()
