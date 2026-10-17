#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot3d_html.py — Visualización interactiva de zeta(1/2+it)

Lee los CSV generados por pfi01.py y produce UN HTML autocontenido con
tres paneles Plotly:

  1. Trayectoria 3D (X, Y, t), color = log10|zeta(1/2+it)|
  2. Proyección 2D (X, Y), color = t
  3. Z(t) real vs. t, con los gamma_n marcados

Los gamma_n calculados por mpmath aparecen como marcas cian, de modo que
se puede verificar visualmente el alineamiento.

v1.1: panel_2d ahora dibuja "lines+markers" en vez de solo "markers".
Sin la línea, una espiral que vuelve una y otra vez por regiones parecidas
del plano (algo esperado, no un error) es visualmente indistinguible de
puntos sueltos al azar. El color por t sigue yendo en los marcadores
(Scattergl.line.color no acepta array; marker.color si).

Requisitos:
  pip install plotly numpy        # Python puro, sin compilación (apto Termux)

Uso:
  python3 plot3d_html.py 1000
  termux-open trayectoria_N1000.html
"""

import csv
import os
import sys

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- parámetros -------------------------------------------------------------
N        = int(sys.argv[1]) if len(sys.argv) > 1 else 300
OUT_DIR  = f"salida_pfi01_N{N}"
MAX_SEG  = 8000            # downsample para el panel 3D
OFFLINE  = True            # True: embebe plotly.js (~3 MB, abre sin internet)
                            # False: usa CDN (archivo chico, requiere conexión)
OUT_HTML = f"trayectoria_N{N}.html"


# --- carga -------------------------------------------------------------------
def cargar_trayectoria(ruta):
    with open(ruta) as f:
        rows = list(csv.DictReader(f))
    t = np.array([float(r["t"])   for r in rows])
    X = np.array([float(r["X_t"]) for r in rows])
    Y = np.array([float(r["Y_t"]) for r in rows])
    Z = np.array([float(r["Z_t"]) for r in rows])
    return t, X, Y, Z


def cargar_ceros(ruta):
    with open(ruta) as f:
        return np.array([float(r["t_n"]) for r in csv.DictReader(f)])


# --- paneles ------------------------------------------------------------------
def panel_3d(t, X, Y, gammas):
    R    = np.sqrt(X**2 + Y**2)
    logR = np.log10(np.clip(R, 1e-3, None))
    tmax = float(t.max())

    if len(t) > MAX_SEG:
        idx = np.linspace(0, len(t) - 1, MAX_SEG).astype(int)
        t, X, Y, logR = t[idx], X[idx], Y[idx], logR[idx]

    traza = go.Scatter3d(
        x=X, y=Y, z=t, mode="lines",
        line=dict(color=logR, colorscale="Plasma", width=2,
                  colorbar=dict(title="log10|zeta|", x=1.00,
                                len=0.55, thickness=12, y=0.78)),
        name="Gamma(t)",
        hovertemplate="t=%{z:.3f}<br>X=%{x:.4f}<br>Y=%{y:.4f}<extra></extra>",
    )
    marcas = gammas[gammas <= tmax][:100]
    eje_t = go.Scatter3d(
        x=np.zeros_like(marcas), y=np.zeros_like(marcas), z=marcas,
        mode="markers", marker=dict(size=3, color="cyan"),
        name="gamma_n",
        hovertemplate="gamma=%{z:.4f}<extra></extra>",
    )
    lim = float(max(np.abs(X).max(), np.abs(Y).max()) or 1.0)
    return traza, eje_t, lim, tmax


def panel_2d(t, X, Y):
    if len(t) > MAX_SEG * 3:
        idx = np.linspace(0, len(t) - 1, MAX_SEG * 3).astype(int)
        t, X, Y = t[idx], X[idx], Y[idx]
    # Scattergl.line.color NO acepta arrays (solo Scatter3d lo hace);
    # marker.color si admite un array para la escala de color. La linea
    # gris solo da continuidad visual para distinguir "espiral que se
    # cruza consigo misma" de "puntos sueltos".
    traza = go.Scattergl(
        x=X, y=Y, mode="lines+markers",
        line=dict(color="rgba(130,130,130,0.25)", width=1),
        marker=dict(
            color=t, colorscale="Turbo", size=2.5,
            colorbar=dict(title="t", x=1.00, len=0.55,
                          thickness=12, y=0.22),
        ),
        name="(X,Y)",
        hovertemplate="X=%{x:.4f}<br>Y=%{y:.4f}<extra></extra>",
    )
    origen = go.Scatter(
        x=[0], y=[0], mode="markers",
        marker=dict(size=14, color="red", symbol="x"),
        name="origen (cero)",
        hovertemplate="origen<extra></extra>",
    )
    lim = float(max(np.abs(X).max(), np.abs(Y).max()) or 1.0)
    return traza, origen, lim


def panel_Z(t, Z, gammas):
    traza = go.Scatter(
        x=t, y=Z, mode="lines",
        line=dict(color="steelblue", width=1.0),
        name="Z(t)",
        hovertemplate="t=%{x:.3f}<br>Z=%{y:.4f}<extra></extra>",
    )
    marcadores = gammas[gammas <= t.max()][:500]
    ceros_z = go.Scatter(
        x=marcadores, y=np.zeros_like(marcadores),
        mode="markers", marker=dict(size=5, color="cyan", symbol="circle"),
        name="gamma_n",
        hovertemplate="gamma=%{x:.4f}<extra></extra>",
    )
    cero_h = go.Scatter(
        x=[float(t.min()), float(t.max())], y=[0, 0],
        mode="lines", line=dict(color="black", width=0.8, dash="dot"),
        showlegend=False, hoverinfo="skip",
    )
    return traza, ceros_z, cero_h


# --- ensamblado ----------------------------------------------------------------
def construir_html(t, X, Y, Z, gammas):
    traza3d, eje_t, lim, tmax = panel_3d(t, X, Y, gammas)
    tr2d, origen, lim2d        = panel_2d(t, X, Y)
    trZ, ceros_z, cero_h       = panel_Z(t, Z, gammas)

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "scene"}, {"type": "xy"}],
               [{"type": "xy", "colspan": 2}, None]],
        subplot_titles=(
            f"Trayectoria 3D — (X, Y, t) · N_ZEROS={N}",
            "Proyección (X, Y) · color = t",
            "Z(t) real vs. t · marcas cian en cada gamma_n",
        ),
        vertical_spacing=0.10,
        horizontal_spacing=0.10,
    )

    fig.add_trace(traza3d, row=1, col=1)
    fig.add_trace(eje_t,   row=1, col=1)
    fig.add_trace(tr2d,    row=1, col=2)
    fig.add_trace(origen,  row=1, col=2)
    fig.add_trace(cero_h,  row=2, col=1)
    fig.add_trace(trZ,     row=2, col=1)
    fig.add_trace(ceros_z, row=2, col=1)

    fig.update_scenes(
        xaxis_title="X(t) = Re zeta(1/2+it)",
        yaxis_title="Y(t) = Im zeta(1/2+it)",
        zaxis_title="t",
        xaxis=dict(range=[-lim, lim]),
        yaxis=dict(range=[-lim, lim]),
        zaxis=dict(range=[0, tmax]),
        aspectmode="manual",
        aspectratio=dict(x=1, y=1, z=2.2),
        camera=dict(eye=dict(x=1.5, y=1.5, z=0.9)),
    )
    fig.update_xaxes(title_text="X(t)", range=[-lim2d, lim2d], row=1, col=2)
    fig.update_yaxes(title_text="Y(t)", range=[-lim2d, lim2d],
                      scaleanchor="x", scaleratio=1, row=1, col=2)
    fig.update_xaxes(title_text="t", row=2, col=1)
    fig.update_yaxes(title_text="Z(t)", row=2, col=1)

    fig.update_layout(
        height=1050,
        margin=dict(l=50, r=120, t=110, b=60),
        showlegend=True,
        legend=dict(orientation="h", y=-0.06, x=0.5, xanchor="center"),
        template="plotly_white",
        title=dict(
            text=(f"<b>zeta(1/2+it)</b> — N_ZEROS={N}, "
                  f"{len(t)} muestras en t in [{t.min():.2f}, {t.max():.2f}]<br>"
                  f"<sub>El pulso fractal del infinito · visualización interactiva · "
                  f"Benjamín Cabeza Durán</sub>"),
            x=0.5, xanchor="center",
        ),
    )

    return fig.to_html(
        include_plotlyjs=(True if OFFLINE else "cdn"),
        full_html=True,
        config=dict(responsive=True, displaylogo=False, scrollZoom=True),
    )


# --- main ----------------------------------------------------------------------
def main():
    ruta_traj  = os.path.join(OUT_DIR, "trayectoria.csv")
    ruta_ceros = os.path.join(OUT_DIR, "ceros.csv")
    if not os.path.exists(ruta_traj):
        sys.exit(f"falta {ruta_traj} — corre pfi01.py con N_ZEROS={N} primero")

    print(f"plot3d_html.py — leyendo {OUT_DIR}/")
    t, X, Y, Z = cargar_trayectoria(ruta_traj)
    gammas     = cargar_ceros(ruta_ceros)
    print(f"  {len(t)} muestras, t in [{t.min():.2f}, {t.max():.2f}]")
    print(f"  {len(gammas)} ceros, ultimo gamma = {gammas[-1]:.3f}")

    html = construir_html(t, X, Y, Z, gammas)
    with open(OUT_HTML, "w") as f:
        f.write(html)
    kb = os.path.getsize(OUT_HTML) / 1024
    print(f"  guardado: {OUT_HTML} ({kb:.0f} KB)")
    print(f"\nAbre con:  termux-open {OUT_HTML}")


if __name__ == "__main__":
    main()
