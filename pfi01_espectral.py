#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pfi01_espectral.py — Estadistica espectral (RMT) sobre los ceros ya calculados

Lee salida_pfi01_N{N}/ceros.csv (lo que ya produjo pfi01.py: NO vuelve a
tocar zeta, es barato) y calcula:

  - P(s): histograma de espaciados normalizados vs. GUE y Poisson
  - <r>: razon de espaciados consecutivos (NO requiere unfolding)
  - [opcional, off por defecto] Delta_3 y Sigma^2 (largo alcance)

Estado de la discusion (ver manuscrito, seccion "Datos, soporte y
visualizacion"): con N=5000, <r>=0.616 esta ~4 sigma por encima del
valor asintotico GUE (0.5996) -- NO es un error numerico (la direccion
es la opuesta a la que producirian errores de redondeo) y es consistente
en direccion y magnitud aproximada con la correccion de altura finita de
Bogomolny-Keating (1996) / Berry-Keating (1999): a T finito, <r> converge
a GUE como ~1/log^2(T), y nuestras dos ventanas de t lo confirman
cualitativamente. Por eso Delta_3/Sigma^2 se tratan aqui como
confirmacion OPCIONAL, no como paso obligatorio.

Uso:
  python3 pfi01_espectral.py 5000              # P(s) y <r> (rapido)
  python3 pfi01_espectral.py 5000 --delta3     # + Delta_3/Sigma^2
"""

import csv
import os
import sys

import numpy as np
import mpmath as mp

mp.mp.dps = 20


def cargar_ceros(out_dir):
    ruta = os.path.join(out_dir, "ceros.csv")
    with open(ruta) as f:
        rows = list(csv.DictReader(f))
    return np.array([float(r["t_n"]) for r in rows])


# ---------------------------------------------------------------------------
# P(s) vs GUE / Poisson (requiere unfolding local, densidad ~log(t/2pi)/2pi)
# ---------------------------------------------------------------------------
def p_gue(s):
    return (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)


def p_poisson(s):
    return np.exp(-s)


def espaciados_normalizados(ts):
    deltas = np.diff(ts)
    t_medio = 0.5 * (ts[:-1] + ts[1:])
    densidad_local = np.log(t_medio / (2 * np.pi)) / (2 * np.pi)
    return deltas * densidad_local  # = delta / (1/densidad_local)


def analisis_Ps(ts, out_dir):
    s = espaciados_normalizados(ts)
    bins = np.linspace(0, 3, 13)
    hist, edges = np.histogram(s, bins=bins, density=True)
    centros = (edges[:-1] + edges[1:]) / 2

    filas = list(zip(centros, hist, p_gue(centros), p_poisson(centros)))
    with open(os.path.join(out_dir, "espectral_Ps.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s", "empirico", "GUE", "Poisson"])
        w.writerows(filas)

    frac = (s < 0.3).mean()
    s_fino = np.linspace(0, 0.3, 200000)
    p_gue_03 = np.trapezoid(p_gue(s_fino), s_fino)
    se = np.sqrt(p_gue_03 * (1 - p_gue_03) / len(s))
    print(f"P(s): n={len(s)}  media(s)={s.mean():.4f} (debe rondar 1.0)")
    print(f"  fraccion s<0.3 = {frac:.4f}  (GUE={p_gue_03:.4f}, "
          f"Poisson={1-np.exp(-0.3):.4f}, desviacion={(p_gue_03-frac)/se:.1f} sigma)")


# ---------------------------------------------------------------------------
# <r>: razon de espaciados consecutivos (sin unfolding)
# ---------------------------------------------------------------------------
def analisis_r(ts):
    deltas = np.diff(ts)
    r = np.minimum(deltas[:-1], deltas[1:]) / np.maximum(deltas[:-1], deltas[1:])
    n = len(r)
    media, sd = r.mean(), r.std(ddof=1)
    r0 = r - media
    rho1 = np.dot(r0[:-1], r0[1:]) / np.dot(r0, r0)
    se = (sd / np.sqrt(n)) * np.sqrt(max(1 + 2 * rho1, 1e-6))  # corregido por autocorr. lag-1
    gue = 0.5996
    print(f"<r>: n={n}  media={media:.4f}  sd={sd:.4f}  rho1={rho1:.3f}")
    print(f"  SE (corregido por autocorrelacion) = {se:.5f}  "
          f"desviacion = {(media-gue)/se:.2f} sigma  (GUE={gue}, Poisson=0.386, GOE=0.5307)")
    return media


# ---------------------------------------------------------------------------
# OPCIONAL: Delta_3 y Sigma^2 (largo alcance) -- confirmacion, no obligatorio
#
# Formula GUE verificada por triangulacion en 3 fuentes independientes
# (todas coinciden exactamente):
#   Delta_3_GUE(L) = ( ln(2*pi*L) + gamma - 5/4 ) / (2*pi^2)
#                  = ( ln(L) + 1.1651 ) / (2*pi^2)
# OJO: una version anterior de esta formula circulo con el prefactor de
# GOE (1/pi^2) y luego con la constante de Sigma^2 (gamma+1-pi^2/8) mal
# traspuesta a Delta_3 -- ninguna de esas dos coincidia ni siquiera entre
# si. Esta es la que de verdad aparece igual en multiples fuentes.
# Poisson: Delta_3(L) = L/15.
# ---------------------------------------------------------------------------
GAMMA_EULER = 0.5772156649


def delta3_gue(L):
    return (np.log(L) + 1.1651) / (2 * np.pi**2)


def unfold_riemann_von_mangoldt(ts):
    """x_i = N(t_i), con N(t) = (t/2pi)log(t/2pi) - t/2pi + 7/8 (parte suave,
    sin el termino oscilatorio S(t): eso es justo la fluctuacion que se
    quiere medir, no se debe restar)."""
    t = ts / (2 * np.pi)
    return t * np.log(t) - t + 7.0 / 8.0


def analisis_delta3_sigma2(ts, out_dir, L_vals=(5, 10, 20, 35, 50), puntos_grilla=400):
    # Delta_3 es una integral sobre epsilon de (escalera - A*eps - B)^2, NO una
    # suma sobre los m puntos muestreados con igual peso (eso fue el bug de la
    # primera version: subestimaba Delta_3 porque ignoraba el ancho de cada
    # escalon). Se aproxima evaluando la escalera en una grilla fina dentro de
    # cada ventana y ajustando A,B por minimos cuadrados sobre esa grilla.
    x = unfold_riemann_von_mangoldt(ts)
    filas = []
    for L in L_vals:
        d3_vals, s2_vals = [], []
        paso = L / 4.0
        a = x[0]
        while a + L <= x[-1]:
            ventana = x[(x >= a) & (x < a + L)]
            if len(ventana) >= 3:
                m = len(ventana)
                s2_vals.append(m)

                grilla = np.linspace(a, a + L, puntos_grilla)
                escalera = np.searchsorted(ventana, grilla, side="right")
                A_mat = np.vstack([grilla - a, np.ones_like(grilla)]).T
                coef, *_ = np.linalg.lstsq(A_mat, escalera, rcond=None)
                resid = escalera - (A_mat @ coef)
                # mean(resid**2) sobre una grilla uniforme en [a,a+L] YA
                # aproxima (1/L)*integral(resid^2 d-eps) -- esa es la
                # definicion de Delta_3, no hay que dividir por L otra vez
                # (ese doble conteo fue el bug real de la primera version).
                d3_vals.append(np.mean(resid**2))
            a += paso
        if d3_vals:
            d3 = float(np.mean(d3_vals))
            s2 = float(np.var(s2_vals))
            filas.append((L, d3, delta3_gue(L), L / 15, s2, len(d3_vals)))

    with open(os.path.join(out_dir, "espectral_delta3_sigma2.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L", "Delta3_empirico", "Delta3_GUE", "Delta3_Poisson", "Sigma2_empirico", "n_ventanas"])
        w.writerows(filas)

    print(f"{'L':>5} {'Delta3 emp':>11} {'Delta3 GUE':>11} {'Delta3 Poi':>11} {'Sigma2 emp':>11} {'n_vent':>7}")
    for L, d3, d3g, d3p, s2, nv in filas:
        print(f"{L:5.0f} {d3:11.4f} {d3g:11.4f} {d3p:11.4f} {s2:11.4f} {nv:7d}")


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    hacer_delta3 = "--delta3" in sys.argv
    out_dir = f"salida_pfi01_N{N}"

    ts = cargar_ceros(out_dir)
    print(f"{len(ts)} ceros cargados de {out_dir}/ceros.csv, t hasta {ts[-1]:.1f}\n")

    analisis_Ps(ts, out_dir)
    print()
    analisis_r(ts)

    if hacer_delta3:
        print("\n[opcional] Delta_3 / Sigma^2 (confirmacion, no obligatorio):")
        analisis_delta3_sigma2(ts, out_dir)


if __name__ == "__main__":
    main()
