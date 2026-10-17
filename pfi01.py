#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PFI-01 — Laboratorio computacional de "El pulso fractal del infinito"
Programa de investigación de Benjamín Cabeza Durán sobre la Hipótesis de Riemann.

Qué calcula (todo se guarda en CSV para que cada afirmación del manuscrito
pueda reconstruirse desde los datos):
  1. ceros.csv            - primeros N_ZEROS ceros no triviales (n, t_n)
  2. separaciones.csv     - huecos entre ceros consecutivos (Delta_n)
  3. trayectoria.csv      - Gamma(t)=zeta(1/2+it): X(t), Y(t), R(t), fase, Z(t)
  4. dimension_caja.csv   - conteo de cajas N(eps) para estimar dimension fractal
  5. autocorrelacion.csv  - autocorrelacion normalizada de Z(t)
  6. espectro_entropia.csv- espectro de potencia de Z(t) + entropia espectral
  7. gue_comparacion.csv  - espaciados normalizados vs. la prediccion GUE

Requisitos: pip install mpmath numpy
  (en Termux, si numpy falla por compilacion: pkg install python-numpy)

Uso:
  python3 pfi01.py

IMPORTANTE - alcance de este script:
  Esto NO demuestra ni refuta nada sobre la Hipotesis de Riemann. Es la
  version PFI-01, pensada para correr en un telefono. Los resultados
  (sobre todo la comparacion GUE) son ilustrativos salvo que N_ZEROS sea
  de miles. Cada numero que produce este script es un dato a interpretar,
  no una conclusion.

  v0.2: la resolucion de la trayectoria ahora se fija como densidad
  (DENSIDAD_T, puntos por unidad de t) en vez de un total fijo de puntos.
  En la v0.1, el total estaba fijado en 2000 sin importar t_max, asi que
  al subir N_ZEROS la densidad real caia sin avisar — dos corridas con
  distinto N_ZEROS no eran comparables entre si (ni la dimension de caja
  ni la entropia espectral, que dependen de la resolucion). Si ya corriste
  la v0.1 con algun N_ZEROS, esos resultados no son comparables con los
  que produzca esta version.
"""

import csv
import os
import time

import mpmath as mp
import numpy as np

# ---------------------------------------------------------------------------
# Parametros (ajustables)
# ---------------------------------------------------------------------------
N_ZEROS = 1000              # cuantos ceros no triviales calcular
DENSIDAD_T = 15.0          # puntos de muestreo por unidad de t (NO total de
                           # puntos): esto es lo que hay que mantener fijo
                           # para comparar [0,100] vs [0,1000] etc. de forma
                           # honesta. Antes la resolucion total estaba fija
                           # en 2000 puntos sin importar el rango de t, asi
                           # que al subir N_ZEROS la densidad caia en
                           # silencio y contaminaba la comparacion.
MAX_LAG_AUTOCORR = 150     # lags maximos para la autocorrelacion
EPSILONS_CAJA = [1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256]  # para dimension de caja
FACTOR_DENSIFICADO = 20   # sub-puntos interpolados entre cada par de muestras
                          # consecutivas antes de contar cajas (ver nota en
                          # dimension_caja: sin esto, a eps chico se cuentan
                          # muestras sueltas, no la curva)
DPS = 20                   # digitos de precision de mpmath
OUT_DIR = f"salida_pfi01_N{N_ZEROS}"  # cada escala de N_ZEROS va a su propia
                                       # carpeta, para poder comparar
                                       # [0,100] vs [0,1000] vs [0,10^4]...
                                       # sin pisar corridas anteriores

mp.mp.dps = DPS


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
def guardar_csv(nombre, cabecera, filas):
    ruta = os.path.join(OUT_DIR, nombre)
    with open(ruta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cabecera)
        w.writerows(filas)
    print(f"  guardado: {ruta}  ({len(filas)} filas)")


# ---------------------------------------------------------------------------
# 1-2. Ceros y separaciones
# ---------------------------------------------------------------------------
def calcular_ceros(n_zeros):
    print(f"[1/7] Calculando los primeros {n_zeros} ceros no triviales...")
    ceros = []
    for n in range(1, n_zeros + 1):
        rho = mp.zetazero(n)
        ceros.append((n, float(rho.imag)))
    guardar_csv("ceros.csv", ["n", "t_n"], ceros)
    return ceros


def calcular_separaciones(ceros):
    print("[2/7] Calculando separaciones entre ceros consecutivos...")
    filas = []
    for i in range(len(ceros) - 1):
        n, t_n = ceros[i]
        _, t_n1 = ceros[i + 1]
        filas.append((n, t_n, t_n1, t_n1 - t_n))
    guardar_csv("separaciones.csv", ["n", "t_n", "t_n+1", "delta_n"], filas)
    return filas


# ---------------------------------------------------------------------------
# 3. Trayectoria Gamma(t) = zeta(1/2+it) y Z(t) de Riemann-Siegel
# ---------------------------------------------------------------------------
def calcular_trayectoria(t_max, n_puntos):
    print(f"[3/7] Muestreando Gamma(t) y Z(t) en [0, {t_max:.1f}] "
          f"({n_puntos} puntos; puede tardar un rato)...")
    ts = np.linspace(1e-4, t_max, n_puntos)
    filas = []
    X, Y, Zs = [], [], []
    fase_prev = None
    fase_acum = 0.0
    dos_pi = float(2 * mp.pi)
    t0 = time.time()
    for i, t in enumerate(ts):
        s = mp.mpc(0.5, float(t))
        z = mp.zeta(s)
        x, y = float(z.real), float(z.imag)
        r = float(mp.sqrt(z.real * z.real + z.imag * z.imag))
        phi = float(mp.atan2(z.imag, z.real))
        theta = float(mp.siegeltheta(float(t)))
        # Z(t) = Re[ e^{i theta(t)} * zeta(1/2+it) ]; su parte imaginaria
        # deberia ser ~0 (es la identidad del apendice, seccion 3).
        Zt = x * float(mp.cos(theta)) - y * float(mp.sin(theta))

        if fase_prev is not None:
            dphi = phi - fase_prev
            if dphi > mp.pi:
                dphi -= dos_pi
            elif dphi < -mp.pi:
                dphi += dos_pi
            fase_acum += dphi
        else:
            fase_acum = phi
        fase_prev = phi

        filas.append((float(t), x, y, r, phi, fase_acum, Zt))
        X.append(x)
        Y.append(y)
        Zs.append(Zt)

        if (i + 1) % 500 == 0:
            print(f"    {i+1}/{n_puntos} puntos "
                  f"({time.time()-t0:.1f}s transcurridos)")

    guardar_csv(
        "trayectoria.csv",
        ["t", "X_t", "Y_t", "R_t", "fase", "fase_desenrollada", "Z_t"],
        filas,
    )
    return np.array(ts), np.array(X), np.array(Y), np.array(Zs)


# ---------------------------------------------------------------------------
# 4. Dimension de caja (box-counting) de (X(t), Y(t))
# ---------------------------------------------------------------------------
def dimension_caja(X, Y, epsilons, factor_densificado=FACTOR_DENSIFICADO):
    print("[4/7] Estimando dimension de caja de la trayectoria (X,Y)...")
    # Problema a evitar: con solo los puntos ya calculados, a partir de
    # cierto eps chico cada celda ocupada tiende a contener como mucho una
    # muestra -> N(eps) deja de medir la curva y empieza a medir cuantos
    # puntos muestreamos. La senal es una pendiente local que se desploma
    # en los eps mas finos. Por eso se interpola linealmente entre muestras
    # consecutivas (no hace falta volver a evaluar zeta: son solo rectas
    # entre puntos ya calculados) antes de contar cajas.
    idx = np.arange(len(X))
    idx_fino = np.linspace(0, len(X) - 1, (len(X) - 1) * factor_densificado + 1)
    Xf = np.interp(idx_fino, idx, X)
    Yf = np.interp(idx_fino, idx, Y)

    xmin, xmax = Xf.min(), Xf.max()
    ymin, ymax = Yf.min(), Yf.max()
    escala = max(xmax - xmin, ymax - ymin) or 1.0
    Xn = (Xf - xmin) / escala
    Yn = (Yf - ymin) / escala
    total_puntos = len(Xn)

    filas = []
    puntos_log = []
    avisos = []
    for eps in epsilons:
        celdas = np.unique(
            np.stack([np.floor(Xn / eps), np.floor(Yn / eps)], axis=1), axis=0
        )
        n_eps = len(celdas)
        filas.append((eps, n_eps))
        puntos_log.append((np.log(1.0 / eps), np.log(n_eps)))
        if n_eps > 0.3 * total_puntos:
            avisos.append(eps)

    guardar_csv("dimension_caja.csv", ["epsilon", "N_epsilon"], filas)
    print("    eps        N(eps)     log(1/eps)   log N(eps)")
    for (eps, n_eps), (lx, ly) in zip(filas, puntos_log):
        marca = "  <- limitado por muestreo" if eps in avisos else ""
        print(f"    {eps:<10.6f} {n_eps:<10d} {lx:<12.3f} {ly:<.3f}{marca}")

    # ajusta la pendiente solo con los eps que NO estan marcados como
    # limitados por muestreo (si hay al menos 3 puntos limpios)
    limpios = [(lx, ly) for (eps, _), (lx, ly) in zip(filas, puntos_log)
               if eps not in avisos]
    fuente = limpios if len(limpios) >= 3 else puntos_log
    xs = np.array([p[0] for p in fuente])
    ys = np.array([p[1] for p in fuente])
    pendiente, intercepto = np.polyfit(xs, ys, 1)

    if avisos:
        print(f"    ({len(avisos)} de {len(epsilons)} valores de eps estan "
              f"limitados por muestreo (N(eps)>30% de los {total_puntos} "
              f"puntos densificados) y se excluyeron del ajuste)")
    print(f"    pendiente log N(eps) vs log(1/eps) = {pendiente:.3f}  "
          f"(descriptiva; NO es una prueba de fractalidad)")
    return pendiente


# ---------------------------------------------------------------------------
# 5. Autocorrelacion de Z(t)
# ---------------------------------------------------------------------------
def autocorrelacion(Z, max_lag):
    print("[5/7] Calculando autocorrelacion de Z(t)...")
    z = Z - np.mean(Z)
    n = len(z)
    var = np.dot(z, z) / n
    filas = []
    for lag in range(min(max_lag, n - 1) + 1):
        c = np.dot(z[: n - lag], z[lag:]) / (n - lag) / var
        filas.append((lag, c))
    guardar_csv("autocorrelacion.csv", ["lag", "autocorrelacion"], filas)


# ---------------------------------------------------------------------------
# 6. Espectro de potencia y entropia espectral de Z(t)
# ---------------------------------------------------------------------------
def espectro_y_entropia(Z, dt):
    print("[6/7] Calculando espectro de potencia y entropia espectral...")
    z = Z - np.mean(Z)
    espectro = np.abs(np.fft.rfft(z)) ** 2
    freqs = np.fft.rfftfreq(len(z), d=dt)
    p = espectro / espectro.sum()
    p_nz = p[p > 0]
    H = float(-np.sum(p_nz * np.log(p_nz)))
    H_max = float(np.log(len(p_nz)))
    H_norm = H / H_max if H_max > 0 else float("nan")

    filas = [(float(f), float(pw)) for f, pw in zip(freqs, espectro)]
    guardar_csv("espectro_entropia.csv", ["frecuencia", "potencia"], filas)
    print(f"    entropia espectral normalizada = {H_norm:.4f}  "
          f"(0 = muy concentrado, 1 = tipo ruido blanco)")
    return H, H_norm


# ---------------------------------------------------------------------------
# 7. Comparacion con estadistica GUE (espaciados normalizados)
# ---------------------------------------------------------------------------
def comparacion_gue(ceros):
    print("[7/7] Comparando espaciados normalizados con la prediccion GUE...")
    filas = []
    for i in range(len(ceros) - 1):
        _, t_n = ceros[i]
        _, t_n1 = ceros[i + 1]
        delta = t_n1 - t_n
        # densidad media local de ceros: dN/dT ~ (1/2pi) log(T/2pi)
        # => espaciado medio local ~ 2*pi / log(T/2pi)
        t_medio = 0.5 * (t_n + t_n1)
        densidad_local = float(mp.log(t_medio / (2 * mp.pi)) / (2 * mp.pi))
        espaciado_medio = 1.0 / densidad_local if densidad_local > 0 else float("nan")
        s = delta / espaciado_medio
        p_gue = (32.0 / np.pi ** 2) * s ** 2 * np.exp(-4.0 * s ** 2 / np.pi)
        filas.append((t_medio, delta, espaciado_medio, s, p_gue))
    guardar_csv(
        "gue_comparacion.csv",
        ["t_medio", "delta_n", "espaciado_medio_local", "s_normalizado", "P_GUE(s)"],
        filas,
    )
    if len(ceros) - 1 < 500:
        print(f"    AVISO: solo hay {len(ceros)-1} espaciados. "
              f"Para una comparacion GUE con algun valor estadistico "
              f"hacen falta miles (sube N_ZEROS). Esto es un smoke test.")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"PFI-01 — dps={DPS}, N_ZEROS={N_ZEROS}, "
          f"densidad={DENSIDAD_T} puntos/unidad de t\n")

    t0 = time.time()
    ceros = calcular_ceros(N_ZEROS)
    calcular_separaciones(ceros)

    t_max = ceros[-1][1] + 5.0
    n_puntos = max(500, int(DENSIDAD_T * t_max))
    dt = t_max / n_puntos
    print(f"    (t_max={t_max:.1f} -> {n_puntos} puntos a densidad fija; "
          f"esto escala con N_ZEROS, a diferencia de la version anterior)")
    ts, X, Y, Z = calcular_trayectoria(t_max, n_puntos)

    dimension_caja(X, Y, EPSILONS_CAJA)
    autocorrelacion(Z, MAX_LAG_AUTOCORR)
    espectro_y_entropia(Z, dt)
    comparacion_gue(ceros)

    print(f"\nListo en {time.time()-t0:.1f}s. "
          f"7 archivos CSV en ./{OUT_DIR}/")


if __name__ == "__main__":
    main()
