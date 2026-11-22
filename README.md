# El pulso fractal del infinito — PFI-01

Programa de investigación de Benjamín Cabeza Durán: lectura geométrica y
espectral de la Hipótesis de Riemann. El ensayo y el apéndice matemático
están en [`docs/manuscrito.md`](docs/manuscrito.md); esta carpeta en
general es el pipeline computacional que sostiene esa lectura.

## Dos lecturas distintas

- **Lectura geométrica** (trayectoria Γ(t), dimensión de caja, autocorrelación,
  entropía espectral): descriptiva, depende de la resolución de muestreo, sin
  predicción teórica cerrada contra la cual comparar.
- **Lectura espectral** (espaciados normalizados, ⟨r⟩, Δ₃, Σ²): comparación
  contra predicciones exactas de GUE y Poisson. Con N_ZEROS≥5000 el acuerdo
  con GUE tiene peso estadístico y las desviaciones son cuantificables en σ.

## Contenido

- `docs/manuscrito.md` — el ensayo, el apéndice matemático y los resultados
  espectrales, en un solo documento versionado (v0.3).
- `pfi01.py` — ceros, separaciones, trayectoria Γ(t)=ζ(1/2+it), dimensión de
  caja, autocorrelación, entropía espectral, comparación GUE básica.
  Parámetros ajustables al principio del archivo (N_ZEROS, DENSIDAD_T, DPS).
- `pfi01_espectral.py` — estadística espectral (RMT) sobre los ceros que ya
  calculó `pfi01.py`: no vuelve a tocar zeta, es barato. P(s) y ⟨r⟩ por
  defecto; Δ₃/Σ² con `--delta3` (opcional, confirmatorio — ver más abajo).
- `plot3d_html.py` — lee los CSV de `pfi01.py` y arma un HTML con 3 paneles
  Plotly (trayectoria 3D, proyección 2D, Z(t) con los ceros marcados).
- `salida_pfi01_N1000/`, `salida_pfi01_N5000/` — salidas de `pfi01.py` para
  cada escala. N=5000 (4999 separaciones) es la que sostiene los hallazgos
  espectrales de abajo.
- `trayectoria_N1000.html` — visualización 3D/2D/Z(t) para esa corrida.

## Uso

```
pip install mpmath numpy plotly
python3 pfi01.py                       # ajusta N_ZEROS arriba del archivo
python3 pfi01_espectral.py 5000        # P(s) + <r>
python3 pfi01_espectral.py 5000 --delta3   # + Delta_3/Sigma^2
python3 plot3d_html.py 5000
termux-open trayectoria_N5000.html
```

## Alcance

Nada de esto demuestra ni refuta la Hipótesis de Riemann.

En la lectura geométrica, `pfi01.py` y `plot3d_html.py` son instrumentación:
generan datos y los dejan ver, no conclusiones. La dimensión de caja es
descriptiva y depende de la resolución de muestreo.

En la lectura espectral, con N=5000 la comparación GUE dejó de ser smoke
test: dos estadísticos de corto alcance (P(s<0.3) y ⟨r⟩) muestran un exceso
sobre el valor asintótico GUE de ~5.4σ y ~4.2σ respectivamente, y un tercero
de largo alcance (Δ₃) apunta en la misma dirección (espectro más rígido que
GUE). No lo interpretamos como anomalía: es consistente en dirección y orden
de magnitud con la corrección de altura finita a la convergencia GUE descrita
por Bogomolny–Keating (1996) y Berry–Keating (1999) — a T finito el espectro
es más rígido que el límite asintótico, y la corrección decae con T. Dos
ventanas de t (antes/después de t=1500) se mueven en esa dirección. Un
preprint reciente (Alarcón 2026, sin revisión por pares) da un ajuste
numérico específico para ⟨r⟩(T); lo mencionamos solo como referencia
tentativa, no como fuente del mecanismo.

## Historial

- v0.1 → v0.2 (`pfi01.py`): la resolución de la trayectoria pasó de un total
  fijo de puntos a una densidad fija (puntos por unidad de t), porque el
  total fijo contaminaba la comparación entre corridas con distinto
  N_ZEROS, y además el conteo de cajas se quedaba sin puntos a epsilon
  chico (contaba muestras sueltas, no la curva) — se corrigió interpolando
  entre muestras antes de contar cajas.
- `plot3d_html.py` v1.1: el panel de proyección 2D pasó de `markers` a
  `lines+markers` — sin línea, una espiral que vuelve sobre regiones
  parecidas del plano es indistinguible de ruido.
- N=1000 → N=5000: la comparación GUE pasó de smoke test a hallazgo con
  peso estadístico (ver "Alcance"). En el camino se corrigieron dos errores
  de revisión ajena sin verificar (una desviación estándar de ⟨r⟩ asumida
  en vez de medida, que inflaba la significancia a ~14σ cuando la cifra
  medida da ~4σ; y una fórmula de Δ₃ con el prefactor de GOE aplicado a
  GUE) y dos bugs propios en la primera implementación de Δ₃ (ajuste sobre
  los puntos muestreados en vez de sobre la integral continua; una división
  por L de más). Ninguno de los cuatro se detectó por inspección — los
  cuatro salieron de recalcular con los datos en mano.
- `salida_pfi01_N5000/gue_comparacion.csv` quedó con solo el encabezado en
  la primera entrega (el análisis se hizo por consola y nunca se escribió
  a archivo) — corregido, ahora trae las 4999 filas.
- v0.3 (`docs/manuscrito.md`): se agregó la sección V con los resultados
  espectrales de N=5000 y se creó `docs/` para que el ensayo viva junto al
  pipeline en vez de aparte.
# pulso_fractal_infinito
