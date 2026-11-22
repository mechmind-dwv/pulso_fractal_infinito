# analisis_avanzado

Scripts de análisis sobre las tablas de ceros de Odlyzko (v0.4). Ejecutar desde esta carpeta.

| Script | Qué hace |
|---|---|
| `ver_matrix.py` | P(s), Σ²(L) y R₂(x) de los ceros frente a GUE simulado, Poisson y Montgomery |
| `turing_alturas.py` | Conteo tipo Turing (d_n ≈ S(T)) + control de borrado + estadística por bloques de altura |
| `selberg.py` | Var S(T) y curtosis frente a ln ln T |
| `selberg_alturas.py` | Lo mismo en las ventanas de 10⁴ ceros a 10¹², 10²¹ y 10²² (mpmath) |
| `espectral_100k.py` | Histograma P(s) y ⟨r⟩ de los 100k ceros |
| `ceros_100k.txt` | Primeros 100.000 ceros (Odlyzko) |

## Material exploratorio, sin respaldo estadístico

`cazar_atractor.py`, `atractor_3d.py`, `atractor_masivo.py`, `atractor_normalizado.py`, `agujero_negro_real.py`, `red_cosmica.py`, `red_cosmica_1M.py`, `dimension_atractor.py`, `contraste_poisson.py`: visualizaciones del espacio de fases de Γ(t) y de la inmersión de Takens de los espaciados.

**No** muestran un atractor ni una "red cósmica": con controles (espaciados barajados, GUE simulado) los resultados son indistinguibles. Se conservan como instrumentación visual. Detalle en el README principal, sección "Qué NO se sostiene".

`criterio_li_1M.py` apunta a `ceros_1millon.txt`, que no está en el repo; la positividad de λₙ para n ≤ 100 es un resultado ya conocido.

## Datos

Ver la sección "Datos" del README principal para las descargas.
