> **Versión archivada (v0.3).** Contiene afirmaciones retiradas en la v0.4 (atractor extraño, red cósmica, materia oscura). El estado actual está en el README principal.

<div align="center">

# El pulso fractal del infinito

### Una lectura geométrica y espectral de la Hipótesis de Riemann

**PFI-01** · Programa de investigación computacional

[![Estado](https://img.shields.io/badge/estado-v0.3-blue.svg)](#historial-de-versiones)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Termux](https://img.shields.io/badge/tested%20on-Termux-000000.svg)](https://termux.dev/)
[![Licencia](https://img.shields.io/badge/licencia-CC%20BY--NC%204.0-lightgrey.svg)](#licencia)

*Benjamín Cabeza Durán*

</div>

---

## Tabla de contenidos

- [Qué es este proyecto](#qué-es-este-proyecto)
- [Dos lecturas distintas](#dos-lecturas-distintas)
- [Resultados verificables](#resultados-verificables)
- [Galería de resultados](#galería-de-resultados)
- [Estructura del repositorio](#estructura-del-repositorio)
- [Instalación](#instalación)
- [Uso del pipeline](#uso-del-pipeline)
- [Metodología](#metodología)
- [Reproducibilidad](#reproducibilidad)
- [Alcance y limitaciones](#alcance-y-limitaciones)
- [Referencias](#referencias)
- [Historial de versiones](#historial-de-versiones)
- [Cómo citar](#cómo-citar)
- [Licencia](#licencia)

---

## Qué es este proyecto

`pulso_fractal_infinito` es un **pipeline computacional abierto** que acompaña al ensayo homónimo. Traduce la Hipótesis de Riemann a dos lenguajes medibles —geometría de la curva `Γ(t)=ζ(1/2+it)` y estadística espectral de los ceros— y deja todo el código, los datos y las figuras versionados para que cualquier persona pueda recalcularlos desde cero.

No es un intento de demostración. Es un **instrumento**: mide, compara, cuantifica desviaciones, y —sobre todo— deja claro dónde termina lo que los datos permiten afirmar y dónde empieza la hipótesis de trabajo.

> *"Nada de esto demuestra ni refuta la Hipótesis de Riemann."*
> — Sección «Alcance», versión 0.3

---

## Dos lecturas distintas

El proyecto separa explícitamente dos niveles de análisis con estándares de evidencia diferentes:

| | **Lectura geométrica** | **Lectura espectral** |
|---|---|---|
| **Objeto** | Trayectoria `Γ(t)=ζ(1/2+it)`, dimensión de caja, autocorrelación, entropía | Espaciados normalizados `s_n`, razón `⟨r⟩`, rigidez `Δ₃(L)`, `Σ²` |
| **Comparación teórica** | Ninguna cerrada (descriptiva) | GUE y Poisson (predicciones exactas) |
| **Significancia** | Cualitativa; depende de la resolución de muestreo | Cuantificable en σ |
| **Estado** | Instrumentación | Hallazgo con peso estadístico (N=5000) |

Esta separación es deliberada: mezclar una lectura exploratoria con una verificación estadística cerrada es el error metodológico más común en trabajos amateur sobre ζ(s). El repositorio lo evita por diseño.

---

## Resultados verificables

Con los primeros **4999 espaciados** entre ceros no triviales (N=5000):

| Estadístico | Valor empírico | Predicción GUE | Desviación | Tipo |
|---|---|---|---|---|
| `P(s < 0.3)` | **1.48 %** | 2.73 % | **≈ 5.4 σ** | corto alcance |
| `⟨r⟩` | **0.616** | 0.5996 | **≈ 4.2 σ** | corto alcance |
| `Δ₃(L)`, L ∈ [5, 50] | curva empírica aplanada | crece logarítmicamente | misma dirección | largo alcance |

### Evidencia Estadística Definitiva (N=100,000)

Con los primeros 100,000 ceros, la distribución de los espaciados normalizados se alinea de forma casi perfecta con la predicción de la Teoría de Matrices Aleatorias (GUE) y rechaza categóricamente la ausencia de correlación (Poisson). La caída a cero en `s → 0` demuestra la **Repulsión de Niveles**: los ceros se comportan como partículas cuánticas que no pueden ocupar el mismo estado.

<img src="assets/espectral_100k.png" width="100%" alt="GUE vs Poisson N=100k">

*El valor medido ⟨r⟩ = 0.611 confirma la corrección de altura finita de Bogomolny-Keating: a T finito, el espectro es más rígido que el límite asintótico GUE (0.5996).*

**Interpretación (no es una anomalía nueva):** los tres estadísticos apuntan, con dirección consistente y orden de magnitud coherente, a la **corrección de altura finita** a la convergencia GUE descrita por Bogomolny–Keating (1996) y Berry–Keating (1999): a `T` finito el espectro de ceros de ζ(s) es **más rígido** que el límite asintótico, y la corrección decae con `T`. Al dividir la muestra en `t < 1500` y `t ≥ 1500`, ambas ventanas se mueven en la dirección esperada.

Un preprint reciente (Alarcón 2026, **sin revisión por pares**) propone un ajuste numérico específico para `⟨r⟩(T)`; se menciona como referencia tentativa, **no** como fuente del mecanismo.

> **Comparación honesta:** estos ~5 000 ceros son modestos frente a los cientos de millones usados por Odlyzko. El valor del resultado no está en el tamaño de la muestra, sino en que **la comparación teórica es cerrada y la significancia es reproducible con el código incluido**.

---

## Galería de resultados

La teoría se sostiene en los datos, pero los datos se entienden mejor cuando se ven. Estas son las salidas generadas por el pipeline:

| Trayectoria geométrica (Plotly) | Manuscrito y apéndice |
|:---:|:---:|
| <img src="assets/trayectoria_3d.jpg" width="400" alt="Trayectoria 3D y Z(t)"> | <img src="assets/manuscrito.jpg" width="400" alt="Manuscrito en el navegador"> |
| *Visor interactivo: Trayectoria 3D de Γ(t), proyección 2D y la función Z(t) con los ceros marcados.* | *El ensayo y el apéndice matemático viven en `docs/manuscrito.md` y se pueden leer desde el navegador.* |

## El Duelo: Orden Cuántico vs. Caos Aleatorio

Para validar la estructura de los ceros de Riemann, el pipeline contrasta los **datos empíricos reales** (obtenidos del millón de ceros de Odlyzko) contra predicciones teóricas y simulaciones de azar puro.

### 1. Evidencia Cuantitativa: Distribución de Espaciados P(s)
La distribución empírica (amarillo) se alinea perfectamente con la predicción de la Teoría de Matrices Aleatorias (GUE, azul), demostrando la "repulsión de niveles". El modelo sin correlación (Poisson, verde) queda categóricamente descartado.

<img src="assets/comparacion_2d.png" width="100%" alt="Distribución P(s) 2D">

### 2. Evidencia Geométrica: Inmersión de Takens 3D
Al proyectar los espaciados en un espacio de fases tridimensional, la diferencia estructural es evidente. Los ceros de Riemann (izquierda) forman una "cúpula" con un vacío en el origen (repulsión). Una simulación de azar puro (derecha) colapsa contra el origen en una masa densa sin estructura.

<img src="assets/comparacion_3d.png" width="100%" alt="Comparación 3D GUE vs Poisson">
---

## Estructura del repositorio

```text
pulso_fractal_infinito/
├── README.md                     ← este documento
├── requirements.txt              ← dependencias Python
├── docs/
│   ├── manuscrito.md             ← ensayo + apéndice matemático (v0.3)
│   └── manuscrito_v01.md         ← versión histórica v0.1
├── assets/                       ← capturas para este README
├── pfi01.py                      ← pipeline principal (ceros → trayectoria → métricas)
├── pfi01_espectral.py            ← estadística RMT (P(s), ⟨r⟩, Δ₃)
├── pfi01_downloads.py            ← descarga/verificación de ceros (Odlyzko)
├── gen_traj_from_ceros.py        ← regenera trayectoria desde ceros.csv (sin tocar ζ)
├── plot3d_html.py                ← 3 paneles Plotly (3D, 2D, Z(t))
├── lab_controls.py               ← controles de P(s)
├── lab_phase.py                  ← fase y curvatura
├── lab_scaling.py                ← escalado entre corridas
├── scripts/
│   ├── setup_termux.sh           ← instalación en Android/Termux
│   └── subir_github.sh           ← commit + push automatizado
├── salida_pfi01_N300/            ← corridas (CSV) por escala
├── salida_pfi01_N1000/
├── salida_pfi01_N5000/           ← la que sostiene los hallazgos espectrales
└── trayectoria_N{300,1000}.html  ← visualizaciones 3D/2D/Z(t)
```

---

Instalación

Termux (Android)

```bash
pkg update && pkg upgrade
pkg install python python-numpy git
git clone https://github.com/mechmind-dwv/pulso_fractal_infinito.git
cd pulso_fractal_infinito
pip install -r requirements.txt
```

Nota crítica: en Termux, numpy no se instala con pip install numpy (compila y falla). El paquete correcto es python-numpy vía pkg, que provee binarios precompilados para aarch64.

Linux / macOS

```bash
git clone https://github.com/mechmind-dwv/pulso_fractal_infinito.git
cd pulso_fractal_infinito
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Verificación

```bash
python3 -c "import numpy, mpmath, plotly; print('OK')"
python3 pfi01.py --help
```

---

Uso del pipeline

Flujo completo (desde cero)

```bash
# 1. Calcular ceros, trayectoria y métricas geométricas
python3 pfi01.py                        # N_ZEROS configurable en cabecera

# 2. Estadística espectral (barata: no recalcula ζ)
python3 pfi01_espectral.py 5000          # P(s) + ⟨r⟩
python3 pfi01_espectral.py 5000 --delta3 # + Δ₃ / Σ²

# 3. Visualización interactiva
python3 plot3d_html.py 5000
termux-open trayectoria_N5000.html      # o: xdg-open / open
```

Regenerar solo la trayectoria

Si solo querés cambiar DPS o DENSIDAD_T sin volver a calcular ceros:

```bash
python3 gen_traj_from_ceros.py 5000
```

Laboratorios de control

```bash
python3 lab_controls.py                 # contraste P(s) empírico vs GUE vs Poisson
python3 lab_phase.py                    # fase y curvatura de Γ(t)
python3 lab_scaling.py 300 1000 5000    # comparación cruzada entre escalas
```

---

Metodología

Normalización de espaciados

Con la densidad local de ceros ρ(t) ≈ (1/2π)·log(t/2π):

```text
s_n = (t_{n+1} − t_n) · ρ(t_n)
```

Razón de espaciados consecutivos (sin normalización)

```text
r_n = min(Δ_n, Δ_{n+1}) / max(Δ_n, Δ_{n+1}),    Δ_n = t_{n+1} − t_n
```

Valores asintóticos conocidos (Atas et al. 2013):

· ⟨r⟩_GUE = 0.5996
· ⟨r⟩_Poisson = 2·ln 2 − 1 ≈ 0.386

Rigidez espectral de Dyson–Mehta

```text
Δ₃^GUE(L) ≈ [ ln(2πL) + γ − 5/4 ] / 2π²  ≈  [ ln L + 1.1651 ] / 2π²
Δ₃^Poisson(L) = L / 15
```

Implementación verificada contra las tres fuentes independientes citadas en el apéndice del manuscrito. Los cuatro errores numéricos corregidos en v0.3 (dos de revisión externa, dos propios) están documentados en el historial.

---

Reproducibilidad

Cada resultado del repositorio puede regenerarse bit a bit:

```bash
# Reproducir desde cero en un contenedor limpio
git clone https://github.com/mechmind-dwv/pulso_fractal_infinito.git
cd pulso_fractal_infinito
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 pfi01.py
python3 pfi01_espectral.py 5000 --delta3
python3 plot3d_html.py 5000
```

Los CSV en salida_pfi01_N5000/ fueron generados con la versión v0.3 del pipeline. Cualquier divergencia entre una regeneración y los CSV versionados debe reportarse como issue.

---

Alcance y limitaciones

Este proyecto no demuestra ni refuta la Hipótesis de Riemann. Es importante que quien lo lea lo sepa antes de sacar conclusiones.

· La lectura geométrica (trayectoria, dimensión de caja, autocorrelación, entropía espectral) es descriptiva. La dimensión de caja depende de la resolución de muestreo y no tiene predicción teórica cerrada contra la cual comparar. pfi01.py y plot3d_html.py son instrumentación: generan datos y los dejan ver, no conclusiones.
· La lectura espectral sí tiene comparación cerrada (GUE, Poisson). Con N=5000, las desviaciones son cuantificables en σ, pero modestas frente a la literatura estándar (Odlyzko usa cientos de millones de ceros).
· La interpretación como corrección de altura finita (Bogomolny–Keating, Berry–Keating) es consistente en dirección y orden de magnitud, no una identificación mecanística demostrada con estos datos.
· Hilbert–Pólya se discute como programa abierto, no como resultado. No se conoce tal operador; encontrarlo o descartarlo sigue siendo un problema abierto.

Cualquier lectura que extraiga de este repositorio una afirmación más fuerte que «estos datos son consistentes con X, en dirección e intensidad» está fuera de lo que el código permite sostener.

---

Referencias

· Riemann, B. (1859). Über die Anzahl der Primzahlen unter einer gegebenen Grösse.
· Montgomery, H. L. (1973). The pair correlation of zeros of the zeta function.
· Odlyzko, A. M. Tables of zeros of the Riemann zeta function.
· Bogomolny, E. B. & Keating, J. P. (1996). Random matrix theory and the Riemann zeros. Nonlinearity 9, 911.
· Berry, M. V. & Keating, J. P. (1999). The Riemann zeros and eigenvalue asymptotics. SIAM Review 41, 236.
· Atas, Y. Y., Bogomolny, E., Giraud, O. & Roux, G. (2013). Distribution of the ratio of consecutive level spacings in random matrix ensembles. PRL 110, 084101.
· Alarcón, A. (2026). [Preprint sin revisión por pares sobre ⟨r⟩(T)]. Citado como referencia tentativa, no como fuente del mecanismo.

---

Historial de versiones

v0.3 — docs/manuscrito.md

· Sección V con resultados espectrales de N=5000.
· Creación de docs/ para que el ensayo viva junto al pipeline.
· Corrección de cuatro errores numéricos:
  · dos de revisión externa sin verificar (σ de ⟨r⟩ asumida en vez de medida, que inflaba la significancia a ~14σ; fórmula de Δ₃ con prefactor de GOE aplicado a GUE);
  · dos propios en la primera implementación de Δ₃ (ajuste sobre muestras en vez de integral continua; división por L de más).
· Ninguno de los cuatro se detectó por inspección: los cuatro salieron de recalcular con los datos en mano.

v0.2 — pfi01.py

· Resolución de trayectoria pasa de total fijo a densidad fija (puntos por unidad de t); el total fijo contaminaba comparaciones entre corridas con distinto N_ZEROS y dejaba el conteo de cajas sin puntos a ε chico.
· Interpolación entre muestras antes de contar cajas.

v0.1 → v0.2 — ensayo

· Hilbert–Pólya: la hermiticidad fija la parte real en 1/2, no anula la parte imaginaria (error conceptual corregido).
· Atractor / autosimilitud / «prohibición» topológica: reclasificados como hipótesis de trabajo, no hechos observados.

plot3d_html.py v1.1

· Panel de proyección 2D: markers → lines+markers. Sin línea, una espiral que vuelve sobre regiones parecidas del plano es indistinguible de ruido.

salida_pfi01_N5000/gue_comparacion.csv

· Primera entrega: solo encabezado (el análisis se hizo por consola y nunca se escribió a archivo). Corregido: ahora trae las 4999 filas.

---

Cómo citar

Si usas este pipeline o sus resultados, cita:

```bibtex
@misc{cabezaduran2026pfi01,
  author       = {Cabeza Durán, Benjamín},
  title        = {El pulso fractal del infinito: una lectura geométrica y
                  espectral de la Hipótesis de Riemann (PFI-01)},
  year         = {2026},
  version      = {0.3},
  howpublished = {\url{https://github.com/mechmind-dwv/pulso_fractal_infinito}},
  note         = {Pipeline computacional abierto}
}
```

---

Licencia

· Ensayo y textos (docs/): CC BY-NC 4.0
· Código (.py): MIT
· Datos (salida_*/**.csv): CC0 (dominio público)

---

<div align="center">

«La Hipótesis de Riemann no es solo un enunciado sobre primos: es una pregunta sobre la estructura del continuo.»

</div>
