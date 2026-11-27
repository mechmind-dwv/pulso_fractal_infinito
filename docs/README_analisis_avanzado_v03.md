> **Versión archivada (v0.3).** Contiene afirmaciones retiradas en la v0.4 (atractor extraño, red cósmica, materia oscura). El estado actual está en el README principal.

<div align="center">

# El pulso fractal del infinito

### Una lectura geométrica y espectral de la Hipótesis de Riemann

**PFI-01** · Programa de investigación computacional

[![Estado](https://img.shields.io/badge/estado-v0.4-blue.svg)](#historial-de-versiones)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Termux](https://img.shields.io/badge/tested%20on-Termux-000000.svg)](https://termux.dev/)
[![Licencia](https://img.shields.io/badge/licencia-CC%20BY--NC%204.0-lightgrey.svg)](#licencia)

*Benjamín Cabeza Durán*

</div>

---

## Qué es este proyecto

`pulso_fractal_infinito` es un **pipeline computacional abierto** que traduce la Hipótesis de Riemann a dos lenguajes medibles: la geometría de la curva `Γ(t)=ζ(1/2+it)` y la estadística espectral de los ceros. El objetivo no es demostrar la hipótesis, sino **cartografiar la estructura oculta del caos** que gobierna los números primos.

> *"Nada de esto demuestra ni refuta la Hipótesis de Riemann."*

---

## Resultados verificables y Análisis del Caos (N=100,000)

Además de la estadística espectral estándar (GUE vs Poisson), el proyecto ha explorado la **dinámica de los ceros** mediante el análisis de atractores y mapas de retorno:

| Métrica | Resultado | Significado Físico |
|---|---|---|
| **Dimensión Fractal (D₂)** | `≈ 1.2256` | Evidencia de un **Atractor Extraño** en la sección de Poincaré. |
| **Espaciados P(s<0.3)** | `1.48 %` | Desviación de ~5.4σ sobre GUE. Espectro más rígido. |
| **⟨r⟩** | `0.616` | Desviación de ~4.2σ. Firma de caos cuántico. |

### Galería de Atractores (100k ceros de Odlyzko)

La proyección de los espaciados normalizados revela estructuras idénticas a las simulaciones de materia oscura y agujeros negros:

| Atractor (Mapa de Retorno) | Red Cósmica (Inmersión de Takens) |
|:---:|:---:|
| <img src="assets/atractor_agujero_negro.jpg" width="400" alt="Agujero Negro"> | <img src="assets/red_cosmica.jpg" width="400" alt="Red Cosmica"> |
| *Repulsión de niveles. Los ceros actúan como partículas cuánticas.* | *Filamentos y vacíos estadísticos. La huella de la materia oscura.* |

---

## Estructura del repositorio

```text
pulso_fractal_infinito/
├── README.md                     ← este documento
├── requirements.txt              ← dependencias Python
├── docs/
│   ├── manuscrito.md             ← ensayo + apéndice matemático (v0.4)
│   └── manuscrito_v01.md         ← versión histórica v0.1
├── pfi01.py                      ← pipeline principal (ceros → trayectoria → métricas)
├── pfi01_espectral.py            ← estadística RMT (P(s), ⟨r⟩, Δ₃)
├── pfi01_downloads.py            ← descarga/verificación de ceros (Odlyzko)
├── gen_traj_from_ceros.py        ← regenera trayectoria desde ceros.csv
├── plot3d_html.py                ← 3 paneles Plotly (3D, 2D, Z(t))
├── analisis_avanzado/            ← NUEVOS EXPERIMENTOS
│   ├── ceros_100k.txt            ← Datos crudos de Odlyzko
│   ├── cazar_atractor.py         ← Sección de Poincaré 2D
│   ├── atractor_3d.py            ← Espacio de fases 3D
│   ├── dimension_atractor.py     ← Cálculo de D₂
│   ├── agujero_negro_real.py     ← Mapa de retorno 100k ceros
│   ├── atractor_normalizado.py   ← Mapas GUE
│   └── red_cosmica.py            ← Inmersión de Takens (Red Cósmica)
├── scripts/
│   ├── setup_termux.sh           ← instalación en Android/Termux
│   └── subir_github.sh           ← commit + push automatizado
├── salida_pfi01_N300/            ← corridas (CSV) por escala
├── salida_pfi01_N1000/
├── salida_pfi01_N5000/
└── trayectoria_N{300,1000}.html  ← visualizaciones 3D/2D/Z(t)
```

---

Instalación y Uso Rápido

Termux (Android)

```bash
pkg update && pkg upgrade
pkg install python python-numpy python-pandas git
pip install -r requirements.txt
```

Ejecutar el análisis de atractores

```bash
cd analisis_avanzado
python3 red_cosmica.py
python3 dimension_atractor.py
```

---

Alcance y limitaciones

· Lectura geométrica: Descriptiva. No predice la Hipótesis, pero cartografía su dinámica.
· Lectura espectral: Comparación cerrada contra GUE y Poisson. Con N=5000/100k es estadísticamente significativa.
· Atractores y Red Cósmica: El cálculo de la dimensión D_2 demuestra que no es ruido blanco, sino un caos determinista. Las analogías visuales con la cosmología son una herramienta pedagógica, no una afirmación física de que los ceros sean galaxias.

---

Cómo citar

```bibtex
@misc{cabezaduran2026pfi01,
  author       = {Cabeza Durán, Benjamín},
  title        = {El pulso fractal del infinito (PFI-01)},
  year         = {2026},
  version      = {0.4},
  howpublished = {\url{https://github.com/mechmind-dwv/pulso_fractal_infinito}},
}
```

---

Licencia

· Ensayo y textos (docs/): CC BY-NC 4.0
· Código (.py): MIT
· Datos (salida_*/**.csv): CC0 (dominio público)

<div align="center">
*«La Hipótesis de Riemann no es solo un enunciado sobre primos: es una pregunta sobre la estructura del continuo.»*
</div>
