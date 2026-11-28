# analisis_avanzado

Scripts de análisis sobre las tablas de ceros de Odlyzko (v0.5). Ejecutar desde esta carpeta.
Registro completo de ensayos, con resultados, controles y retractaciones: `../docs/ENSAYOS.md`.

| Script | Qué hace | Ensayo |
|---|---|---|
| `ver_matrix.py` | P(s), Σ²(L) y R₂(x) de los ceros frente a GUE simulado, Poisson y Montgomery | E01 |
| `turing_alturas.py` | Conteo tipo Turing (d_n ≈ S(T)) con control de borrado | E02 |
| `selberg.py`, `selberg_alturas.py` | Var S(T) y curtosis frente a ln ln T, hasta 10²² | E03 |
| `momentos_zeta.py` | Momentos de \|ζ\| frente a Keating–Snaith (`validar` o `Tmax`) | E04 |
| `momento4_polinomio.py` | Cuarto momento frente al polinomio de Ivić–Conrey | E04 |
| `orbitas_primos.py` | Fórmula explícita: picos en x = ln n y factor de forma K(τ) | E05 |
| `syk.py`, `stats_syk.py` | SYK q=4 (⟨r⟩, K(τ), P(s), Δ₃) frente a la clase RMT | E06 |
| `xp_caja.py` | Hamiltoniano xp en una caja: peine rígido, no GUE | E11 |
| `espectral_100k.py` | Histograma P(s) y ⟨r⟩ de los 100k ceros | E01 |
| `psi_desde_ceros.py` | Reconstrucción de ψ(x) con K ceros (fórmula explícita) y control con γ al azar | E07 |
| `psi_escalada.py`, `barrido_K.py` | (ψ(x) − x)/√x hasta 10⁷ frente a la reconstrucción y a la cota condicional a RH; RMS según K | E18 |
| `weil_positividad.py` | Matriz de Weil con gaussianas, solo primos frente a ceros, y control sintético | E12 |
| `jacobi_inverso.py` | Matriz de Jacobi (Lanczos, 250 dígitos) desde 80 ceros frente a GUE | E10 |
| `neff_cue.py` | Ceros frente a CUE de tamaño efectivo N_eff (Bogomolny et al. 2006) | E19 |
| `r2_bk.py` | R₂(s) de los ceros frente a GUE y a la primera corrección de Bogomolny–Keating; amplitud según la altura | E20 |
| `bk_completa.py` | Fórmula completa de Bogomolny–Keating para R₂(s) (diagonal + fuera de diagonal, sin expandir) frente a los ceros | E21 |
| `pico_tau1.py` | Pico K(τ = 1) ≈ 41: identidad K(1) = \|Σ e^{2πi dₙ}\|²/N con el conteo de Turing | E22 |

## Material exploratorio, sin respaldo estadístico

`cazar_atractor.py`, `atractor_3d.py`, `atractor_masivo.py`, `atractor_normalizado.py`, `agujero_negro_real.py`,
`red_cosmica.py`, `red_cosmica_1M.py`, `dimension_atractor.py`, `contraste_poisson.py`: visualizaciones del espacio
de fases de Γ(t) y de la inmersión de Takens de los espaciados.

**No** muestran un atractor ni una "red cósmica": con controles (espaciados barajados, GUE simulado) los
resultados son indistinguibles. Se conservan como instrumentación visual (E09).

`criterio_li.py` y `criterio_li_1M.py`: el criterio de Li es circular cuando se impone ρ = 1/2 + it (E08).

## Datos

No se versionan. Ver la sección "Datos" del README principal para las descargas.

### E23

| Script | Qué hace |
|---|---|
| `connes_zst.py` | Forma de Weil con primos ≤ X en doble precision; `python3 connes_zst.py 13 40 ceros_100k.txt [control]` |
| `connes_mp.py` | Misma construccion con mpmath y barrido en precision multiple; `python3 connes_mp.py 13 40 120` (~30-50 s) |

### E24 y E25

| Script | Qué hace |
|---|---|
| `matroide_li.py` | Relaciones enteras entre las ordenadas de los ceros (PSLQ); `python3 matroide_li.py 60 8 80 100000` |
| `grafo_controles.py` | Grafo cuántico K_V con longitudes ln p y controles (azar, ln n, S simétrica); `python3 grafo_controles.py 6 60 40` |

`connes_mp.py` cuenta ahora los ceros consecutivos reproducidos y reparte los intervalos con polo (arregla ceros perdidos junto a un polo 2πn/L). Con X = 50 y N = 95: `python3 connes_mp.py 50 95 160` (~10 min en móvil).

### E26

| Script | Qué hace |
|---|---|
| `ps_bbl.py` | P(s) de los ceros frente a GUE (núcleo seno) y a CUE de N finito (determinante de Fredholm por Nyström); `python3 ps_bbl.py ceros_100k.txt` |
| `var_alturas.py` | var(s) y ⟨r⟩ por bloques de altura y ajuste del exceso; `python3 var_alturas.py ceros_1millon.txt 100000` |
| `var_alto.py` | var(s) y ⟨r⟩ de las ventanas de Odlyzko (10¹², 10²¹, 10²²) frente a los modelos de E26; `python3 var_alto.py .` |
