# Ensayo 3 — ¿Cuánta aritmética hace falta para obtener los ceros?

Registro: E14 en `../docs/ENSAYOS.md`. Scripts solo con numpy. Ejecutar con `bash run_all.sh`.

## Pregunta

Un operador de Hilbert–Pólya tendría que producir a la vez (i) estadística GUE y (ii) los primos como
órbitas periódicas. Buscar una matriz con esas dos propiedades es investigación abierta y no la intentamos aquí.
Lo que sí se puede medir con los datos de que disponemos es una pregunta previa y más modesta:

- **A.** ¿Cuántos enteros n ≤ M hacen falta para localizar los ceros a una altura t, y qué estadística tiene el
  conjunto de ceros que se obtiene si se trunca antes?
- **B.** ¿Se ven los primos en los ceros reales y no en un espectro GUE genérico de la misma densidad?

## Método

**Parte A** (`ceros_desde_enteros.py`). Se usa la suma principal de Riemann–Siegel truncada:
Z_M(t) = 2 · Σ_{n≤M} cos(θ(t) − t ln n)/√n. Para M = ⌊√(t/2π)⌋ (el corte de Riemann–Siegel) es la aproximación
estándar; para M menor faltan términos. Se buscan los cambios de signo en una ventana de 1500 de ancho y se
comparan con los ceros reales: fracción encontrada (a menos de 0.25 espaciados), error rms en espaciados y ⟨r⟩.
M = 1 da Z = 2cos θ(t), es decir, la red regular de la densidad media (⟨r⟩ = 1).

**Parte B** (`contraste_trazas.py`). Se calcula S(x) = Σ w(t)cos(x t) con ventana gaussiana en torno a T0.
La fórmula explícita predice, en x = ln n, S = −Λ(n)/(2π√n) · ∫₀^∞ w dt. Se compara con un espectro GUE
(matriz 2·K×2·K) desplegado y llevado a la misma densidad media con la inversa de N_suave(t).

## Resultados (ceros_100k.txt)

Parte A, ventana [10 000, 11 500], N_RS = 41, 1779 ceros reales, ⟨r⟩ real = 0.6147:

| M | ceros encontrados | error rms / espaciado | ⟨r⟩ |
|---|---|---|---|
| 1 | 0.589 | 0.260 | 1.0000 |
| 5 | 0.629 | 0.338 | 0.8087 |
| 20 | 0.687 | 0.448 | 0.6165 |
| 41 (= N_RS) | **0.968** | **0.176** | 0.6008 |
| 82 (2·N_RS) | 0.595 | 0.569 | 0.4650 |

Ventana [60 000, 61 500], N_RS = 98, ⟨r⟩ real = 0.6106: con M = 98 se encuentra el 0.988 de los ceros,
error rms 0.107 y ⟨r⟩ = 0.6102.

Parte B, T0 = 40 000, σ = 100, 837 ceros (primeros primos):

| n | predicción de Weil | ceros reales | GUE |
|---|---|---|---|
| 2 | −19.55 | −19.51 | −1.63 |
| 3 | −25.30 | −25.24 | −0.79 |
| 5 | −28.71 | −28.66 | −1.71 |
| 7 | −29.34 | −29.28 | −2.22 |
| 11 | −28.84 | −28.77 | −4.06 |

Control con 23 valores genéricos de x: |S| típico 0.01 (ceros) frente a 2.17 y máximo 10.7 (GUE).

## Lectura

1. La localización de los ceros salta de ≈0.69 a ≈0.97–0.99 al pasar de M = 20 a M = N_RS. Hace falta toda la
   aritmética hasta √(t/2π); con la mitad no basta.
2. **Un ⟨r⟩ cercano a GUE no prueba que los ceros estén bien situados.** Con M = 20 el ⟨r⟩ (0.6165) ya es GUE-like
   pero solo se encuentra el 69 % de los ceros. El estadístico local no distingue un conjunto correcto de uno
   plausible.
3. Por encima de N_RS la suma sin términos de corrección empeora (⟨r⟩ cae a 0.47). No se debe leer como que
   "más primos es peor": es la aproximación truncada sin corrección, no la función.
4. Los picos en x = ln n están en los ceros reales (cociente con la predicción ≥ 0.996) y no en un GUE de igual
   densidad. Esa es la diferencia estructural entre los ceros y cualquier ensemble genérico, SYK incluido.

## Limitaciones

- Una ventana por altura en la parte A; los números exactos dependen de la ventana.
- La parte B reproduce la fórmula explícita (ya conocida); lo informativo es el contraste con el GUE.
- El GUE de comparación es una sola realización, con desplegado polinómico.
- No se construye ningún operador. La relación entre el corte √(t/2π) y el tiempo de Heisenberg de Berry–Keating
  es una hipótesis de trabajo; aquí no se mide.

## Siguiente paso

Una matriz cuyo espectro reproduzca a la vez las trazas Σ ln p/√p y la estadística GUE sigue abierta. Una
versión medible: ajustar los parámetros de un ensemble para que S(x) tenga picos en ln p y comprobar si ⟨r⟩ y Δ₃
se mantienen.

## Parte C — GUE con los primos metidos a mano (E15)

`gue_con_primos.py`. Pregunta: si a un GUE se le impone la señal de los primos, ¿se obtienen los ceros?
Se desplegan niveles GUE (u) y se llevan a t = N_mod⁻¹(u), con N_mod(t) = N_suave(t) + S_P(t) y
S_P(t) = −(1/π) Σ_{n≤P} w_n Λ(n)/(√n ln n) sin(t ln n), con peso de Fejér w_n = 1 − ln n/ln(P+1).
Mismas 8 realizaciones GUE en todas las filas (T0 = 40 000, σ = 100, 837 niveles). Los ceros reales dan
⟨r⟩ = 0.6122, picos −19.51, −25.24, −28.66, −29.28 (n = 2, 3, 5, 7) y |S| genérico 0.01.

| P | ⟨r⟩ | S(2) | S(3) | S(5) | S(7) | \|S\| genérico |
|---|---|---|---|---|---|---|
| 1 (GUE puro) | 0.5991 ± 0.0032 | −0.11 | 0.50 | 2.03 | −0.64 | 2.62 |
| 13 | 0.5965 ± 0.0031 | −14.59 | −14.72 | −9.50 | −8.56 | 2.62 |
| 50 | 0.5906 ± 0.0029 | −16.30 | −18.17 | −15.35 | −15.81 | 2.70 |
| 200 | 0.5800 ± 0.0026 | −17.18 | −20.00 | −18.46 | −19.74 | 2.68 |

Lectura: (1) los picos aparecen con la amplitud de la fórmula explícita (con el peso de Fejér, dentro de ±2 de
ruido). (2) A medida que se meten más primos, ⟨r⟩ **baja** de 0.599 a 0.580 en vez de acercarse al 0.612 de los
ceros: la modulación lisa no es la manera en que los primos entran en las fluctuaciones. (3) El nivel de |S| en x
genéricos no cambia (≈ 2.6), frente a 0.01 en los ceros: el factor ~260 es el "silencio" entre los picos, que en los
ceros es una identidad exacta y en un GUE con primos añadidos no aparece. Estadística local GUE más picos en los
primos no bastan para caracterizar a los ceros. Limitaciones: una construcción concreta (modulación de la densidad
con peso de Fejér), una ventana, un tamaño de GUE.

## Parte C — GUE con los primos metidos a mano (E15)

`gue_con_primos.py`. Pregunta: si a un GUE se le impone la señal de los primos, ¿se obtienen los ceros?
Se desplegan niveles GUE (u) y se llevan a t = N_mod⁻¹(u), con N_mod(t) = N_suave(t) + S_P(t) y
S_P(t) = −(1/π) Σ_{n≤P} w_n Λ(n)/(√n ln n) sin(t ln n), con peso de Fejér w_n = 1 − ln n/ln(P+1).
Mismas 8 realizaciones GUE en todas las filas (T0 = 40 000, σ = 100, 837 niveles). Los ceros reales dan
⟨r⟩ = 0.6122, picos −19.51, −25.24, −28.66, −29.28 (n = 2, 3, 5, 7) y |S| genérico 0.01.

| P | ⟨r⟩ | S(2) | S(3) | S(5) | S(7) | \|S\| genérico |
|---|---|---|---|---|---|---|
| 1 (GUE puro) | 0.5991 ± 0.0032 | −0.11 | 0.50 | 2.03 | −0.64 | 2.62 |
| 13 | 0.5965 ± 0.0031 | −14.59 | −14.72 | −9.50 | −8.56 | 2.62 |
| 50 | 0.5906 ± 0.0029 | −16.30 | −18.17 | −15.35 | −15.81 | 2.70 |
| 200 | 0.5800 ± 0.0026 | −17.18 | −20.00 | −18.46 | −19.74 | 2.68 |

Lectura: (1) los picos aparecen con la amplitud de la fórmula explícita (con el peso de Fejér, dentro de ±2 de
ruido). (2) A medida que se meten más primos, ⟨r⟩ **baja** de 0.599 a 0.580 en vez de acercarse al 0.612 de los
ceros: la modulación lisa no es la manera en que los primos entran en las fluctuaciones. (3) El nivel de |S| en x
genéricos no cambia (≈ 2.6), frente a 0.01 en los ceros: el factor ~260 es el "silencio" entre los picos, que en los
ceros es una identidad exacta y en un GUE con primos añadidos no aparece. Estadística local GUE más picos en los
primos no bastan para caracterizar a los ceros. Limitaciones: una construcción concreta (modulación de la densidad
con peso de Fejér), una ventana, un tamaño de GUE.
