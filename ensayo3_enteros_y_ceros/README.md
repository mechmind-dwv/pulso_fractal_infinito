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
