# Registro de ensayos (v0.5)

Cada ensayo: pregunta, resultado frente a teoría, control, veredicto. Lo no verificado está marcado.

| Nº | Ensayo | Veredicto |
|---|---|---|
| E01 | Estadística local y de largo alcance frente a GUE | Sólido |
| E02 | Conteo de ceros (Turing) hasta 1.13e6 | Sólido, numérico |
| E03 | Selberg: varianza de S(T) y de ln\|ζ\| | Parcial |
| E04 | Momentos de \|ζ\| frente a Keating–Snaith e Ivić–Conrey | Sólido |
| E05 | Fórmula explícita: primos desde los ceros | Sólido |
| E06 | SYK frente a los ceros (clase de simetría) | Sólido en clase GUE |
| E07 | Reconstrucción de ψ(x) con K ceros | Sólido |
| E08 | Criterio de Li | Circular, sin valor probatorio |
| E09 | Atractor 3D, red cósmica, toro, materia oscura | Refutado como evidencia |
| E10 | Operador inverso (Jacobi) desde los ceros | Negativo |
| E11 | xp discretizado en una caja | Negativo |
| E12 | Positividad de Weil con gaussianas | Parcial |

## Resultados clave

- **E01**: ⟨r⟩ de los ceros a t≈1e6 = 0.6056 (GUE 0.5996). Σ²(L) satura (Berry). Control: GUE simulado.
- **E02**: N(t) coherente con 2 001 052 ceros hasta t = 1 132 490.7. Comprobación numérica tipo Turing, no prueba;
  la verificación rigurosa llega a ≈3·10¹² (Platt–Trudgian, arXiv:2004.09765).
- **E03**: Var S(T) ~ (1/2π²) ln ln T. Var ln\|ζ\|: pendiente ajustada 0.555 frente a 0.500 (abierto).
- **E04** (T = 1e8, medido/KS): k=0.5 → 1.003; k=1 → 1.010 (M₁ local = L + 2γ); k=1.5 → 1.063; k=2 → 1.310.
  El cuarto momento coincide con el polinomio de Ivić–Conrey a 0.993–1.004 hasta T = 1e6.
- **E05**: picos de la suma Σ w(tₙ)cos(x tₙ) en x = ln n solo para potencias de primo, cociente 0.9999.
  Control con 39 x genéricos: \|S\| ≤ 3.0 frente a decenas de miles.
- **E06**: SYK q=4, N=20 (GSE) ⟨r⟩ = 0.6738 ± 0.0058 (0.6744); N=22 (GUE) ⟨r⟩ = 0.6045 ± 0.0042 (0.5996);
  P(s) y Δ₃(L) de N=22 coinciden con GUE (Δ₃(18) = 0.2032 frente a 0.2055; GOE daría 0.286).
  Misma clase de simetría que los ceros; no mismo origen (SYK no tiene estructura aritmética).
- **E07**: error RMS de ψ(x) en [2, 60] con 10 / 100 / 1000 / 10⁴ / 10⁵ ceros = 0.743 / 0.373 / 0.148 / 0.071 / 0.022.
  Control con 1000 γ al azar: 2.532.
- **E10**: la matriz de Jacobi de 80 ceros no se distingue de la de un espectro GUE de igual densidad.
- **E11**: xp en una caja da niveles equiespaciados (⟨r⟩ = 1.000), no GUE, y densidad constante.
  La discretización en rejilla está dominada por efectos de la rejilla.
- **E12**: la forma cuadrática de Weil calculada solo con primos coincide con la de los ceros a 1e-10 y es
  positiva en las familias probadas (autovalor mínimo +0.133 con s=3; +3.5e-7 con s=1.5). Una familia finita
  solo comprueba una condición necesaria.

## Retractaciones y errores corregidos

- Atractor extraño, "red cósmica" y Takens 3D: retirados como evidencia.
- Criterio de Li como prueba: retirado (circular).
- Constante de GOE usada por error en Σ² de GUE: corregida.
- Media acumulada en vez de local en el chequeo de M₁: corregida.
- Ruido de K(τ) en SYK subestimado (±0.05–0.1; real ±0.3–0.5).
- "Berry–Keating exige romper T porque los ceros son GUE": es una sugerencia (Schumayer–Hutchinson,
  arXiv:1101.3116), no un teorema.

## Abiertos

- Pendiente 0.555 frente a 0.500 en Var ln\|ζ\|.
- Desviaciones de K(τ) a τ pequeño y pico en τ = 1 en los ceros (atribuido a un efecto de red, sin verificar).
- Máximo de ln\|ζ\| en ventanas frente a ln ln T − ¾ ln ln ln T (Fyodorov–Hiary–Keating).
- Operador de Hilbert–Pólya: sin candidato aceptado. Sin leer aún: Connes, billares de microondas con ferritas,
  cuerpos de funciones y la réplica arXiv:1705.06767.

## Literatura leída

- Schumayer–Hutchinson, arXiv:1101.3116: los ceros son GUE a corto alcance; sugiere un sistema que rompa T.
- Sierra, arXiv:1601.01797: con H = xp la simetría temporal se rompe; no lo deriva de GUE.
- Bender–Brody–Müller, arXiv:1608.03679: operador no hermitiano; realidad de autovalores heurística.
- Bellissard, arXiv:1704.02644: objeta autoadjunción y pertenencia a L².

## Ensayos añadidos después del registro inicial

- **E13 — Correlación de pares de SYK.** R₂(x) de SYK N=22 frente a 1 − (sin πx/πx)²: RMS 0.020 (8
  realizaciones). Control N=16 (GOE): RMS 0.089 y R₂(0.125) = 0.184 frente a 0.065. Misma clase GUE que los ceros.
- **E14 — Aritmética necesaria para localizar los ceros** (carpeta `ensayo3_enteros_y_ceros/`). Con los primeros
  M enteros de la suma de Riemann–Siegel, la fracción de ceros localizados pasa de ≈0.69 (M = N_RS/2) a 0.97–0.99
  (M = ⌊√(t/2π)⌋) a alturas 1e4 y 6e4. Un ⟨r⟩ GUE-like aparece ya con M = 20 sin que los ceros estén bien situados.
  Los picos en x = ln n están en los ceros reales (cociente ≥ 0.996 con la predicción de Weil) y no en un GUE de
  igual densidad. Una ventana por altura; no se construye ningún operador.

## Ensayos añadidos después del registro inicial

- **E13 — Correlación de pares de SYK.** R₂(x) de SYK N=22 frente a 1 − (sin πx/πx)²: RMS 0.020 (8
  realizaciones). Control N=16 (GOE): RMS 0.089 y R₂(0.125) = 0.184 frente a 0.065. Misma clase GUE que los ceros.
- **E14 — Aritmética necesaria para localizar los ceros** (carpeta `ensayo3_enteros_y_ceros/`). Con los primeros
  M enteros de la suma de Riemann–Siegel, la fracción de ceros localizados pasa de ≈0.69 (M = N_RS/2) a 0.97–0.99
  (M = ⌊√(t/2π)⌋) a alturas 1e4 y 6e4. Un ⟨r⟩ GUE-like aparece ya con M = 20 sin que los ceros estén bien situados.
  Los picos en x = ln n están en los ceros reales (cociente ≥ 0.996 con la predicción de Weil) y no en un GUE de
  igual densidad. Una ventana por altura; no se construye ningún operador.
- **E15 — GUE con primos metidos a mano** (`ensayo3_enteros_y_ceros/gue_con_primos.py`). Un GUE desplegado y
  modulado con la serie de primos truncada (P = 1…200, peso de Fejér) reproduce los picos en x = ln n con la
  amplitud de la fórmula explícita (dentro del ruido), pero ⟨r⟩ baja de 0.599 a 0.580 (ceros: 0.612) y |S| en x
  genéricos se queda en ≈ 2.6 frente a 0.01 de los ceros. Estadística local GUE más picos en los primos no
  caracterizan los ceros. Una construcción, una ventana.
- **E15 — GUE con primos metidos a mano** (`ensayo3_enteros_y_ceros/gue_con_primos.py`). Un GUE desplegado y
  modulado con la serie de primos truncada (P = 1…200, peso de Fejér) reproduce los picos en x = ln n con la
  amplitud de la fórmula explícita (dentro del ruido), pero ⟨r⟩ baja de 0.599 a 0.580 (ceros: 0.612) y |S| en x
  genéricos se queda en ≈ 2.6 frente a 0.01 de los ceros. Estadística local GUE más picos en los primos no
  caracterizan los ceros. Una construcción, una ventana.
- **E16 — El silencio como precisión de posición** (`ensayo3_enteros_y_ceros/silencio.py`). |S| en x genéricos:
  0.014 en los ceros reales; con ruido de posición ε espaciados sube a 0.06 (ε = 10⁻²), 0.51 (0.1) y 2.1 (0.4),
  mientras los picos apenas cambian. Ceros de la suma de Riemann–Siegel principal: 1.40 con M = N_RS (79) frente
  a 0.014. Un modelo que reproduzca el silencio ha de fijar cada cero a mucho menos del 1 % de un espaciado.
- **E17 — Riemann–Siegel completa** (`ensayo3_enteros_y_ceros/rs_completa.py`). Con N = ⌊√(t/2π)⌋ términos y la
  corrección C₀ se localizan todos los ceros de las ventanas T0 = 5 000, 40 000 y 70 000 con error máximo
  1.5·10⁻⁴, 4.6·10⁻⁵ y 3.3·10⁻⁵ espaciados, y el silencio |S| coincide con el de los ceros reales (0.010–0.014).
  Sin C₀: 2–12 ceros perdidos y |S| = 0.7–1.4. Propiedad conocida de ζ, cuantificada con la métrica de E16.
- **E17 — Riemann–Siegel completa** (`ensayo3_enteros_y_ceros/rs_completa.py`). Con N = ⌊√(t/2π)⌋ términos y la
  corrección C₀ se localizan todos los ceros de las ventanas T0 = 5 000, 40 000 y 70 000 con error máximo
  1.5·10⁻⁴, 4.6·10⁻⁵ y 3.3·10⁻⁵ espaciados, y el silencio |S| coincide con el de los ceros reales (0.010–0.014).
  Sin C₀: 2–12 ceros perdidos y |S| = 0.7–1.4. Propiedad conocida de ζ, cuantificada con la métrica de E16.
