# Registro de ensayos (v0.10)

Cada ensayo: pregunta, resultado frente a teoría, control, veredicto. Lo no verificado está marcado. El registro más detallado (con los números de cada ejecución) está en el cuaderno de ensayos; aquí van los resultados que sostiene el repositorio.

| Nº | Ensayo | Veredicto |
|---|---|---|
| E01 | Estadística local y de largo alcance frente a GUE | Sólido |
| E02 | Conteo de ceros (Turing) hasta 1.13e6 | Sólido, numérico |
| E03 | Selberg: varianza de S(T) y de ln\|ζ\| | Parcial |
| E04 | Momentos de \|ζ\| frente a Keating–Snaith e Ivić–Conrey | Sólido |
| E05 | Fórmula explícita: primos desde los ceros | Sólido |
| E06 | SYK frente a los ceros (clase de simetría) | Sólido en clase GUE; no mismo origen |
| E07 | Reconstrucción de ψ(x) con K ceros | Sólido |
| E08 | Criterio de Li | Circular, sin valor probatorio |
| E09 | Atractor 3D, red cósmica, toro, materia oscura | Refutado como evidencia |
| E10 | Operador inverso (Jacobi) desde los ceros | Negativo |
| E11 | xp discretizado en una caja | Negativo |
| E12 | Positividad de Weil con gaussianas | Parcial |
| E13 | Correlación de pares de SYK | Sólido en clase GUE; no mismo origen |
| E14 | Aritmética necesaria para localizar los ceros | Positivo |
| E15 | GUE con primos metidos a mano | Negativo |
| E16 | El silencio como precisión de posición | Positivo (métrica) |
| E17 | Riemann–Siegel completa | Positivo (propiedad conocida) |
| E18 | Oscilación de ψ(x) − x con 100 ceros | Positivo |
| E19 | Ceros frente a CUE de tamaño efectivo | Parcial |
| E20 | Corrección de Bogomolny–Keating a R₂(s), primer orden | Negativo/parcial |
| E21 | Fórmula completa de Bogomolny–Keating para R₂(s) | Positivo (resultado conocido, reproducido) |
| E22 | Pico K(τ = 1) | Positivo (identidad exacta) |
| E23 | Los ceros desde los primos (forma de Weil, Connes–Consani–Moscovici) | Positivo; reproduce las cifras del artículo |
| E24 | Matroide de dependencia racional de los ceros | Nulo y esperado |
| E25 | Grafo cuántico con longitudes ln p | Nulo y esperado (GUE por simetría) |
| E26 | P(s) de los ceros a t bajo frente a GUE y CUE finito | Parcial |

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
- **E13** (R₂ de SYK): R₂(x) de SYK N=22 frente a 1 − (sin πx/πx)²: RMS 0.020 (8 realizaciones). Control N=16 (GOE):
  RMS 0.089 y R₂(0.125) = 0.184 frente a 0.065. Misma clase GUE que los ceros.
- **E14** (`ensayo3_enteros_y_ceros/`): con los primeros M enteros de la suma de Riemann–Siegel, la fracción de ceros
  localizados pasa de ≈0.69 (M = N_RS/2) a 0.97–0.99 (M = ⌊√(t/2π)⌋) a alturas 1e4 y 6e4. Un ⟨r⟩ GUE-like aparece ya
  con M = 20 sin que los ceros estén bien situados. Los picos en x = ln n están en los ceros reales (cociente ≥ 0.996
  con la predicción de Weil) y no en un GUE de igual densidad. Una ventana por altura; no se construye ningún operador.
- **E15** (`gue_con_primos.py`): un GUE desplegado y modulado con la serie de primos truncada (P = 1…200, peso de
  Fejér) reproduce los picos en x = ln n con la amplitud de la fórmula explícita (dentro del ruido), pero ⟨r⟩ baja de
  0.599 a 0.580 (ceros: 0.612) y \|S\| en x genéricos se queda en ≈ 2.6 frente a 0.01 de los ceros. Estadística local
  GUE más picos en los primos no caracterizan los ceros. Una construcción, una ventana.
- **E16** (`silencio.py`): \|S\| en x genéricos: 0.014 en los ceros reales; con ruido de posición ε espaciados sube a
  0.06 (ε = 10⁻²), 0.51 (0.1) y 2.1 (0.4), mientras los picos apenas cambian. Ceros de la suma de Riemann–Siegel
  principal: 1.40 con M = N_RS (79) frente a 0.014. Un modelo que reproduzca el silencio ha de fijar cada cero a mucho
  menos del 1 % de un espaciado.
- **E17** (`rs_completa.py`): con N = ⌊√(t/2π)⌋ términos y la corrección C₀ se localizan todos los ceros de las
  ventanas T0 = 5 000, 40 000 y 70 000 con error máximo 1.5·10⁻⁴, 4.6·10⁻⁵ y 3.3·10⁻⁵ espaciados, y el silencio \|S\|
  coincide con el de los ceros reales (0.010–0.014). Sin C₀: 2–12 ceros perdidos y \|S\| = 0.7–1.4. Propiedad conocida
  de ζ, cuantificada con la métrica de E16.
- **E18** (`psi_escalada.py`, `barrido_K.py`): D(x) = (ψ(x) − x)/√x hasta 10⁷. RMS(D − fórmula explícita con 100 ceros)
  = 0.078 frente a 0.325 con 100 γ al azar; baja a 0.002 con 10⁶ ceros. Máximo de \|D\| con x ≥ 101: 0.757, el 94.7 %
  de la cota condicional a RH ln x (ln x − ln ln x)/(8π) (arXiv:2312.05628, leída, no rederivada). Cumplir una cota
  condicional a RH no prueba RH.
- **E19** (`neff_cue.py`): la correlación de espaciados a 10¹², 10²¹ y 10²² (−0.332 ± 0.008, −0.317 ± 0.009,
  −0.305 ± 0.011) sigue a CUE de tamaño efectivo N_eff = ln(E/2π)/√(12Λ), Λ = 1.57314 (Bogomolny–Bohigas–Leboeuf–
  Monastra 2006, arXiv:math/0602270), y tiende a −0.308 (GUE). A t bajo (N_eff ≈ 2–3) la fórmula no se aplica.
- **E20** (`r2_bk.py`): la corrección de Bogomolny–Keating a primer orden sobrecorrige R₂(s) a t ≲ 10⁶: la amplitud
  medida es 0.35–0.48 de la prevista y sube con la altura. A 10¹²–10²² queda bajo el ruido de una ventana. No refuta BK.
- **E21** (`bk_completa.py`): la fórmula completa de BK sin expandir (solo primos y ζ(1+iε), sin parámetros) reproduce
  R₂ a t bajo: χ²/bin 1.4 frente a 7.8 de GUE (t ≈ 4·10⁴) y 3.7 frente a 31.6 (t ≈ 6.5·10⁵). Resultado conocido en
  esencia; aquí reproducido. Queda un residuo a t ≈ 6.5·10⁵.
- **E22** (`pico_tau1.py`): K(τ = 1) ≈ 41 es la identidad exacta K(1) = \|Σ e^{2πi dₙ}\|²/N con dₙ el conteo de
  Turing; refleja que Var S(t) ≈ 0.11 es pequeña (Selberg), no un artefacto de discretización.
- **E23** (`connes_zst.py`, `connes_mp.py`): forma cuadrática de Weil QW_λ con solo los primos ≤ X en la base V_n; las
  raíces reales de Σ cₙ/(z − 2πn/L) (L = ln X) del autovector de menor autovalor reproducen los ceros de Riemann. Con
  120 dígitos: X = 13, N = 40 → 28 raíces, todas ceros, 0 espurias; X = 30, N = 60 → 33 de 33 ceros (30 a < 1e-30);
  X = 50, N = 190 → 138 de 140 (los 2 que faltan estaban junto a un polo y los descartaba el barrido; corregido en la
  v2: X = 50, N = 95 da 53 de 53). Con N = 120 reproduce las cifras del artículo de Connes–Consani–Moscovici
  (arXiv:2511.22755): X = 13, primer cero 2.44e-55 y cero 50 2.04e-3; X = 9, primer cero 1.6e-34 y cero 20 2.4e-7.
  Control (pesos barajados): falla. Es mi ruta (raíces de la función racional), no su operador de rango uno; no verifiqué
  el Teorema 1.1 ni la autoadjunción. Reformulación de RH, no prueba. A X = 100, N = 200 hace falta más de 300 dígitos.
- **E24** (`matroide_li.py`): relaciones lineales enteras entre las ordenadas γₙ (PSLQ): 0 en 113 subconjuntos de 8
  (80 dígitos, coeficientes ≤ 1e5) y 109 de 12 (120 dígitos, coeficientes ≤ 1e4). Controles: azar de precisión
  completa sin relación; relación plantada recuperada. Matroide libre dentro del alcance del método; no prueba LI.
- **E25** (`grafo_controles.py`): grafo cuántico completo K₆ con aristas de longitud ln p. ⟨r⟩ con 40 realizaciones:
  ln p 0.6004 ± 0.0017, longitudes al azar 0.5990 ± 0.0016, ln n 0.6000 ± 0.0017 (GUE 0.5996); con S simétrica 0.5310 ±
  0.0014 (GOE 0.5307). La clase la fija la simetría de la dispersión, no la aritmética. Recuento de niveles validado
  con la ley de Weyl (1.000).
- **E26** (`ps_bbl.py`, `var_alturas.py`): a t ≈ 4·10⁴ la P(s) de los ceros es más estrecha que GUE (χ²/bin 27.7, suelo
  del método 1.7; var(s) 0.1607 frente a 0.1801). Un CUE de N finito (N = 4–5) explica solo una parte (χ²/bin 19;
  var 0.168–0.172); la continuación a N_eff = 2 no normaliza. Hasta t ≈ 1.1·10⁶, el exceso de ⟨r⟩ sobre GUE es
  compatible con 1/ln²(t/2π) (⟨r⟩∞ = 0.5998); la asíntota de var(s) no queda determinada con ese rango. El fichero
  `ceros_1millon.txt` contiene 2 001 052 ceros (idéntico a `ceros_2M.txt.gz` descomprimido). A t ≈ 10¹²–10²² (una ventana de 10⁴ ceros por altura, `var_alto.py`) var(s) = 0.1767, 0.1795, 0.1758 (± 0.003 cada una; media ponderada 0.1772 ± 0.0016, GUE a −1.7σ) y ⟨r⟩ = 0.6000, 0.5972, 0.6019: compatible con GUE y con 1/ln² para ⟨r⟩; descarta la asíntota 0.1715 (ajuste 1/ln³, p = 0.003) pero no distingue las demás. Acotado, no cerrado.

## Retractaciones y errores corregidos

- Atractor extraño, "red cósmica" y Takens 3D: retirados como evidencia.
- Criterio de Li como prueba: retirado (circular).
- Constante de GOE usada por error en Σ² de GUE: corregida.
- Media acumulada en vez de local en el chequeo de M₁: corregida.
- Ruido de K(τ) en SYK subestimado (±0.05–0.1; real ±0.3–0.5).
- "Berry–Keating exige romper T porque los ceros son GUE": es una sugerencia (Schumayer–Hutchinson,
  arXiv:1101.3116), no un teorema.
- E18: un primer cociente 1.971 contra ln²x/8π era un artefacto (la cota solo vale desde x ≈ 73.2); una estimación de
  \|S\| ≈ 0.7 con M = N_RS resultó ser 1.40 al medirla.
- E23: el barrido de raíces de `connes_mp.py` descartaba los ceros a menos de ~0.01 de un polo 2πn/L; corregido
  (intervalos con polo partidos en dos). Una versión pegada de otra fuente omitía el término arquimediano W02: error de
  ensamblado, no un resultado.
- E24: el primer control al azar usaba reales de 53 bits, que sí tienen relación entera; corregido con precisión completa.
- E25: tres scripts previos de grafos (`grafo_primos.py`, `grafo_traza.py`, `grafo_traza2.py`) se descartaron: detección
  de niveles por umbral manual sin validar con Weyl, sin control de longitudes y con un test de picos que no discrimina.

## Abiertos

- Pendiente 0.555 frente a 0.500 en Var ln\|ζ\|.
- Desviaciones de K(τ) a τ pequeño (el pico en τ = 1 queda explicado por E22).
- Máximo de ln\|ζ\| en ventanas frente a ln ln T − ¾ ln ln ln T (Fyodorov–Hiary–Keating).
- Operador de Hilbert–Pólya: sin candidato aceptado. E23 reproduce los ceros desde los primos con una forma de Weil,
  pero es una reformulación de RH, no un operador demostrado.
- Predicción completa de P(s) a t bajo (E26) y asíntota de var(s); residuo de R₂ a t ≈ 6.5·10⁵ (E21).
- Si la razón κ de BK tiende a 1 con la altura (la señal cae bajo el ruido a gran altura).
- Análogo en cuerpos de funciones de las correcciones de BK.
- Sin leer a fondo: Connes (1998, trazas en el espacio de clases de adeles), Katz–Sarnak, billares de microondas con
  ferritas, la réplica arXiv:1705.06767 y la sección numérica completa de arXiv:2511.22755.

## Literatura leída

- Schumayer–Hutchinson, arXiv:1101.3116: los ceros son GUE a corto alcance; sugiere un sistema que rompa T.
- Sierra, arXiv:1601.01797: con H = xp la simetría temporal se rompe; no lo deriva de GUE.
- Bender–Brody–Müller, arXiv:1608.03679: operador no hermitiano; realidad de autovalores heurística.
- Bellissard, arXiv:1704.02644: objeta autoadjunción y pertenencia a L².
- Bogomolny–Bohigas–Leboeuf–Monastra, arXiv:math/0602270: N_eff y CUE (E19–E21).
- Cotas explícitas de ψ(x) − x, arXiv:2312.05628: leída, no rederivada (E18).
- Connes–Consani–Moscovici, arXiv:2511.22755 («Zeta Spectral Triples»): leídos el Teorema 1.1 y las cifras numéricas;
  fórmulas cerradas solo parcialmente (E23).
