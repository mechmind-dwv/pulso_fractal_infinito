#!/usr/bin/env python3
"""Pasa README.md de v0.4 a v0.5 (inserta secciones; no borra nada). Uso: python3 actualiza_readme_v05.py [README.md]"""
import sys
f = sys.argv[1] if len(sys.argv) > 1 else "README.md"
t = open(f, encoding="utf-8").read()
if "v0.5" in t:
    sys.exit("Ya es v0.5; no hago nada.")

t = t.replace("**v0.4**", "**v0.5**", 1)

NUEVO = '''## Novedades de la v0.5: del «qué» al «cómo» (ensayos E10–E18)

Los ceros ya están medidos; aquí se pregunta de dónde sale esa estructura. Cada ensayo lleva dato, teoría exacta y un control. Detalle en `docs/ENSAYOS.md`, `ensayo3_enteros_y_ceros/README.md` y `analisis_avanzado/README.md`.

| Ensayo | Resultado | Tipo |
|---|---|---|
| SYK N=22 frente a los ceros | Misma clase GUE en ⟨r⟩, P(s), Δ₃ y R₂ (RMS 0.020; control N=16 GOE: 0.089). Comparten clase estadística, **no origen** | Positivo |
| `xp` en una caja (Berry–Keating) | Da un peine rígido (⟨r⟩ = 1), no GUE. Parte del artefacto de rejilla está documentado | Negativo |
| Jacobi inverso con 80 ceros | Indistinguible de GUE (rms 3.02 frente a 3.70 ± 0.38). No identifica un operador | Negativo |
| Positividad de Weil | La matriz construida solo con primos coincide con la de los ceros a 1e-10; autovalor mínimo > 0 | Positivo (comprobación, no prueba) |
| ψ(x) desde los ceros | Error RMS de reconstrucción: 0.743 con 10 ceros, 0.022 con 10⁵; control con γ al azar ≈ 2.5 | Positivo |
| Enteros → ceros (`ensayo3_enteros_y_ceros/`) | Con M = N_RS = ⌊√(t/2π)⌋ enteros se fijan el 96.8 % de los ceros a t = 10⁴ | Positivo |
| GUE + primos insertados | Falla: ⟨r⟩ baja a 0.580 y el «silencio» de la suma sobre ceros es 2.68 frente a 0.01 en los ceros reales | Negativo |
| Riemann–Siegel con C₀ | Reproduce las posiciones de los ceros a ≲ 1.5×10⁻⁴ espaciados y el silencio (0.010–0.014, igual que los reales). Propiedad conocida, aquí cuantificada | Positivo |
| E18: x^(−1/2)(ψ(x) − x) hasta 10⁷ | Ver abajo | Positivo |

### E18. La oscilación de los primos con 100 ceros

`psi_escalada.py` (solo numpy; matplotlib únicamente para dibujar). D(x) = (ψ(x) − x)/√x, con ψ exacta por criba hasta 10⁷.

| | valor |
|---|---|
| rango de D | [−0.770, 0.673] |
| RMS(D − fórmula explícita con 100 ceros) | 0.078 |
| RMS(D − control con 100 γ al azar) | 0.325 |
| RMS(D) | 0.243 |
| máx. \\|D\\| con x ≥ 101 | 0.757 |

La cota condicional a RH, ln x (ln x − ln ln x)/(8π) para x ≥ 101 (Schoenfeld y Grenié–Molteni; leída en arXiv:2312.05628, no rederivada), vale 8.55 en x = 10⁷; el máximo medido es el 94.7 % de la cota, en x ≈ 101, y queda muy por debajo después. Cumplir una cota **condicional** a RH no prueba RH; hasta 10⁷ no podía ser de otro modo.

Una nota de método: una primera ejecución dio un cociente 1.971 contra ln²x/8π. Era un artefacto: esa cota solo vale desde x ≈ 73.2 y la rejilla empezaba en x = 20.

**Corrección de afirmaciones anteriores:** el ruido de K(τ) en SYK es ±0.3–0.5 (antes escribí ±0.05–0.1); una estimación mía de |S| ≈ 0.7 para Riemann–Siegel con M = N_RS resultó ser 1.40 al medirla.

---

## Contexto externo (25 nov 2026)

El 20 de noviembre de 2026 se anunció una prueba formalizada en Lean de la Hipótesis de Riemann (en realidad de la generalizada), atribuida a Levent Alpöge y descrita como hallada con ayuda de Claude. Según lo publicado, compila, usa tres axiomas estándar y la han aceptado varios verificadores; quedan pendientes una ejecución «paranoica» del kernel y el Comparator. Está en revisión por la comunidad. Este repositorio **no la da por cerrada ni por refutada**, y no la usa ni depende de ella: no he podido consultar la declaración original del autor, y los recuentos de ficheros y líneas que circulan no coinciden entre fuentes. Si se confirma, responde al «qué» y a RH, pero no hace innecesarias las preguntas de este programa (de dónde sale la estructura GUE de los ceros), que son de otro tipo.

---

'''
marca = "## Resumen de la v0.4"
assert marca in t
t = t.replace(marca, NUEVO + marca, 1)
t = t.replace("## Resumen de la v0.4", "## Resumen de la v0.4 (base, sin cambios)", 1)

t = t.replace("## Reproducibilidad\n", "## Reproducibilidad\n\nLos ensayos nuevos se ejecutan así (los scripts de ψ(x) están en `analisis_avanzado/`; los de enteros y ceros, con `run_all.sh`, en `ensayo3_enteros_y_ceros/`):\n\n```bash\ncd analisis_avanzado\npython3 psi_escalada.py 1e7 ceros_100k.txt 100   # E18 -> psi_escalada.png\ncd ../ensayo3_enteros_y_ceros && bash run_all.sh\n```\n\nLos de v0.4:\n", 1)
open(f, "w", encoding="utf-8").write(t)
print("README actualizado a v0.5")
