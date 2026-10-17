# El pulso fractal del infinito
## Una lectura geométrica y espectral de la Hipótesis de Riemann
Autor: Benjamín Cabeza Durán  
Versión: 0.1 — Programa de investigación

---

## 1. Planteamiento

[texto de §1 tal cual lo escribiste, copiado acá]

## 2. El espacio complejo como espacio dinámico

[§2]

## 3. Radio y fase

[§3]

## 4. El evento singular: el cero

[§4]

## 5. La transformación de Riemann–Siegel

[§5]

## 6. De la metáfora a una hipótesis falsable

[§6, PFI-1]

## 7. ¿Es realmente fractal?

[§7]

> **Nota metodológica (v0.1)**: un valor `1 < D_B < 2` no implica fractalidad.
> La pendiente `D_B ≈ 1.803` medida con N=1000 es **descriptiva** y depende
> de la densidad de muestreo. Se reporta como observación, no como
> propiedad matemática.

## 8. El segundo nivel: autosimilitud

[§8]

## 9. El tercer nivel: los ceros como eventos espectrales

[§9]

## 10. El puente hacia Hilbert–Pólya

[§10, con la corrección: `Re(ρ_n)=1/2`, NO `Im(ρ_n)=0`]

## 11. Hipótesis PFI-2

[§11]

## 12. El experimento computacional

Tabla de observables por instante t:

| t | Z(t) | \|ζ(1/2+it)\| | φ(t) | κ(t) | v(t) |
|---:|---:|---:|---:|---:|---:|
| t_n | 0 | 0 | — | — | — |
| t_n + ε | · | · | · | · | · |

Implementación: `pfi01.py`, `lab_phase.py`.

## 13. El criterio decisivo

> No buscaremos solamente patrones. Buscaremos contraejemplos.

Controles implementados: `lab_controls.py` (GUE, Poisson, gaps barajados).
Control pendiente (v0.2): función L de Dirichlet con carácter no trivial.

## 14. Arquitectura conceptual

Primos → ζ(s) → Ξ(t) → Z(t) → Γ(t) → ceros → estadística → Ĥ  
Γ(t) → fase → amplitud → escala → geometría → invariantes

Objetivo: GEOMETRÍA ∩ ESPECTRO.

---

## Resultados numéricos v0.1 (a completar)

Ver `docs/resultado_v01.md`. Se completa **después** de correr los tests
de `docs/criterios_falsacion.md`, con los números tal como salen.
