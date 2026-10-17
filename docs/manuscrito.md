# El pulso fractal del infinito

### Una lectura geométrica y espectral de la Hipótesis de Riemann

*Benjamín Cabeza Durán — Versión 0.3*

---

## I. El enigma de los bloques fundamentales

Los números primos son el esqueleto indivisible de la aritmética. Pese a su
simplicidad —mayores que uno, divisibles solo por sí mismos y por la
unidad—, su distribución ha desafiado la intuición humana durante siglos:
parecían chispas repartidas al azar sobre un lienzo infinito. En 1859,
Bernhard Riemann intuyó que, tras esa aparente anarquía, subyacía una
simetría oculta cuya clave está en el plano complejo.

## II. La línea crítica

La función zeta, $\zeta(s)$, traduce la distribución de los primos a una
geometría en el plano complejo. Sobre la línea crítica, $\Re(s)=\tfrac12$,
la trayectoria $\Gamma(t)=\zeta(\tfrac12+it)$ no es una sucesión estática de
puntos sino un movimiento continuo que gira en torno al origen; cada paso de
esta curva por el punto $(0,0)$ marca un cero no trivial.

La Hipótesis de Riemann afirma que todos esos ceros ocurren exactamente
sobre esa línea. La pregunta que este programa explora no es si eso es así
—eso ya lo conjetura la hipótesis— sino si la propia geometría de
$\Gamma(t)$ contiene una estructura (¿un atractor?, ¿una forma de
autosimilitud?) capaz de forzarlo. Esa estructura, si existe, todavía no se
ha medido ni demostrado: es la hipótesis de trabajo, no un hecho ya
observado.

## III. Del caos numérico al rigor determinista

Verificar computacionalmente miles de millones de ceros no basta para
demostrar la hipótesis: un contraejemplo podría aparecer más allá de lo
calculado. El programa de Hilbert–Pólya propone otra vía: si existiera un
operador autoadjunto cuyo espectro reprodujera las alturas $t_n$ de los
ceros, la prueba dejaría de depender de revisar caso por caso y pasaría a
depender de una propiedad estructural del espacio en que vive la función.

Es una vía real, pero sigue siendo un programa: hoy no se conoce tal
operador, y no puede afirmarse todavía que una geometría fractal «prohíba»
ceros fuera de la línea. Lo que sí puede formularse con precisión es la
pregunta.

## IV. Conclusión de la lectura geométrica

La Hipótesis de Riemann no es solo un enunciado sobre primos: es una
pregunta sobre la estructura del continuo. Leer la función zeta como un
pulso con fase, amplitud y —quizá— estructura de escala tiende un puente
entre la teoría de números, la física y la geometría de sistemas dinámicos.
Si algún día se demuestra que esa estructura obliga a los ceros a
permanecer sobre la línea crítica, se habrá encontrado la regla que
organiza el aparente caos de los primos. Hasta entonces, es una hipótesis
de trabajo — no un hecho establecido.

## V. Resultados: estadística espectral (N=5000)

Separado de la lectura geométrica anterior —que sigue siendo exploratoria—,
el pipeline (`pfi01.py` + `pfi01_espectral.py`) compara los espaciados entre
ceros contra predicciones exactas de teoría de matrices aleatorias (RMT).
Con los primeros 5000 ceros no triviales:

- **P(s<0.3) = 1.48%**, frente al 2.73% que predice GUE y ~26% que
  predeciría un proceso sin correlación (Poisson) — desviación ≈5.4σ.
- **⟨r⟩ = 0.616** (razón de espaciados consecutivos, sin normalizar por
  densidad), frente al 0.5996 asintótico de GUE — desviación ≈4.2σ tras
  corregir por la autocorrelación entre pares consecutivos.
- **Δ₃(L)**, el estadístico de rigidez de largo alcance, confirma la misma
  dirección: el espectro es más rígido que GUE en L∈[5,50], y la curva
  empírica se aplana mientras la predicción GUE sigue creciendo
  logarítmicamente.

Los tres estadísticos —dos de corto alcance, uno de largo alcance— apuntan
en la misma dirección: más rigidez que el límite asintótico de GUE. No lo
interpretamos como una anomalía nueva: es consistente, en dirección y orden
de magnitud, con la corrección de altura finita a la convergencia GUE
descrita por Bogomolny–Keating (1996) y Berry–Keating (1999) — a T finito
el espectro de ceros de zeta es más rígido que en el límite T→∞, y esa
corrección decae con T. Al dividir la muestra en t<1500 y t≥1500, ambas
ventanas se mueven en la dirección esperada (acercándose a GUE al aumentar
t).

Esto es, hasta ahora, el único resultado del proyecto que no es solo
instrumentación: tiene comparación teórica cerrada y significancia
cuantificable — modesta frente a los cientos de millones de ceros de la
literatura estándar (Odlyzko), pero real.

---

## Apéndice matemático

### 1. Símbolos y fundamentos

$s=\sigma+it$, con $\sigma=\Re(s)$ y $t=\Im(s)$. Para $\Re(s)>1$,

$$\zeta(s)=\sum_{n=1}^{\infty}\frac{1}{n^s},$$

extendida a todo el plano por continuación analítica (salvo un polo simple
en $s=1$). La línea crítica es $\Re(s)=\tfrac12$; la hipótesis afirma que
todo cero no trivial cumple $\boxed{\sigma=\tfrac12}$.

### 2. Geometría: rotación y fase

Sobre la línea crítica, $\Gamma(t)=\zeta(\tfrac12+it)=R(t)e^{i\phi(t)}$, con
$R(t)=|\zeta(\tfrac12+it)|$ y $\phi(t)=\arg\zeta(\tfrac12+it)$. Un cero
ocurre exactamente cuando $R(t_n)=0$.

### 3. Las funciones Ξ(t) y Z(t)

La función completada
$\xi(s)=\tfrac12 s(s-1)\pi^{-s/2}\Gamma(\tfrac{s}{2})\zeta(s)$ satisface
$\xi(s)=\xi(1-s)$. Definiendo $\Xi(t):=\xi(\tfrac12+it)$, se tiene
$\Xi(t)\in\mathbb{R}$ para $t\in\mathbb{R}$ y $\Xi(t)=\Xi(-t)$. La función
de Riemann–Siegel,

$$Z(t)=e^{i\theta(t)}\zeta(\tfrac12+it),$$

también es real para $t$ real, con $|Z(t)|=|\zeta(\tfrac12+it)|$ y
$\zeta(\tfrac12+it)=Z(t)e^{-i\theta(t)}$: una amplitud real multiplicada
por una fase.

### 4. Hilbert–Pólya

Si existiera un operador autoadjunto $\hat H$ cuyo espectro $\{t_n\}$
coincidiera con las alturas de los ceros, la hermiticidad garantizaría que
esos $t_n$ son reales. Como cada cero se escribe $\rho_n=\tfrac12+it_n$,
eso es justo lo que fijaría su parte real en $1/2$:
$\boxed{\Re(\rho_n)=\tfrac12}$. No es que la parte imaginaria se anule —
nunca lo hace: $t_n$ es la altura del cero, y crece sin límite (el primero
está en $t\approx14.13$). Tal operador no se conoce; encontrarlo, o
descartar que exista, sigue siendo un problema abierto.

### 5. Estadística espectral (RMT)

Para dos ceros consecutivos con alturas $t_n<t_{n+1}$, el espaciado
normalizado usa la densidad local $\rho(t)\approx\frac{1}{2\pi}\log(t/2\pi)$:

$$s_n = (t_{n+1}-t_n)\cdot\rho(t_n).$$

La razón de espaciados consecutivos, que no requiere normalizar,

$$r_n = \frac{\min(\Delta_n,\Delta_{n+1})}{\max(\Delta_n,\Delta_{n+1})},
\qquad \Delta_n = t_{n+1}-t_n,$$

tiene valor asintótico conocido: $\langle r\rangle_{\text{GUE}}=0.5996$,
$\langle r\rangle_{\text{Poisson}}=2\ln2-1\approx0.386$ (Atas et al. 2013).

La rigidez espectral de Dyson–Mehta, verificada por triangulación en tres
fuentes independientes,

$$\Delta_3^{\text{GUE}}(L) \approx \frac{\ln(2\pi L)+\gamma-\tfrac54}{2\pi^2}
\approx \frac{\ln L + 1.1651}{2\pi^2}, \qquad
\Delta_3^{\text{Poisson}}(L) = \frac{L}{15}.$$

---

## Notas de revisión

**v0.1 → v0.2.** Dos correcciones a la versión original del ensayo:

1. En Hilbert–Pólya, la hermiticidad fija la parte real del cero en 1/2 —
   no anula su parte imaginaria, que nunca es cero.
2. Atractor, autosimilitud y la «prohibición» topológica de ceros fuera de
   la línea se presentan como hipótesis por comprobar, no como hechos ya
   observados.

**v0.2 → v0.3.** Se añade la sección V con los resultados espectrales de
N=5000. En el camino se corrigieron cuatro errores numéricos — dos de
revisión externa sin verificar (una desviación estándar de ⟨r⟩ asumida en
vez de medida, que inflaba la significancia de ~4σ a ~14σ; una fórmula de
Δ₃ con el prefactor de GOE aplicado a GUE) y dos propios en la primera
implementación de Δ₃ (ajuste sobre los puntos muestreados en vez de sobre
la integral continua; una división por L de más). Ninguno se detectó por
inspección: los cuatro salieron de recalcular con los datos en mano. Un
preprint reciente (Alarcón 2026, sin revisión por pares) propone un ajuste
numérico específico para la convergencia de ⟨r⟩ con T; se menciona como
referencia tentativa, no como fuente del mecanismo, que corresponde a
Bogomolny–Keating (1996) y Berry–Keating (1999).
