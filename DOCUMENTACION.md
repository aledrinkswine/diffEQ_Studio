# DiffEq Studio — Documentación

> Aplicación de escritorio para Windows que resuelve **ecuaciones diferenciales ordinarias paso a paso**, las verifica y las grafica. Está pensada para el curso de **Cálculo 4 (Ecuaciones Diferenciales)**.

---

## Índice

1. [¿Para qué sirve?](#1-para-qué-sirve)
2. [Características principales](#2-características-principales)
3. [Guía de uso](#3-guía-de-uso)
4. [Fundamento matemático de cada módulo](#4-fundamento-matemático-de-cada-módulo)
5. [Cómo se hizo](#5-cómo-se-hizo)
6. [Arquitectura del código](#6-arquitectura-del-código)
7. [Decisiones de diseño importantes](#7-decisiones-de-diseño-importantes)
8. [Pruebas](#8-pruebas)
9. [Compilar el ejecutable](#9-compilar-el-ejecutable)
10. [Limitaciones conocidas](#10-limitaciones-conocidas)
11. [Trabajo futuro](#11-trabajo-futuro)

---

## 1. ¿Para qué sirve?

Las calculadoras simbólicas normales dan **la respuesta**, pero no **el procedimiento**. DiffEq Studio está hecho para estudiar y para verificar tareas. Para cada ecuación:

1. **Identifica el tipo** de ecuación (separable, lineal, Bernoulli, homogénea, exacta, coeficientes constantes…).
2. **Explica el método** y muestra cada paso intermedio, como se haría en el cuaderno: factor integrante, ecuación característica, Wronskiano, etc.
3. **Verifica** la solución sustituyéndola en la ecuación original.
4. **Resuelve el problema de valor inicial** si se dan condiciones iniciales.
5. **Grafica** la solución, la familia de soluciones o el campo de direcciones.

Además incluye la **transformada de Laplace** (directa, inversa, PVI y convolución), **sistemas de EDO** con retrato de fase, un laboratorio de **métodos numéricos** (Euler, Heun y RK4) y un **campo de direcciones interactivo**.

Usos típicos:
- Comprobar ejercicios resueltos a mano.
- Ver el procedimiento de un método que todavía no se domina.
- Visualizar el comportamiento de las soluciones (equilibrios, asíntotas, oscilaciones).
- Comparar la precisión de los métodos numéricos.

---

## 2. Características principales

| Módulo | Contenido |
|---|---|
| **Primer orden** | Integración directa · Variables separables (con soluciones de equilibrio) · Lineal con factor integrante · Bernoulli · Homogénea (y = v·x) · Exacta · Exacta con factor integrante μ(x) o μ(y) |
| **Orden superior** | Lineales con coeficientes constantes de **cualquier orden**: ecuación característica, discriminante, raíces reales, repetidas y complejas, conjunto fundamental, Wronskiano y **variación de parámetros** para la parte no homogénea |
| **Laplace** | Transformada directa (tabla, linealidad, 1.er y 2.º teorema de traslación, derivada de la transformada, impulso) · inversa con fracciones parciales y completando el cuadrado · PVI con escalón u(t−a) e impulso δ(t−a) · convolución y su teorema |
| **Sistemas** | Lineales X' = AX (+ F(t)) de cualquier tamaño: valores y vectores propios reales, complejos y repetidos (vector generalizado), matriz exponencial, variación de parámetros, clasificación del punto crítico · No lineales: equilibrios, jacobiano y linealización · Retrato de fase interactivo con nulclinas y vectores propios |
| **Métodos numéricos** | Euler · Euler mejorado (Heun) · Runge-Kutta 4 · tabla por iteración · error absoluto frente a la solución exacta · gráfica comparativa |
| **Campo de direcciones** | Campo de pendientes de y' = F(x, y) · clic para trazar la curva solución por cualquier punto · isoclina y' = 0 (equilibrios) |
| **Ayuda** | Sintaxis de entrada y resumen de los módulos |

Además:
- **Verificación automática** de cada solución.
- **52 ejemplos precargados** que se resuelven con un clic.
- **Cálculo en segundo plano**: la ventana no se congela mientras resuelve.
- Gráficas con zoom, desplazamiento y opción de guardar como imagen (barra de Matplotlib).
- Un único `.exe` que no requiere instalar Python.

---

## 3. Guía de uso

### 3.1 Abrir la aplicación
Haz doble clic en `dist/DiffEqStudio.exe`. La primera apertura tarda unos segundos porque el ejecutable descomprime sus librerías internas.

### 3.2 Sintaxis de entrada

| Qué | Cómo se escribe | Ejemplo |
|---|---|---|
| Derivadas | `y'`, `y''`, `y'''`, `dy/dx` | `y'' + 4y = 0` |
| Multiplicación | `*` opcional entre número y letra | `3x*y`, `2y'` |
| Potencias | `^` | `y^3`, `x^(1/2)` |
| Exponencial | `e^(...)` o `exp(...)` | `e^(-2x)` |
| Funciones | `sin cos tan sen tg ln log sqrt arctan` | `tan(x)` |
| Forma diferencial | `M dx + N dy = 0` | `(2x*y + 3) dx + (x^2 - 1) dy = 0` |
| Condiciones iniciales | `y(x0)=valor, y'(x0)=valor` | `y(0)=1, y'(0)=0` |

> Para un producto de dos variables hay que usar `*`: escribe `x*y`, no `xy`.

### 3.3 Pestañas *Primer orden* y *Orden superior*
1. Escribe la ecuación (o elige una de la lista **Ejemplos**).
2. Opcionalmente escribe las condiciones iniciales.
3. Presiona **Resolver ▶** o Enter.
4. A la izquierda aparece el procedimiento y a la derecha la gráfica.

Qué se grafica:
- **Primer orden:** el campo de direcciones con la solución particular, o una familia de curvas si no hay condición inicial.
- **Orden superior:** la solución particular, o varias soluciones de la familia (C₁, C₂ con distintos valores).

### 3.4 Pestaña *Laplace*
Arriba eliges uno de los cuatro modos:

| Modo | Qué escribir | Ejemplo |
|---|---|---|
| Transformada ℒ{f(t)} | f(t) | `e^(2t) sin(3t)`, `u(t-2) e^(t-2)`, `delta(t-1)` |
| Inversa ℒ⁻¹{F(s)} | F(s) | `(2s + 1)/(s^2 + 5s + 6)`, `e^(-2s)/(s^2 + 1)` |
| PVI con Laplace | ecuación + condiciones en t = 0 | `y'' + 4y = delta(t - pi)` con `y(0)=0, y'(0)=1` |
| Convolución | f(t) y g(t) | `e^t` y `sin(t)` |

El escalón unitario se escribe `u(t-a)` (también `H(t-a)`) y el impulso `delta(t-a)`. La gráfica muestra f(t), y(t) o las tres funciones de la convolución. Los impulsos se dibujan como flechas.

### 3.5 Pestaña *Sistemas*
Escribe una ecuación por línea, por ejemplo:
```
x' = x + 2y
y' = 3x + 2y
```
Las condiciones iniciales se escriben como `x(0)=1, y(0)=0`. También se acepta `dx/dt = …`, y **Ctrl+Enter** resuelve desde el cuadro del sistema.

Para sistemas 2×2 autónomos se dibuja el **retrato de fase** con:
- el flujo del campo,
- las nulclinas x' = 0 e y' = 0 (líneas punteadas),
- las rectas de los vectores propios,
- los equilibrios (puntos negros),
- la trayectoria con X(0).

**Haz clic** en el retrato para trazar la trayectoria que pasa por ese punto. Si hay condiciones iniciales, abajo se grafican x(t), y(t).

### 3.6 Pestaña *Métodos numéricos*
Escribe F(x, y) en `y' =`, el punto inicial (x₀, y₀), el paso `h` y el `x final`, y presiona **Calcular**. La tabla muestra cada iteración de los tres métodos. Si existe solución exacta, también muestra su valor y el error de Euler y de RK4.

### 3.7 Pestaña *Campo de direcciones*
Escribe F(x, y) y la ventana (x min/max, y min/max) y presiona **Dibujar**. **Haz clic** en cualquier punto de la gráfica para trazar la solución que pasa por él, hacia adelante y hacia atrás. La línea punteada naranja es la isoclina y' = 0. **Borrar curvas** redibuja el campo limpio.

---

## 4. Fundamento matemático de cada módulo

### 4.1 Primer orden

La ecuación se lleva primero a la **forma normal** y' = F(x, y) despejando y'. Después se prueban los métodos **en este orden** y se usa el primero que aplica:

| # | Método | Criterio de detección | Procedimiento mostrado |
|---|---|---|---|
| 1 | Integración directa | F no depende de y | y = ∫F dx + C |
| 2 | Variables separables | F = g(x)·h(y) (`separatevars`) | ∫dy/h(y) = ∫g(x)dx, despeje de y, soluciones de equilibrio h(y)=0 |
| 3 | Lineal | ∂F/∂y no depende de y | y' + P(x)y = Q(x), μ = e^{∫P dx}, (μy)' = μQ, integración y despeje |
| 4 | Bernoulli | F = −P(x)y + Q(x)yⁿ, n ≠ 0, 1 | v = y^{1−n} convierte la ecuación en lineal en v; se resuelve y se regresa a y |
| 5 | Homogénea | F(tx, ty) = F(x, y) | y = vx ⇒ v + xv' = F(1, v), separación de variables, regreso a v = y/x |

**Ecuaciones exactas** (cuando se escriben como `M dx + N dy = 0`):
1. Se prueba el criterio ∂M/∂y = ∂N/∂x.
2. Si no se cumple, se busca un factor integrante:
   - μ(x) = e^{∫(M_y − N_x)/N dx} si ese cociente depende solo de x.
   - μ(y) = e^{∫(N_x − M_y)/M dy} si ese otro cociente depende solo de y.
3. Se construye la función potencial ψ = ∫M dx + g(y), con g'(y) = N − ∂ψ/∂y.
4. La solución general es ψ(x, y) = C, y se intenta despejar y.

Si ninguna regla aplica, la ecuación se resuelve con `dsolve` de SymPy y se indica el método que SymPy detectó (Riccati, etc.).

### 4.2 Orden superior: coeficientes constantes

Para aₙy⁽ⁿ⁾ + … + a₁y' + a₀y = g(x):

**1) Solución homogénea**
- Se propone y = e^{rx} ⇒ ecuación característica aₙrⁿ + … + a₀ = 0.
- Para orden 2 se muestra el discriminante Δ = b² − 4ac y se clasifica el caso.
- Cada raíz aporta funciones al conjunto fundamental:

| Raíz | Funciones |
|---|---|
| Real r, multiplicidad m | e^{rx}, x·e^{rx}, …, x^{m−1}e^{rx} |
| Compleja α ± βi, multiplicidad m | x^j e^{αx}cos βx, x^j e^{αx}sin βx, para j = 0…m−1 |

- y_h = C₁y₁ + … + Cₙyₙ

**2) Solución particular — variación de parámetros**

Se usa este método porque funciona para **cualquier** g(x), por ejemplo tan x o sec x, y no solo para las formas de coeficientes indeterminados:
- Matriz y determinante wronskiano W.
- Se resuelve el sistema W·u' = (0, …, 0, g/aₙ) y se integra cada uᵢ'.
- y_p = Σ uᵢyᵢ. Los términos que ya pertenecen a y_h se eliminan.

**3) Solución general** y = y_h + y_p

Las ecuaciones de orden superior que **no** tienen coeficientes constantes (Cauchy-Euler, series, etc.) se resuelven con SymPy y se muestra el nombre del método.

### 4.3 Problema de valor inicial
Se forma el sistema yᵏ(x₀) = valor para cada condición y se despejan las constantes C₁…Cₙ. Si la solución es implícita (ψ = C), la constante se calcula sustituyendo el punto inicial. Si la solución tiene varias ramas (por ejemplo ±√…), se elige la que cumple la condición.

### 4.4 Transformada de Laplace

**Directa.** Por linealidad, cada término c·g(t) se transforma por separado y se indica qué regla de la tabla se aplica:

| Forma de g(t) | Regla |
|---|---|
| 1, tⁿ, e^(at), sin bt, cos bt, sinh bt, cosh bt | Tabla básica |
| e^(at)·f(t) | 1.er teorema de traslación: F(s − a) |
| tⁿ·f(t) | Derivada de la transformada: (−1)ⁿ F⁽ⁿ⁾(s) |
| u(t − a)·g(t) | 2.º teorema de traslación: e^(−as)·ℒ{g(t + a)} |
| δ(t − a) | e^(−as) |

**Inversa.**
1. Se separan los términos con e^(−as) y se agrupan por su desplazamiento a.
2. Cada grupo se descompone en **fracciones parciales** (`apart`).
3. A cada fracción se le aplica la regla que corresponde:
   - Factor lineal (s − a)ⁿ: tⁿ⁻¹e^(at)/(n−1)!.
   - Cuadrático irreducible: se **completa el cuadrado** y se usa e^(αt)cos βt o e^(αt)sin βt.
4. Los grupos con e^(−as) se desplazan: g(t − a)·u(t − a).

**PVI.**
1. Se transforman ambos lados con ℒ{y⁽ᵏ⁾} = sᵏY − sᵏ⁻¹y(0) − … − y⁽ᵏ⁻¹⁾(0).
2. Se despeja Y(s) de la ecuación algebraica.
3. Se aplica la inversa paso a paso.

La solución se verifica:
- con `checkodesol` si no hay discontinuidades;
- numéricamente en varios puntos lejos de los saltos si hay u(t−a) o δ(t−a), porque la verificación simbólica no es fiable con esas funciones.

Las condiciones iniciales se comprueban siempre.

**Convolución.** Se calcula (f ∗ g)(t) = ∫₀ᵗ f(τ)g(t−τ)dτ y se comprueba el teorema ℒ{f ∗ g} = F(s)·G(s).

### 4.5 Sistemas de EDO

**Lineales X' = AX + F(t).** La app detecta si el sistema es lineal con coeficientes constantes (∂fᵢ/∂xⱼ constantes) y obtiene A y F(t).

1. Ecuación característica det(A − λI) = 0. En 2×2 se muestran también la traza, el determinante y tr² − 4·det.
2. Soluciones fundamentales según el tipo de valor propio:

| Valor propio | Soluciones |
|---|---|
| Real λ con vector v | v·e^(λt) |
| Complejo α + βi con v = a + bi | e^(αt)(a cos βt − b sin βt), e^(αt)(a sin βt + b cos βt) |
| Repetido con un solo vector | Vector generalizado (A − λI)w = v ⇒ e^(λt)(t·v + w) |
| Multiplicidad mayor deficiente | Matriz exponencial e^(At) |

3. Si hay F(t), se aplica **variación de parámetros**: Φ(t) = [X₁ … Xₙ], u' = Φ⁻¹F y X_p = Φ·∫u' dt.
4. Verificación X' − AX − F = 0 y cálculo de las constantes con X(t₀).

**Clasificación del punto crítico (2×2):**

| Valores propios | Tipo | Estabilidad |
|---|---|---|
| Reales, mismo signo, distintos | Nodo | estable si λ < 0 |
| Reales, signos opuestos | Punto silla | inestable |
| Repetidos | Nodo propio (estrella) o impropio | estable si λ < 0 |
| Complejos, α ≠ 0 | Foco (espiral) | estable si α < 0 |
| Imaginarios puros | Centro | estable, no asintóticamente |
| Algún λ = 0 | Degenerado | recta de equilibrios |

**No lineales.**
1. Se buscan los equilibrios resolviendo f(X) = 0.
2. Se calcula el **jacobiano** J.
3. Cada equilibrio se clasifica con los valores propios de J en ese punto (linealización, teorema de Hartman-Grobman).
4. La app advierte que un centro o λ = 0 no son concluyentes para el sistema no lineal.

Las trayectorias se calculan numéricamente con RK4.

### 4.6 Métodos numéricos

Con paso h y xₙ₊₁ = xₙ + h:

| Método | Fórmula | Orden del error global |
|---|---|---|
| Euler | yₙ₊₁ = yₙ + h·f(xₙ, yₙ) | O(h) |
| Heun (Euler mejorado) | k₁ = f(xₙ, yₙ), k₂ = f(xₙ + h, yₙ + hk₁), yₙ₊₁ = yₙ + h(k₁ + k₂)/2 | O(h²) |
| RK4 | yₙ₊₁ = yₙ + h(k₁ + 2k₂ + 2k₃ + k₄)/6 | O(h⁴) |

La solución exacta, cuando existe, se calcula con `dsolve` para medir el error.

### 4.7 Campo de direcciones
En una malla de puntos (x, y) se dibuja un segmento de pendiente F(x, y), normalizado para que todos tengan el mismo largo. La curva al hacer clic se integra con RK4 en ambas direcciones. La integración se detiene si la solución diverge.

---

## 5. Cómo se hizo

### 5.1 Tecnologías

| Herramienta | Versión | Para qué |
|---|---|---|
| Python | 3.14 | Lenguaje |
| **SymPy** | 1.14 | Álgebra simbólica: parseo, derivadas, integrales, `solve`, `dsolve`, `checkodesol`, *pretty printing* |
| **NumPy** | 2.5 | Cálculo numérico vectorizado |
| **Matplotlib** | 3.11 | Gráficas embebidas en la ventana |
| **Tkinter / ttk** | (incluido en Python) | Interfaz gráfica |
| **PyInstaller** | 6.22 | Empaquetado en un único `.exe` |

Se eligió **Tkinter** porque viene con Python (no añade dependencias) y produce ejecutables más pequeños que Qt. Se eligió **SymPy** porque es el único motor simbólico libre y completo en Python.

### 5.2 Proceso de desarrollo

1. **Definición del alcance (versión mínima):** primer orden, orden superior con coeficientes constantes, métodos numéricos y campo de direcciones.
2. **Motor matemático (`solvers.py`):** la parte central. SymPy da la respuesta, pero **no los pasos**, así que cada método se programó a mano replicando el procedimiento del cuaderno. Cada paso se guarda como texto o como expresión simbólica.
3. **Pruebas del motor** con los 19 ejemplos (`test_solvers.py`) antes de hacer la interfaz. Se corrigieron decimales en condiciones iniciales (convertidos a fracciones exactas) y nombres de métodos.
4. **Métodos numéricos (`numerics.py`):** Euler, Heun, RK4, trazado de curvas y campo de pendientes.
5. **Interfaz (`app.py`):** pestañas, área de procedimiento con formato, gráficas embebidas y cálculo en un hilo secundario.
6. **Prueba de humo de la interfaz (`test_gui.py`):** abre la app, resuelve todos los ejemplos y guarda las gráficas como imágenes para revisarlas. Así se detectaron y corrigieron dos fallos visuales:
   - Las asíntotas verticales se dibujaban como líneas.
   - La familia de `y'' + y = tan(x)` no se graficaba porque aparecía un logaritmo de un número negativo.
7. **Empaquetado** con PyInstaller y verificación de que el `.exe` abre su ventana.

**Segunda parte (Laplace y sistemas):**

8. **`laplace.py`:** transformada directa con la regla de la tabla para cada término, inversa con fracciones parciales, PVI y convolución. Se corrigieron dos problemas:
   - `expand` movía e^(−as) al denominador como 1/e^(as). Ahora cada e^(cs) se reemplaza por un símbolo antes de expandir.
   - `limit` y `checkodesol` se quedaban colgados con funciones escalón. Ahora las condiciones se evalúan por sustitución y la ecuación se verifica numéricamente lejos de los saltos.
9. **`systems.py`:** lectura de sistemas escritos ecuación por ecuación, valores y vectores propios con todos los casos, variación de parámetros matricial, clasificación del punto crítico y linealización de sistemas no lineales. Las ecuaciones con matrices se muestran con `Eq(..., evaluate=False)`, porque SymPy evaluaba «símbolo = matriz» como `False`.
10. **Interfaz:** se creó la clase base `StepsTab` (panel de pasos + gráfica) que ahora comparten *Primer orden*, *Orden superior*, *Laplace* y *Sistemas*. Se agregaron el retrato de fase interactivo y las series de tiempo, con un corte automático cuando la solución diverge.
11. **Pruebas** de los 52 ejemplos (`test_laplace.py` y `test_gui.py` ampliado) y recompilación del `.exe`.

---

## 6. Arquitectura del código

```
diffeq_studio/
├── app.py            Interfaz gráfica (pestañas, gráficas, hilos)
├── solvers.py        Motor simbólico: parseo + métodos paso a paso (EDO de una ecuación)
├── laplace.py        Transformada de Laplace: directa, inversa, PVI, convolución
├── systems.py        Sistemas de EDO: valores propios, variación de parámetros, linealización
├── numerics.py       Euler / Heun / RK4 (escalar y sistemas), campo de direcciones
├── examples.py       Ejemplos precargados y valores por defecto
├── test_solvers.py   Prueba del motor con todos los ejemplos
├── test_laplace.py   Prueba de Laplace y sistemas
├── test_gui.py       Prueba de humo de la interfaz (guarda las gráficas)
├── show_steps.py     Imprime en consola los pasos de una ecuación
├── README.md         Resumen corto
├── DOCUMENTACION.md  Este documento
└── dist/
    └── DiffEqStudio.exe
```

### 6.1 Flujo de una resolución

```
 Usuario escribe "y'' + y = tan(x)"  +  "y(0)=1, y'(0)=0"
        │
        ▼
 app.SolverTab.solve() ──► App.run()  (hilo secundario)
        │
        ▼
 solvers.solve(texto, condiciones)
   ├─ parse_equation()  → normaliza y', y'' → símbolos YP1, YP2 y parsea con SymPy
   ├─ parse_ics()       → [(orden, x0, valor)]
   ├─ _first_order()  ó  _linear_cc()  ó  _exact()   → llena Result.steps
   ├─ _fallback()       → dsolve si ningún método aplica
   ├─ verificación con checkodesol
   └─ _apply_ics()      → constantes y solución particular
        │
        ▼  (cola → hilo de la interfaz)
 SolverTab.show(result)  → escribe pasos con formato
 SolverTab.draw(result)  → gráfica con numerics.py
```

### 6.2 `solvers.py` (motor)

| Elemento | Descripción |
|---|---|
| `parse_equation(text)` | Convierte `y'`, `y''`, `dy/dx` en símbolos internos `YP1`, `YP2`… y parsea con multiplicación implícita (`2x` → `2*x`). Detecta la forma `M dx + N dy`. Devuelve `('ode', (expr, orden))` o `('exact', (M, N))`. |
| `parse_ics(text)` | Expresión regular para `y'(0)=2`. Los valores se convierten a fracciones exactas (`0.5` → `1/2`). |
| `parse_rhs(text)` | Obtiene F(x, y) para los módulos numéricos. Acepta `x + y` o `y' = x + y`. |
| `Result` | Guarda los pasos (`title`, `text`, `math`, `ok`, `warn`), el método, F, las soluciones generales, la solución implícita y la particular. |
| `_separable`, `_linear`, `_bernoulli`, `_homogeneous` | Métodos de primer orden. |
| `_linear_steps(P, Q, var)` | Factor integrante reutilizable (también lo usa Bernoulli con la variable v). |
| `_exact(M, N)` | Exactas y factor integrante. |
| `_linear_cc(expr, n)` | Coeficientes constantes y variación de parámetros. |
| `_solve_implicit(G)` | Intenta despejar y de G(x, y) = C y renombra e^C → C. |
| `_fallback(ode)` | `classify_ode` + `dsolve`, con los nombres de los métodos traducidos. |
| `_apply_ics` | Problema de valor inicial. |
| `solve(text, ics)` | Punto de entrada: orquesta todo lo anterior. |

### 6.3 `laplace.py`

| Elemento | Descripción |
|---|---|
| `t` (positivo), `s`, `τ` | Símbolos. Con t positivo, SymPy no añade u(t) a las inversas. |
| `parse_t`, `parse_s` | Lectura de f(t) y F(s). Reconoce `u()`/`H()` como escalón y `delta()` como impulso. |
| `forward_steps(f, res)` | Linealidad, regla aplicada a cada término y resultado con denominador común. |
| `inverse_steps(F, res)` | Agrupa por e^(−as), aplica fracciones parciales y la regla de cada fracción, y desplaza con u(t−a). |
| `transform`, `inverse`, `ivp`, `convolution` | Funciones de cada modo de la pestaña. Devuelven un `Result` con `res.plot`, la lista de funciones a graficar. |

### 6.4 `systems.py`

| Elemento | Descripción |
|---|---|
| `parse_system(text)` | Lee líneas `x' = …` o `dx/dt = …` y convierte decimales a fracciones. |
| `classify(eigs, A)` | Tipo y estabilidad del punto crítico 2×2. |
| `_linear` | Valores y vectores propios con todos los casos. |
| `_expm` | Respaldo con la matriz exponencial. |
| `_assemble` | Solución homogénea, variación de parámetros, verificación, PVI y clasificación. |
| `_nonlinear` | Equilibrios, jacobiano y linealización. |
| `solve(text, ics)` | Punto de entrada. Deja en `res.sys` los datos para el retrato de fase. |

### 6.5 `numerics.py`

| Función | Descripción |
|---|---|
| `make_f(F)` | Convierte F simbólica en función NumPy vectorizada (`lambdify`) y silencia divisiones por cero. |
| `eval_expr(expr, xs)` | Evalúa y(x). Devuelve NaN donde el valor no es real. Si `lambdify` falla, recurre a `evalf`. |
| `euler`, `heun`, `rk4` | Métodos numéricos de paso fijo. |
| `curve(...)` | Trayectoria RK4 que se detiene al diverger. |
| `slope_field(...)` | Campo de pendientes normalizado + isoclina y' = 0. |
| `make_system`, `rk4_system` | RK4 para sistemas (hacia adelante o atrás); se detiene al diverger. |

### 6.6 `app.py` (interfaz)

| Clase | Descripción |
|---|---|
| `App` | Ventana principal, estilos y pestañas. `run()` ejecuta cálculos en hilos y `_poll()` entrega los resultados al hilo de Tkinter mediante una cola. |
| `PlotPanel` | Figura de Matplotlib con su barra de herramientas. |
| `StepsTab` | Base común: panel de pasos con formato, gráfica, `solve()` en segundo plano y manejo de errores. Cada pestaña define `job()` y `draw()`. |
| `SolverTab` | Pestañas *Primer orden* y *Orden superior* (la misma clase con distintos ejemplos). |
| `LaplaceTab` | Cuatro modos con campos y ejemplos que cambian según el modo. |
| `SystemsTab` | Sistema multilínea, retrato de fase interactivo y series de tiempo. |
| `NumericTab` | Tabla y gráfica de los métodos numéricos. |
| `FieldTab` | Campo de direcciones interactivo. |

---

## 7. Decisiones de diseño importantes

- **Pasos programados a mano en vez de usar solo `dsolve`.** Es lo que da valor académico a la app. `dsolve` se usa solo como respaldo y para la solución exacta en métodos numéricos.
- **Orden de detección en primer orden.** Primero va separable, porque es el método más simple de explicar; luego lineal, Bernoulli y homogénea. Una ecuación puede ser de varios tipos y la app muestra el más sencillo.
- **Variación de parámetros en lugar de coeficientes indeterminados.** Es un procedimiento general que sirve para cualquier g(x).
- **ln|u| en la variación de parámetros.** Al integrar aparecen términos como ln(sin x − 1), que son complejos para todo x real. Se reemplazan por ln|u|, que tiene la misma derivada; así la solución es real y se puede graficar. La verificación se hace sin el valor absoluto.
- **Fracciones exactas en condiciones iniciales.** `y(0)=0.5` se convierte en 1/2 para que el resultado sea simbólico y limpio, no `0.6666…·e^{2x}`.
- **Renombrado de constantes.** e^{C₁} se reemplaza por C₁, como se hace en clase.
- **Asíntotas verticales.** Los tramos que salen muy fuera de la ventana se convierten en NaN, para que Matplotlib no una ramas separadas por una asíntota.
- **Sistemas escritos ecuación por ecuación** (no como matriz), igual que en el cuaderno. La app extrae A y F(t) derivando cada lado derecho, y así distingue sola entre sistemas lineales y no lineales.
- **Vectores propios sin denominadores.** Se multiplican por el mínimo común múltiplo para que queden como en clase: (1, 2) y no (1/2, 1).
- **Retrato de fase centrado** en los equilibrios y la condición inicial, con un tamaño mínimo de ±3. Cuando la solución diverge, la ventana de tiempo se acorta para que se vea su forma.
- **Cálculo en un hilo secundario.** Algunas integrales tardan, y así la ventana no se congela. Tkinter no es seguro entre hilos, por eso los resultados vuelven por una `queue.Queue` revisada con `after()`.

---

## 8. Pruebas

```bash
python test_solvers.py         # 19 ejemplos de primer orden y orden superior; -v muestra los pasos
python test_laplace.py         # 23 ejemplos de Laplace + 10 de sistemas
python test_gui.py carpeta     # abre la app, resuelve los 52 ejemplos y guarda las gráficas
python show_steps.py "y' = x*y" "y(0)=1"   # pasos en consola de una ecuación
```

Resultado actual: los **52 ejemplos** se resuelven sin errores, tanto en el motor como en la interfaz. Todas las soluciones pasan su verificación. Los ejemplos cubren:
- **Primer orden:** separables, lineales, Bernoulli, homogénea, exacta, exacta con μ(x) y logística.
- **Orden superior:** raíces reales, repetidas y complejas, variación de parámetros, tercer orden y Cauchy-Euler.
- **Laplace:** tabla, ambos teoremas de traslación, escalón, impulso, fracciones parciales, PVI y convolución.
- **Sistemas:** silla, nodo, foco, centro, nodo impropio, sistema no homogéneo, 3×3, péndulo amortiguado, especies en competencia y Lotka-Volterra.

---

## 9. Compilar el ejecutable

```bash
pip install sympy numpy matplotlib pyinstaller
python app.py                                                     # ejecutar sin compilar
python -m PyInstaller --onefile --windowed --name DiffEqStudio app.py
```

- `--onefile` genera un único `.exe` en `dist/`.
- `--windowed` evita que se abra una consola negra.
- El ejecutable pesa unos 48 MB porque incluye Python, SymPy, NumPy y Matplotlib.

---

## 10. Limitaciones conocidas

- Las ecuaciones de orden superior **sin coeficientes constantes** (Cauchy-Euler, coeficientes variables, series) se resuelven, pero **sin explicación paso a paso**.
- Las ecuaciones de primer orden que no encajan en los cinco métodos (Riccati, etc.) también usan el respaldo de SymPy.
- Algunas integrales no tienen forma cerrada. En ese caso la solución queda expresada con ∫ y la gráfica puede tardar más.
- Si una solución no se puede despejar, queda **implícita** y la gráfica usa RK4 en lugar de la fórmula.
- `dsolve` puede tardar mucho con ecuaciones muy complicadas. La interfaz sigue respondiendo, pero el resultado demora.
- Las ecuaciones de las pestañas *Primer orden* y *Orden superior* deben escribirse en términos de x e y. En *Laplace* se usa t (también se acepta x).
- El PVI con Laplace requiere coeficientes constantes y condiciones en t = 0.
- En sistemas, el retrato de fase solo se dibuja para sistemas 2×2 autónomos. Para 3×3 o no autónomos se grafican las soluciones en el tiempo.
- Los equilibrios de los sistemas no lineales se buscan de forma exacta. Si `solve` no los encuentra (ecuaciones trascendentes complicadas), solo se muestra el retrato numérico.

---

## 11. Trabajo futuro

Módulos planeados para completar el temario de Cálculo 4:

Ya están hechos: ~~Transformada de Laplace~~ y ~~Sistemas de EDO con retrato de fase~~ (segunda parte).

Pendientes:

1. **Cauchy-Euler paso a paso** y **reducción de orden**.
2. **Series de potencias y Frobenius.**
3. **Aplicaciones con animación:** masa-resorte, circuitos RLC, mezclas, enfriamiento de Newton y crecimiento poblacional.
4. **Coeficientes indeterminados** como método alternativo a la variación de parámetros.
5. **Modo práctica:** ejercicios aleatorios con autocorrección.
6. **Exportar el procedimiento a PDF.**
7. **Series de Fourier y EDP** (calor y onda) por separación de variables.
