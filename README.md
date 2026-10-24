# DiffEq Studio — versión mínima

App de escritorio para Cálculo 4 (ecuaciones diferenciales) que resuelve paso a paso y grafica.

## Módulos
- **Primer orden**: separables, lineales (factor integrante), Bernoulli, homogéneas, exactas (con factor integrante μ(x)/μ(y)).
- **Orden superior**: coeficientes constantes de cualquier orden (ecuación característica, Wronskiano, variación de parámetros). Otros tipos (Cauchy-Euler, etc.) se resuelven con SymPy sin pasos.
- **Laplace**: transformada directa e inversa paso a paso, PVI (con escalón e impulso) y convolución.
- **Sistemas**: valores/vectores propios, variación de parámetros, clasificación del punto crítico, linealización y retrato de fase interactivo.
- **Métodos numéricos**: Euler, Heun y RK4 con tabla y error frente a la solución exacta.
- **Campo de direcciones**: clic sobre la gráfica para trazar la curva solución.

Todas las soluciones se verifican sustituyéndolas en la ecuación.

## Archivos
| Archivo | Contenido |
|---|---|
| `app.py` | Interfaz (Tkinter + Matplotlib) |
| `solvers.py` | Motor simbólico paso a paso (SymPy) |
| `numerics.py` | Euler / Heun / RK4 y campo de direcciones |
| `examples.py` | Ejemplos precargados |
| `test_solvers.py`, `test_gui.py`, `show_steps.py` | Pruebas |

## Ejecutar y compilar
```
pip install sympy numpy matplotlib pyinstaller
python app.py
python -m PyInstaller --onefile --windowed --name DiffEqStudio app.py
```
