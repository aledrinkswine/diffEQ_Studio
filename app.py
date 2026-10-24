"""DiffEq Studio — interfaz gráfica (Tkinter + Matplotlib)."""
import queue
import threading
import tkinter as tk
from tkinter import ttk

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import sympy as sp
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import laplace as LP
import numerics as NM
import solvers as S
import systems as SY
from examples import (CONVOLUTION, FIELD_DEFAULT, FIRST_ORDER, HIGHER_ORDER, INVERSE,
                      IVP_LAPLACE, NUMERIC_DEFAULT, SYSTEMS, TRANSFORM)

ACCENT = '#2b6cb0'
BG = '#f4f6fa'
COLORS = ['#2b6cb0', '#c53030', '#2f855a', '#d69e2e', '#805ad5', '#dd6b20']

HELP = """\
CÓMO ESCRIBIR LAS ECUACIONES

  Derivadas:      y'   y''   y'''    (también dy/dx)
  Productos:      2y   3x*y   x^2 y''      (* opcional entre número y letra)
  Potencias:      y^3   x^(1/2)   e^(2x)
  Funciones:      sin cos tan (o sen, tg)  exp  ln  sqrt  arctan
  Forma M dx + N dy:   (2x*y + 3) dx + (x^2 - 1) dy = 0

  Condiciones iniciales:   y(0)=1, y'(0)=0

MÓDULOS

  Primer orden ........ separables, lineales (factor integrante), Bernoulli,
                        homogéneas y exactas (con factor integrante μ(x) o μ(y)).
  Orden superior ...... coeficientes constantes de cualquier orden: ecuación
                        característica, raíces reales/repetidas/complejas,
                        Wronskiano y variación de parámetros.
                        Otros tipos (Cauchy-Euler, etc.) se resuelven con SymPy.
  Laplace ............. transformada directa (tabla + teoremas de traslación),
                        inversa con fracciones parciales, PVI y convolución.
                        Escalón: u(t-a)   Impulso: delta(t-a)
  Sistemas ............ una ecuación por línea:  x' = x + 2y
                                                 y' = 3x + 2y
                        Lineales: valores/vectores propios (reales, complejos,
                        repetidos), variación de parámetros, clasificación del
                        punto crítico y retrato de fase (clic = trayectoria).
                        No lineales: equilibrios, jacobiano y linealización.
  Métodos numéricos ... Euler, Euler mejorado (Heun) y RK4 con tabla y errores
                        frente a la solución exacta.
  Campo de direcciones  haz clic en la gráfica para trazar la curva solución
                        que pasa por ese punto.

  Todas las soluciones se verifican sustituyéndolas en la ecuación original.
"""


class PlotPanel(ttk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.fig = Figure(figsize=(6, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        self.toolbar.pack(side='bottom', fill='x')
        self.canvas.get_tk_widget().pack(fill='both', expand=True)

    def finish(self, title='', ax=None, xlabel='x', ylabel='y'):
        ax = ax or self.ax
        ax.axhline(0, color='#999', lw=0.8)
        ax.axvline(0, color='#999', lw=0.8)
        ax.grid(alpha=0.3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title, fontsize=10)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8, loc='best')
        self.fig.tight_layout()
        self.canvas.draw_idle()


def _cut(ys, lo, hi):
    """Oculta los tramos fuera de la vista para no unir las asíntotas verticales."""
    span = hi - lo
    return np.where((ys < lo - span) | (ys > hi + span), np.nan, ys)


def _robust_ylim(ax, lines):
    vals = np.concatenate([l.get_ydata()[np.isfinite(l.get_ydata())] for l in lines])
    if vals.size == 0:
        return
    lo, hi = np.percentile(vals, [2, 98])
    if hi - lo < 1e-9:
        lo, hi = lo - 1, hi + 1
    pad = 0.15 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    ax.set_ylim(lo, hi)
    for l in lines:
        l.set_ydata(_cut(l.get_ydata(), lo, hi))


class StepsTab(ttk.Frame):
    """Base: procedimiento paso a paso a la izquierda y gráfica a la derecha."""

    def build_body(self):
        body = ttk.PanedWindow(self, orient='horizontal')
        body.pack(fill='both', expand=True, pady=(10, 0))
        left = ttk.Frame(body)
        self.text = tk.Text(left, wrap='none', font=('Consolas', 11), bg='white',
                            relief='flat', padx=10, pady=8, state='disabled')
        sy = ttk.Scrollbar(left, orient='vertical', command=self.text.yview)
        sx = ttk.Scrollbar(left, orient='horizontal', command=self.text.xview)
        self.text.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side='right', fill='y')
        sx.pack(side='bottom', fill='x')
        self.text.pack(fill='both', expand=True)
        self.text.tag_configure('method', font=('Segoe UI', 13, 'bold'), foreground=ACCENT)
        self.text.tag_configure('title', font=('Segoe UI', 11, 'bold'), foreground='#1a202c',
                                spacing1=8, spacing3=4)
        self.text.tag_configure('ok', foreground='#2f855a', font=('Segoe UI', 10, 'bold'))
        self.text.tag_configure('warn', foreground='#c05621', font=('Segoe UI', 10, 'bold'))
        self.text.tag_configure('error', foreground='#c53030', font=('Segoe UI', 11, 'bold'))
        self.text.tag_configure('math', foreground='#1a365d', lmargin1=24, lmargin2=24)
        body.add(left, weight=1)
        self.plot = PlotPanel(body)
        body.add(self.plot, weight=1)

    def _write(self, chunks):
        t = self.text
        t.config(state='normal')
        t.delete('1.0', 'end')
        for content, tag in chunks:
            t.insert('end', content, tag)
        t.config(state='disabled')

    def fail(self, e):
        self.btn.state(['!disabled'])
        self.status.config(text='')
        msg = str(e) if isinstance(e, S.ParseError) else f"No se pudo resolver: {e}"
        self._write([(msg + '\n', 'error')])

    def show(self, res):
        self.btn.state(['!disabled'])
        self.status.config(text='')
        chunks = [(f"Método: {res.method}\n", 'method')]
        for kind, c in res.steps:
            if kind == 'math':
                chunks.append((S.pretty(c) + '\n\n', 'math'))
            elif kind == 'title':
                chunks.append((f"\n{c}\n", 'title'))
            else:
                chunks.append((c + '\n', kind if kind in ('ok', 'warn') else None))
        self._write(chunks)
        try:
            self.draw(res)
        except Exception as e:
            self.plot.fig.clear()
            self.plot.ax = self.plot.fig.add_subplot(111)
            self.plot.ax.text(0.5, 0.5, f"No se pudo graficar:\n{e}", ha='center',
                              transform=self.plot.ax.transAxes)
            self.plot.canvas.draw_idle()

    def solve(self):
        self.btn.state(['disabled'])
        self.status.config(text='Resolviendo…')
        self.app.run(self.job(), self.show, self.fail)


class SolverTab(StepsTab):
    def __init__(self, master, app, examples, placeholder):
        super().__init__(master, padding=10)
        self.app = app
        self.examples = examples

        top = ttk.Frame(self)
        top.pack(fill='x')
        ttk.Label(top, text='Ecuación:').grid(row=0, column=0, sticky='w')
        self.eq = ttk.Entry(top, font=('Consolas', 13))
        self.eq.grid(row=0, column=1, sticky='ew', padx=6)
        self.eq.insert(0, placeholder[0])
        ttk.Label(top, text='Condiciones iniciales (opcional):').grid(row=1, column=0, sticky='w', pady=(6, 0))
        self.ics = ttk.Entry(top, font=('Consolas', 13))
        self.ics.grid(row=1, column=1, sticky='ew', padx=6, pady=(6, 0))
        self.ics.insert(0, placeholder[1])
        self.btn = ttk.Button(top, text='Resolver ▶', style='Accent.TButton', command=self.solve)
        self.btn.grid(row=0, column=2, rowspan=2, sticky='ns', padx=(4, 0))
        ttk.Label(top, text='Ejemplos:').grid(row=2, column=0, sticky='w', pady=(6, 0))
        self.combo = ttk.Combobox(top, state='readonly', font=('Consolas', 11),
                                  values=[f"{e}    {i}" for e, i in examples])
        self.combo.grid(row=2, column=1, sticky='ew', padx=6, pady=(6, 0))
        self.combo.bind('<<ComboboxSelected>>', self.load_example)
        self.status = ttk.Label(top, text='', foreground=ACCENT)
        self.status.grid(row=2, column=2, sticky='w')
        top.columnconfigure(1, weight=1)
        for w in (self.eq, self.ics):
            w.bind('<Return>', lambda e: self.solve())

        self.build_body()

    def load_example(self, _=None):
        eq, ics = self.examples[self.combo.current()]
        self.eq.delete(0, 'end')
        self.eq.insert(0, eq)
        self.ics.delete(0, 'end')
        self.ics.insert(0, ics)
        self.solve()

    def job(self):
        eq, ics = self.eq.get(), self.ics.get()
        return lambda: S.solve(eq, ics)


    def draw(self, res):
        ax = self.plot.ax
        ax.clear()
        ics = [(k, float(a), float(b)) for k, a, b in res.ics]
        x0, y0 = (ics[0][1], ics[0][2]) if ics else (0.0, 0.0)
        if res.order == 1 and res.F is not None:
            xr, yr = (x0 - 5, x0 + 5), (y0 - 5, y0 + 5)
            f = NM.make_f(res.F)
            NM.slope_field(ax, f, xr, yr)
            if res.particular is not None:
                xs = np.linspace(*xr, 600)
                ax.plot(xs, _cut(NM.eval_expr(res.particular, xs), *yr), color=COLORS[1], lw=2.2,
                        label='Solución particular')
            if ics:
                if res.particular is None:
                    for xe in xr:
                        xs, ys = NM.curve(f, x0, y0, xe, 1e3)
                        ax.plot(xs, ys, color=COLORS[1], lw=2.2,
                                label='Solución (RK4)' if xe == xr[1] else None)
                ax.plot([x0], [y0], 'o', color=COLORS[1])
            else:
                for c in np.linspace(yr[0] + 0.5, yr[1] - 0.5, 9):
                    for xe in xr:
                        xs, ys = NM.curve(f, x0, c, xe, 1e3)
                        ax.plot(xs, ys, color=COLORS[0], lw=1.1, alpha=0.8)
            ax.set_xlim(xr)
            ax.set_ylim(yr)
            self.plot.finish("Campo de direcciones y soluciones")
            return
        xs = np.linspace(x0, x0 + 10, 800)
        lines = []
        if res.particular is not None:
            ys = NM.eval_expr(res.particular, xs)
            lines += ax.plot(xs, ys, color=COLORS[1], lw=2.2, label='Solución particular')
        elif res.general:
            s = res.general[0]
            consts = sorted((c for c in s.free_symbols if c.name.startswith('C')), key=str)
            combos = [(1, 0), (0, 1), (1, 1), (1, -1), (-1, 0.5)]
            for i, combo in enumerate(combos):
                vals = {c: (combo[j] if j < len(combo) else 0) for j, c in enumerate(consts)}
                ys = NM.eval_expr(s.subs(vals), xs)
                lbl = ', '.join(f"{c}={v}" for c, v in vals.items())
                lines += ax.plot(xs, ys, color=COLORS[i % len(COLORS)], lw=1.6, label=lbl)
        else:
            ax.text(0.5, 0.5, 'Sin solución explícita para graficar', ha='center',
                    transform=ax.transAxes)
        if lines:
            _robust_ylim(ax, lines)
        self.plot.finish("Solución particular" if res.particular is not None
                         else "Familia de soluciones")


LAPLACE_MODES = [
    ('transform', 'Transformada ℒ{f(t)}', 'f(t) =', None, [(f, '') for f in TRANSFORM]),
    ('inverse', 'Inversa ℒ⁻¹{F(s)}', 'F(s) =', None, [(F, '') for F in INVERSE]),
    ('ivp', 'PVI con Laplace', 'Ecuación:', 'Condiciones en t = 0:', IVP_LAPLACE),
    ('conv', 'Convolución f ∗ g', 'f(t) =', 'g(t) =', CONVOLUTION),
]


class LaplaceTab(StepsTab):
    def __init__(self, master, app):
        super().__init__(master, padding=10)
        self.app = app
        self.mode = tk.StringVar(value='transform')
        modes = ttk.Frame(self)
        modes.pack(fill='x')
        for key, label, *_ in LAPLACE_MODES:
            ttk.Radiobutton(modes, text=label, value=key, variable=self.mode,
                            command=self.set_mode).pack(side='left', padx=(0, 16))

        top = ttk.Frame(self)
        top.pack(fill='x', pady=(8, 0))
        self.l1 = ttk.Label(top, width=20)
        self.l1.grid(row=0, column=0, sticky='w')
        self.e1 = ttk.Entry(top, font=('Consolas', 13))
        self.e1.grid(row=0, column=1, sticky='ew', padx=6)
        self.l2 = ttk.Label(top, width=20)
        self.l2.grid(row=1, column=0, sticky='w', pady=(6, 0))
        self.e2 = ttk.Entry(top, font=('Consolas', 13))
        self.e2.grid(row=1, column=1, sticky='ew', padx=6, pady=(6, 0))
        self.btn = ttk.Button(top, text='Resolver ▶', style='Accent.TButton', command=self.solve)
        self.btn.grid(row=0, column=2, rowspan=2, sticky='ns', padx=(4, 0))
        ttk.Label(top, text='Ejemplos:').grid(row=2, column=0, sticky='w', pady=(6, 0))
        self.combo = ttk.Combobox(top, state='readonly', font=('Consolas', 11))
        self.combo.grid(row=2, column=1, sticky='ew', padx=6, pady=(6, 0))
        self.combo.bind('<<ComboboxSelected>>', self.load_example)
        self.status = ttk.Label(top, text='', foreground=ACCENT)
        self.status.grid(row=2, column=2, sticky='w')
        top.columnconfigure(1, weight=1)
        for w in (self.e1, self.e2):
            w.bind('<Return>', lambda e: self.solve())
        self.build_body()
        self.set_mode()

    def _cfg(self):
        return next(m for m in LAPLACE_MODES if m[0] == self.mode.get())

    def _fill(self, a, b):
        for e, v in ((self.e1, a), (self.e2, b)):
            e.delete(0, 'end')
            e.insert(0, v)

    def set_mode(self):
        _, _, l1, l2, ex = self._cfg()
        self.l1.config(text=l1)
        if l2:
            self.l2.config(text=l2)
            self.l2.grid()
            self.e2.grid()
        else:
            self.l2.grid_remove()
            self.e2.grid_remove()
        self.combo.config(values=[f"{a}    {b}" for a, b in ex])
        self.combo.set('')
        self._fill(*ex[0])

    def load_example(self, _=None):
        self._fill(*self._cfg()[4][self.combo.current()])
        self.solve()

    def job(self):
        a, b = self.e1.get(), self.e2.get()
        return {'transform': lambda: LP.transform(a), 'inverse': lambda: LP.inverse(a),
                'ivp': lambda: LP.ivp(a, b), 'conv': lambda: LP.convolution(a, b)}[self.mode.get()]

    def draw(self, res):
        ax = self.plot.ax
        ax.clear()
        ts = np.linspace(0, 10, 1000)
        lines, impulses = [], False
        for i, (lbl, e) in enumerate(res.plot):
            color = COLORS[i % len(COLORS)]
            deltas = e.atoms(sp.DiracDelta)
            for d in deltas:
                a = float(sp.solve(d.args[0], LP.t)[0])
                ax.annotate('', xy=(a, 1), xytext=(a, 0),
                            arrowprops=dict(arrowstyle='-|>', color=color, lw=2))
                ax.plot([], [], color=color, label=f"{lbl}: impulso en t = {a:g}")
                impulses = True
            e = e.subs({d: 0 for d in deltas})
            if e == 0:
                continue
            lines += ax.plot(ts, NM.eval_expr(e, ts, LP.t), color=color, lw=2, label=lbl)
        if lines:
            _robust_ylim(ax, lines)
        elif impulses:
            ax.set_ylim(-0.2, 1.3)
        ax.set_xlim(0, 10)
        self.plot.finish(self._cfg()[1].split(' ', 1)[0], xlabel='t', ylabel='')


class SystemsTab(StepsTab):
    def __init__(self, master, app):
        super().__init__(master, padding=10)
        self.app = app
        top = ttk.Frame(self)
        top.pack(fill='x')
        ttk.Label(top, text='Sistema\n(una ecuación por línea):').grid(row=0, column=0, sticky='nw')
        self.eq = tk.Text(top, height=3, font=('Consolas', 13), relief='solid', bd=1, wrap='none')
        self.eq.grid(row=0, column=1, sticky='ew', padx=6)
        self.eq.insert('1.0', SYSTEMS[0][0])
        ttk.Label(top, text='Condiciones iniciales (opcional):').grid(row=1, column=0, sticky='w', pady=(6, 0))
        self.ics = ttk.Entry(top, font=('Consolas', 13))
        self.ics.grid(row=1, column=1, sticky='ew', padx=6, pady=(6, 0))
        self.ics.insert(0, SYSTEMS[0][1])
        self.btn = ttk.Button(top, text='Resolver ▶', style='Accent.TButton', command=self.solve)
        self.btn.grid(row=0, column=2, rowspan=2, sticky='ns', padx=(4, 0))
        ttk.Label(top, text='Ejemplos:').grid(row=2, column=0, sticky='w', pady=(6, 0))
        self.combo = ttk.Combobox(top, state='readonly', font=('Consolas', 11),
                                  values=[f"{e.replace(chr(10), ' ;  ')}      {i}" for e, i in SYSTEMS])
        self.combo.grid(row=2, column=1, sticky='ew', padx=6, pady=(6, 0))
        self.combo.bind('<<ComboboxSelected>>', self.load_example)
        self.status = ttk.Label(top, text='', foreground=ACCENT)
        self.status.grid(row=2, column=2, sticky='w')
        top.columnconfigure(1, weight=1)
        self.ics.bind('<Return>', lambda e: self.solve())
        self.eq.bind('<Control-Return>', lambda e: (self.solve(), 'break')[1])
        self.build_body()
        self.plot.canvas.mpl_connect('button_press_event', self.on_click)
        self.phase = None

    def load_example(self, _=None):
        eq, ics = SYSTEMS[self.combo.current()]
        self.eq.delete('1.0', 'end')
        self.eq.insert('1.0', eq)
        self.ics.delete(0, 'end')
        self.ics.insert(0, ics)
        self.solve()

    def job(self):
        text, ics = self.eq.get('1.0', 'end'), self.ics.get()
        return lambda: SY.solve(text, ics)

    def draw(self, res):
        d = res.sys
        fig = self.plot.fig
        fig.clear()
        self.phase = None
        names = d['names']
        f = NM.make_system(d['rhs'], d['vars'], SY.t)
        x0, t0 = None, 0.0
        if all(nm in d['ics'] for nm in names):
            t0 = float(d['ics'][names[0]][0])
            x0 = np.array([float(d['ics'][nm][1]) for nm in names])
        phase = len(names) == 2 and not any(r.has(SY.t) for r in d['rhs'])
        if phase and x0 is not None:
            gs = fig.add_gridspec(2, 1, height_ratios=[3, 2])
            axp, axt = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        elif phase:
            axp, axt = fig.add_subplot(111), None
        else:
            axp, axt = None, fig.add_subplot(111)
        self.plot.ax = axp or axt
        if axp is not None:
            self._phase(axp, d, f, x0)
        if axt is not None:
            self._time(axt, d, f, x0, t0)

    def _phase(self, ax, d, f, x0):
        pts = list(d['equilibria']) + ([tuple(x0)] if x0 is not None else []) or [(0.0, 0.0)]
        c = np.mean(pts, axis=0)
        R = max(3.0, 1.4 * max(np.hypot(p[0] - c[0], p[1] - c[1]) for p in pts))
        xr, yr = (c[0] - R, c[0] + R), (c[1] - R, c[1] + R)
        X, Yg = np.meshgrid(np.linspace(*xr, 30), np.linspace(*yr, 30))
        U, V = np.nan_to_num(f(0.0, np.array([X, Yg])))
        ax.streamplot(X, Yg, U, V, color='#a0aec0', density=1.3, linewidth=0.8, arrowsize=1)
        Xf, Yf = np.meshgrid(np.linspace(*xr, 200), np.linspace(*yr, 200))
        Uf, Vf = f(0.0, np.array([Xf, Yf]))
        for comp, color, nm in ((Uf, '#d69e2e', d['names'][0]), (Vf, '#2f855a', d['names'][1])):
            if np.nanmin(comp) < 0 < np.nanmax(comp):
                ax.contour(Xf, Yf, comp, [0], colors=color, linewidths=1, linestyles='--')
                ax.plot([], [], '--', color=color, lw=1, label=f"Nulclina {nm}' = 0")
        e0 = d['equilibria'][0] if d['equilibria'] else (0.0, 0.0)
        for i, (val, v) in enumerate(d['dirs']):
            vv = np.array([float(v[0]), float(v[1])])
            vv /= np.linalg.norm(vv)
            ax.plot([e0[0] - 2 * R * vv[0], e0[0] + 2 * R * vv[0]],
                    [e0[1] - 2 * R * vv[1], e0[1] + 2 * R * vv[1]],
                    color=COLORS[4 + i % 2], lw=1.6, label=f"Vector propio (λ = {val})")
        for e in d['equilibria']:
            ax.plot(e[0], e[1], 'o', color='k', ms=7, zorder=5)
        if x0 is not None:
            _, Xs = NM.rk4_system(f, 0.0, x0, 20.0, 3000, 50 * R)
            ax.plot(Xs[:, 0], Xs[:, 1], color=COLORS[1], lw=2.2, label='Trayectoria con X(0)')
            ax.plot(*x0, 'o', color=COLORS[1], zorder=6)
        ax.set_xlim(xr)
        ax.set_ylim(yr)
        self.phase = (ax, f, xr, yr, R)
        title = f"Retrato de fase — {d['class']}" if d['class'] else "Retrato de fase (clic = trayectoria)"
        self.plot.finish(title, ax, d['names'][0], d['names'][1])

    def _time(self, ax, d, f, x0, t0):
        if x0 is None:
            ax.text(0.5, 0.5, 'Escribe condiciones iniciales\npara graficar las soluciones en el tiempo.',
                    ha='center', va='center', transform=ax.transAxes)
            self.plot.finish('', ax, 't', '')
            return
        ts = np.linspace(t0, t0 + 10, 800)
        if d['particular'] is not None:
            series = [NM.eval_expr(comp, ts, SY.t) for comp in d['particular']]
            how = 'solución exacta'
        else:
            ts, Xs = NM.rk4_system(f, t0, x0, 10.0, 2000)
            series = list(Xs.T)
            how = 'RK4'
        # si la solución diverge, se corta la ventana de tiempo para que se vea su forma
        big = np.any([np.abs(ys) > 50 * max(1.0, np.abs(x0).max()) for ys in series], axis=0)
        if big.any() and np.argmax(big) > 20:
            k = np.argmax(big)
            ts, series = ts[:k], [ys[:k] for ys in series]
        lines = []
        for i, (nm, ys) in enumerate(zip(d['names'], series)):
            lines += ax.plot(ts, ys, color=COLORS[i % len(COLORS)], lw=2, label=f"{nm}(t)")
        _robust_ylim(ax, lines)
        self.plot.finish(f"Soluciones en el tiempo ({how})", ax, 't', '')

    def on_click(self, ev):
        if not self.phase or ev.inaxes is not self.phase[0] or self.plot.toolbar.mode:
            return
        ax, f, xr, yr, R = self.phase
        color = COLORS[np.random.randint(len(COLORS))]
        for T in (20.0, -20.0):
            _, Xs = NM.rk4_system(f, 0.0, np.array([ev.xdata, ev.ydata]), T, 3000, 50 * R)
            ax.plot(Xs[:, 0], Xs[:, 1], color=color, lw=1.8)
        ax.plot(ev.xdata, ev.ydata, 'o', color=color)
        ax.set_xlim(xr)
        ax.set_ylim(yr)
        self.plot.canvas.draw_idle()


class NumericTab(ttk.Frame):
    def __init__(self, master, app):
        super().__init__(master, padding=10)
        self.app = app
        top = ttk.Frame(self)
        top.pack(fill='x')
        self.vars = {}
        labels = [("y' =", 34), ('x₀', 7), ('y₀', 7), ('h', 7), ('x final', 7)]
        for i, ((lbl, w), val) in enumerate(zip(labels, NUMERIC_DEFAULT)):
            ttk.Label(top, text=lbl).pack(side='left', padx=(10 if i else 0, 4))
            e = ttk.Entry(top, width=w, font=('Consolas', 12))
            e.insert(0, val)
            e.pack(side='left')
            e.bind('<Return>', lambda ev: self.compute())
            self.vars[lbl] = e
        self.btn = ttk.Button(top, text='Calcular ▶', style='Accent.TButton', command=self.compute)
        self.btn.pack(side='left', padx=10)
        self.status = ttk.Label(self, text='', foreground=ACCENT)
        self.status.pack(fill='x', pady=(6, 0))

        body = ttk.PanedWindow(self, orient='horizontal')
        body.pack(fill='both', expand=True, pady=(8, 0))
        left = ttk.Frame(body)
        cols = ('n', 'x', 'Euler', 'Heun', 'RK4', 'Exacta', 'Error Euler', 'Error RK4')
        self.tree = ttk.Treeview(left, columns=cols, show='headings')
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=40 if c == 'n' else 95, anchor='e')
        sb = ttk.Scrollbar(left, orient='vertical', command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        self.tree.pack(fill='both', expand=True)
        body.add(left, weight=1)
        self.plot = PlotPanel(body)
        body.add(self.plot, weight=1)

    def compute(self):
        v = {k: e.get() for k, e in self.vars.items()}
        self.btn.state(['disabled'])
        self.status.config(text='Calculando…')

        def job():
            F = S.parse_rhs(v["y' ="])
            x0s, y0s = S.parse_value(v['x₀']), S.parse_value(v['y₀'])
            x0, y0 = float(x0s), float(y0s)
            h, xf = float(S.parse_value(v['h'])), float(S.parse_value(v['x final']))
            n = int(round((xf - x0) / h))
            if h <= 0 or n <= 0 or n > 20000:
                raise S.ParseError("Revisa h y x final (se necesitan entre 1 y 20000 pasos).")
            f = NM.make_f(F)
            out = {'F': F, 'x0': x0, 'xf': xf, 'h': h}
            out['xs'], out['euler'] = NM.euler(f, x0, y0, h, n)
            _, out['heun'] = NM.heun(f, x0, y0, h, n)
            _, out['rk4'] = NM.rk4(f, x0, y0, h, n)
            out['exact'] = None
            try:
                ode = sp.Eq(S.yf(S.x).diff(S.x), F.subs(S.Y, S.yf(S.x)))
                ics = {S.yf(sp.nsimplify(x0s, rational=True)): sp.nsimplify(y0s, rational=True)}
                sol = sp.dsolve(ode, S.yf(S.x), ics=ics)
                sol = sol[0] if isinstance(sol, list) else sol
                if not sol.rhs.has(sp.Integral):
                    out['exact'] = sol.rhs
            except Exception:
                pass
            return out

        self.app.run(job, self.show, self.fail)

    def fail(self, e):
        self.btn.state(['!disabled'])
        self.status.config(text=f"Error: {e}", foreground='#c53030')

    def show(self, o):
        self.btn.state(['!disabled'])
        xs = o['xs']
        ex = NM.eval_expr(o['exact'], xs) if o['exact'] is not None else None
        self.tree.delete(*self.tree.get_children())
        fmt = lambda a: f"{a:.6f}" if np.isfinite(a) else '—'
        for i, xv in enumerate(xs):
            row = [i, f"{xv:.4f}", fmt(o['euler'][i]), fmt(o['heun'][i]), fmt(o['rk4'][i])]
            if ex is not None:
                row += [fmt(ex[i]), f"{abs(o['euler'][i] - ex[i]):.2e}", f"{abs(o['rk4'][i] - ex[i]):.2e}"]
            else:
                row += ['—', '—', '—']
            self.tree.insert('', 'end', values=row)
        msg = f"y' = {o['F']}   ·   {len(xs) - 1} pasos con h = {o['h']}"
        if ex is not None:
            msg += f"   ·   Exacta: y = {o['exact']}"
        else:
            msg += "   ·   (sin solución exacta disponible)"
        self.status.config(text=msg, foreground=ACCENT)

        ax = self.plot.ax
        ax.clear()
        if ex is not None:
            xd = np.linspace(xs[0], xs[-1], 400)
            ax.plot(xd, NM.eval_expr(o['exact'], xd), color='k', lw=2, label='Exacta')
        mk = dict(ms=4) if len(xs) <= 60 else dict(ms=0)
        ax.plot(xs, o['euler'], 'o--', color=COLORS[1], lw=1.2, label='Euler', **mk)
        ax.plot(xs, o['heun'], 's--', color=COLORS[3], lw=1.2, label='Heun', **mk)
        ax.plot(xs, o['rk4'], '^-', color=COLORS[0], lw=1.4, label='RK4', **mk)
        self.plot.finish("Comparación de métodos numéricos")


class FieldTab(ttk.Frame):
    def __init__(self, master, app):
        super().__init__(master, padding=10)
        self.app = app
        top = ttk.Frame(self)
        top.pack(fill='x')
        ttk.Label(top, text="y' =").pack(side='left')
        self.F = ttk.Entry(top, width=30, font=('Consolas', 12))
        self.F.insert(0, FIELD_DEFAULT)
        self.F.pack(side='left', padx=4)
        self.F.bind('<Return>', lambda e: self.draw())
        self.lim = {}
        for lbl, val in (('x min', '-5'), ('x max', '5'), ('y min', '-3'), ('y max', '5')):
            ttk.Label(top, text=lbl).pack(side='left', padx=(10, 2))
            e = ttk.Entry(top, width=5, font=('Consolas', 12))
            e.insert(0, val)
            e.pack(side='left')
            e.bind('<Return>', lambda ev: self.draw())
            self.lim[lbl] = e
        self.iso = tk.BooleanVar(value=True)
        ttk.Checkbutton(top, text="Isoclina y' = 0", variable=self.iso,
                        command=self.draw).pack(side='left', padx=10)
        ttk.Button(top, text='Dibujar ▶', style='Accent.TButton', command=self.draw).pack(side='left')
        ttk.Button(top, text='Borrar curvas', command=self.draw).pack(side='left', padx=6)
        self.status = ttk.Label(self, text='Haz clic en la gráfica para trazar la solución que pasa por ese punto.',
                                foreground=ACCENT)
        self.status.pack(fill='x', pady=(6, 0))
        self.plot = PlotPanel(self)
        self.plot.pack(fill='both', expand=True, pady=(6, 0))
        self.plot.canvas.mpl_connect('button_press_event', self.on_click)
        self.f = None
        self.n_curves = 0
        self.after(200, self.draw)

    def draw(self):
        try:
            F = S.parse_rhs(self.F.get())
            xr = (float(self.lim['x min'].get()), float(self.lim['x max'].get()))
            yr = (float(self.lim['y min'].get()), float(self.lim['y max'].get()))
        except Exception as e:
            self.status.config(text=f"Error: {e}", foreground='#c53030')
            return
        self.f, self.xr, self.yr, self.n_curves = NM.make_f(F), xr, yr, 0
        ax = self.plot.ax
        ax.clear()
        NM.slope_field(ax, self.f, xr, yr, n=26, isocline=self.iso.get())
        ax.set_xlim(xr)
        ax.set_ylim(yr)
        self.status.config(text="Haz clic en la gráfica para trazar la solución que pasa por ese punto.",
                           foreground=ACCENT)
        self.plot.finish(f"Campo de direcciones de y' = {F}")

    def on_click(self, ev):
        if self.f is None or ev.inaxes is not self.plot.ax or self.plot.toolbar.mode:
            return
        color = COLORS[self.n_curves % len(COLORS)]
        self.n_curves += 1
        span = self.yr[1] - self.yr[0]
        for xe in self.xr:
            xs, ys = NM.curve(self.f, ev.xdata, ev.ydata, xe, max(abs(self.yr[0]), abs(self.yr[1])) + 10 * span)
            self.plot.ax.plot(xs, ys, color=color, lw=2)
        self.plot.ax.plot([ev.xdata], [ev.ydata], 'o', color=color)
        self.plot.ax.set_xlim(self.xr)
        self.plot.ax.set_ylim(self.yr)
        self.status.config(text=f"Curva por ({ev.xdata:.2f}, {ev.ydata:.2f})", foreground=ACCENT)
        self.plot.canvas.draw_idle()


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('DiffEq Studio — Ecuaciones Diferenciales')
        self.geometry('1360x820')
        self.minsize(1000, 640)
        self.configure(bg=BG)
        self._q = queue.Queue()

        style = ttk.Style(self)
        style.theme_use('clam')
        style.configure('.', background=BG, font=('Segoe UI', 10))
        style.configure('TNotebook.Tab', padding=(16, 7), font=('Segoe UI', 10, 'bold'))
        style.map('TNotebook.Tab', background=[('selected', 'white')], foreground=[('selected', ACCENT)])
        style.configure('Accent.TButton', background=ACCENT, foreground='white',
                        font=('Segoe UI', 10, 'bold'), padding=(14, 6))
        style.map('Accent.TButton', background=[('active', '#2c5282'), ('disabled', '#a0aec0')])
        style.configure('Treeview', font=('Consolas', 10), rowheight=22)
        style.configure('Treeview.Heading', font=('Segoe UI', 9, 'bold'))

        header = tk.Frame(self, bg=ACCENT)
        header.pack(fill='x')
        tk.Label(header, text='DiffEq Studio', bg=ACCENT, fg='white',
                 font=('Segoe UI', 16, 'bold')).pack(side='left', padx=14, pady=8)
        tk.Label(header, text='Cálculo 4 · Ecuaciones diferenciales paso a paso', bg=ACCENT,
                 fg='#bee3f8', font=('Segoe UI', 10)).pack(side='left')

        nb = ttk.Notebook(self)
        nb.pack(fill='both', expand=True, padx=8, pady=8)
        nb.add(SolverTab(nb, self, FIRST_ORDER, FIRST_ORDER[0]), text='Primer orden')
        nb.add(SolverTab(nb, self, HIGHER_ORDER, HIGHER_ORDER[0]), text='Orden superior')
        nb.add(LaplaceTab(nb, self), text='Laplace')
        nb.add(SystemsTab(nb, self), text='Sistemas')
        nb.add(NumericTab(nb, self), text='Métodos numéricos')
        nb.add(FieldTab(nb, self), text='Campo de direcciones')
        help_frame = ttk.Frame(nb, padding=16)
        txt = tk.Text(help_frame, font=('Consolas', 11), relief='flat', bg='white', padx=14, pady=10)
        txt.insert('1.0', HELP)
        txt.config(state='disabled')
        txt.pack(fill='both', expand=True)
        nb.add(help_frame, text='Ayuda')
        self.after(80, self._poll)

    def run(self, func, on_done, on_error):
        """Ejecuta func en segundo plano y entrega el resultado en el hilo de la interfaz."""
        def worker():
            try:
                r = func()
                self._q.put((on_done, r))
            except Exception as e:
                self._q.put((on_error, e))
        threading.Thread(target=worker, daemon=True).start()

    def _poll(self):
        try:
            while True:
                cb, arg = self._q.get_nowait()
                cb(arg)
        except queue.Empty:
            pass
        self.after(80, self._poll)


if __name__ == '__main__':
    App().mainloop()
