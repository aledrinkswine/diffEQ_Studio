"""Métodos numéricos y utilidades de evaluación para las gráficas."""
import numpy as np
import sympy as sp
from matplotlib.collections import LineCollection

from solvers import Y, x


def make_f(F):
    """Convierte F(x, y) simbólica en función numérica vectorizada."""
    f = sp.lambdify((x, Y), F, 'numpy')

    def g(xv, yv):
        with np.errstate(all='ignore'):
            return np.asarray(f(xv, yv), dtype=float) + 0 * np.asarray(xv, dtype=float)
    return g


def eval_expr(expr, xs, var=x):
    """Evalúa y(x) en un arreglo; devuelve NaN donde no es real o no está definida."""
    try:
        f = sp.lambdify(var, expr, 'numpy')
        with np.errstate(all='ignore'):
            v = np.asarray(f(xs), dtype=complex) * np.ones_like(xs)
    except Exception:
        def one(xv):
            try:
                return complex(expr.subs(var, xv).evalf())
            except Exception:
                return np.nan
        v = np.array([one(xv) for xv in xs], dtype=complex)
    return np.where(np.abs(v.imag) < 1e-9, v.real, np.nan)


def euler(f, x0, y0, h, n):
    xs = x0 + h * np.arange(n + 1)
    ys = np.empty(n + 1)
    ys[0] = y0
    for i in range(n):
        ys[i + 1] = ys[i] + h * f(xs[i], ys[i])
    return xs, ys


def heun(f, x0, y0, h, n):
    xs = x0 + h * np.arange(n + 1)
    ys = np.empty(n + 1)
    ys[0] = y0
    for i in range(n):
        k1 = f(xs[i], ys[i])
        k2 = f(xs[i] + h, ys[i] + h * k1)
        ys[i + 1] = ys[i] + h * (k1 + k2) / 2
    return xs, ys


def rk4(f, x0, y0, h, n):
    xs = x0 + h * np.arange(n + 1)
    ys = np.empty(n + 1)
    ys[0] = y0
    for i in range(n):
        xi, yi = xs[i], ys[i]
        k1 = f(xi, yi)
        k2 = f(xi + h / 2, yi + h * k1 / 2)
        k3 = f(xi + h / 2, yi + h * k2 / 2)
        k4 = f(xi + h, yi + h * k3)
        ys[i + 1] = yi + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return xs, ys


def curve(f, x0, y0, x_end, ybound, steps=500):
    """Trayectoria RK4 desde (x0, y0) hasta x_end; se detiene si diverge."""
    xs, ys = rk4(f, x0, y0, (x_end - x0) / steps, steps)
    bad = ~np.isfinite(ys) | (np.abs(ys) > ybound)
    if bad.any():
        k = np.argmax(bad)
        xs, ys = xs[:k + 1], ys[:k + 1]
    return xs, ys


def make_system(rhs, vars_, t):
    """Convierte un sistema simbólico en f(t, X) numérica."""
    f = sp.lambdify((t, *vars_), rhs, 'numpy')

    def g(tv, X):
        with np.errstate(all='ignore'):
            return np.array([np.asarray(c, dtype=float) + 0 * np.asarray(X[0], dtype=float)
                             for c in f(tv, *X)])
    return g


def rk4_system(f, t0, X0, T, steps=2000, bound=1e4):
    """Integra X' = f(t, X) desde t0 durante T (T < 0 integra hacia atrás)."""
    h = T / steps
    ts = t0 + h * np.arange(steps + 1)
    Xs = np.empty((steps + 1, len(X0)))
    Xs[0] = X0
    for i in range(steps):
        ti, Xi = ts[i], Xs[i]
        k1 = f(ti, Xi)
        k2 = f(ti + h / 2, Xi + h * k1 / 2)
        k3 = f(ti + h / 2, Xi + h * k2 / 2)
        k4 = f(ti + h, Xi + h * k3)
        Xs[i + 1] = Xi + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if not np.all(np.isfinite(Xs[i + 1])) or np.max(np.abs(Xs[i + 1])) > bound:
            return ts[:i + 1], Xs[:i + 1]
    return ts, Xs


def slope_field(ax, f, xr, yr, n=22, isocline=False):
    X, Yg = np.meshgrid(np.linspace(*xr, n), np.linspace(*yr, n))
    S = f(X, Yg)
    # dirección en píxeles: así todos los segmentos miden lo mismo en pantalla
    sx, sy = xr[1] - xr[0], yr[1] - yr[0]
    box = ax.get_window_extent()
    W, H = max(box.width, 1), max(box.height, 1)
    vertical = np.isinf(S)
    U = np.where(vertical, 0.0, W / sx)
    V = np.where(vertical, 1.0, np.nan_to_num(S, posinf=0, neginf=0) * H / sy)
    L = np.hypot(U, V)
    half = 0.42 * min(W, H) / n
    dx, dy = half * U / L * sx / W, half * V / L * sy / H
    segs = np.stack([np.stack([X - dx, Yg - dy], -1), np.stack([X + dx, Yg + dy], -1)], -2)
    segs = segs.reshape(-1, 2, 2)[np.isfinite(segs.reshape(-1, 4)).all(axis=1)]
    ax.add_collection(LineCollection(segs, colors='#8796ab', linewidths=1.4, capstyle='round'))
    if isocline:
        Xf, Yf = np.meshgrid(np.linspace(*xr, 200), np.linspace(*yr, 200))
        Sf = f(Xf, Yf)
        if np.nanmin(Sf) < 0 < np.nanmax(Sf):
            ax.contour(Xf, Yf, Sf, [0], colors='#d69e2e', linewidths=1.2, linestyles='--')
