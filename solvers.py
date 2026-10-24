"""Motor matemático de DiffEq Studio: parseo de EDOs y resolución paso a paso."""
import re

import sympy as sp
from sympy.parsing.sympy_parser import (
    convert_xor, function_exponentiation, implicit_application,
    implicit_multiplication, parse_expr, standard_transformations)

x = sp.Symbol('x')
Y = sp.Symbol('y')
yf = sp.Function('y')
YP = [sp.Symbol(f'YP{k}') for k in range(7)]
DX, DY = sp.symbols('dx dy')
C1 = sp.Symbol('C1')

_T = standard_transformations + (convert_xor, implicit_multiplication,
                                 implicit_application, function_exponentiation)
_LOCAL = {
    'x': x, 'y': YP[0], 'e': sp.E, 'E': sp.E, 'pi': sp.pi,
    'ln': sp.log, 'log': sp.log, 'sen': sp.sin, 'tg': sp.tan,
    'arcsen': sp.asin, 'arcsin': sp.asin, 'arccos': sp.acos,
    'arctan': sp.atan, 'arctg': sp.atan, 'raiz': sp.sqrt,
    'dx': DX, 'dy': DY,
}
for _k in range(1, 7):
    _LOCAL[f'YP{_k}'] = YP[_k]

HINT_NAMES = {
    'separable': 'Variables separables',
    '1st_linear': 'Lineal de primer orden',
    'Bernoulli': 'Bernoulli',
    '1st_exact': 'Ecuación exacta',
    '1st_homogeneous_coeff_best': 'Homogénea',
    '1st_homogeneous_coeff_subs_dep_div_indep': 'Homogénea (y = vx)',
    'Riccati_special_minus2': 'Riccati',
    'nth_linear_euler_eq_homogeneous': 'Cauchy-Euler (homogénea)',
    'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters': 'Cauchy-Euler (variación de parámetros)',
    'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients': 'Cauchy-Euler (coef. indeterminados)',
    '2nd_power_series_ordinary': 'Series de potencias (punto ordinario)',
    '2nd_power_series_regular': 'Frobenius (punto singular regular)',
    'nth_order_reducible': 'Reducción de orden',
    'Liouville': 'Liouville',
    'factorable': 'Factorizable',
    'almost_linear': 'Casi lineal',
    'lie_group': 'Grupos de Lie',
}


class ParseError(Exception):
    pass


def _p(s, extra=None):
    try:
        return parse_expr(s, local_dict={**_LOCAL, **(extra or {})}, transformations=_T)
    except Exception as e:
        raise ParseError(f"No pude interpretar: '{s.strip()}'") from e


def _normalize(text):
    s = text.strip().replace('′', "'").replace('’', "'").replace('−', '-')
    for pat, rep in (('d^2y/dx^2', 'YP2'), ('d2y/dx2', 'YP2'), ('dy/dx', 'YP1')):
        s = s.replace(pat, rep)
    for k in range(6, 0, -1):
        s = s.replace('y' + "'" * k, f'YP{k}')
    return s


def parse_value(s, extra=None):
    return _p(_normalize(str(s)), extra)


def parse_equation(text, extra=None):
    """Devuelve ('ode', (expr, orden)) o ('exact', (M, N)) para M dx + N dy = 0."""
    s = _normalize(text)
    if not s:
        raise ParseError("Escribe una ecuación.")
    if s.count('=') > 1:
        raise ParseError("La ecuación debe tener un solo '='.")
    lhs, rhs = (s.split('=') + ['0'])[:2]
    expr = _p(lhs, extra) - _p(rhs if rhs.strip() else '0', extra)
    if expr.has(DX) or expr.has(DY):
        M = sp.simplify(sp.diff(expr, DX).subs(YP[0], Y))
        N = sp.simplify(sp.diff(expr, DY).subs(YP[0], Y))
        if N == 0:
            raise ParseError("Falta el término con dy.")
        return 'exact', (M, N)
    order = max((k for k in range(7) if expr.has(YP[k]) and k > 0), default=0)
    if order < 1:
        raise ParseError("No encontré derivadas de y (usa y', y'' o dy/dx).")
    return 'ode', (expr, order)


def to_ode(expr):
    subs = {YP[k]: yf(x).diff(x, k) for k in range(1, 7)}
    subs[YP[0]] = yf(x)
    return sp.Eq(expr.subs(subs, simultaneous=True), 0)


def parse_rhs(text):
    """Para métodos numéricos/campo de direcciones: devuelve F(x, y) de y' = F(x, y)."""
    s = _normalize(text)
    if '=' not in s and 'YP1' not in s and 'dx' not in s:
        return _p(s).subs(YP[0], Y)
    kind, data = parse_equation(text)
    if kind == 'exact':
        M, N = data
        return sp.simplify(-M / N)
    expr, n = data
    if n != 1:
        raise ParseError("Debe ser una ecuación de primer orden.")
    sols = sp.solve(expr.subs(YP[0], Y), YP[1])
    if len(sols) != 1:
        raise ParseError("No pude despejar y'.")
    return sols[0]


_IC = re.compile(r"y\s*('*)\s*\(\s*([^)]+?)\s*\)\s*=\s*([^,;]+)")


def parse_ics(text):
    s = (text or '').replace('′', "'").replace('’', "'")
    out = [(len(m.group(1)), sp.nsimplify(parse_value(m.group(2)), rational=True),
            sp.nsimplify(parse_value(m.group(3).strip()), rational=True))
           for m in _IC.finditer(s)]
    if s.strip() and not out:
        raise ParseError("Condiciones iniciales: usa el formato y(0)=1, y'(0)=2")
    return out


def dsym(k):
    return Y if k == 0 else sp.Symbol('y' + "'" * k)


def yeq(e):
    return sp.Eq(Y, e)


class Result:
    def __init__(self):
        self.steps = []
        self.method = ''
        self.order = 1
        self.F = None              # y' = F(x, y) para primer orden
        self.general = []          # soluciones explícitas y = ...
        self.implicit = None       # Eq(G(x, y), C1)
        self.particular = None     # solución explícita del PVI
        self.particular_implicit = None
        self.ics = []

    def title(self, t): self.steps.append(('title', t))
    def text(self, t): self.steps.append(('text', t))
    def math(self, e): self.steps.append(('math', e))
    def ok(self, t): self.steps.append(('ok', t))
    def warn(self, t): self.steps.append(('warn', t))


# ---------------------------------------------------------------- primer orden

def _solve_implicit(G, res):
    """Intenta despejar y de G(x, y) = C1."""
    try:
        sols = [s for s in sp.solve(sp.Eq(G, C1), Y) if not s.has(Y)]
    except Exception:
        sols = []
    if not sols:
        res.text("No es posible despejar y: la solución general queda en forma implícita.")
        res.implicit = sp.Eq(G, C1)
        return
    renamed, out = False, []
    for s in sols:
        s2 = sp.expand_power_exp(s)
        if s2.has(sp.exp(C1)):
            s2 = s2.subs(sp.exp(C1), C1)
            renamed = True
        out.append(s2)
    res.text("Despejando y:")
    for s in out:
        res.math(yeq(s))
    if renamed:
        res.text("(Se renombró la constante: e^C₁ → C₁)")
    res.general = out


def _separable(F, res):
    if not F.has(Y):
        res.method = "Integración directa"
        res.title("Método: Integración directa")
        res.text("y' depende solo de x, así que integramos directamente:")
        res.math(yeq(sp.Integral(F, x)))
        sol = sp.integrate(F, x) + C1
        res.math(yeq(sol))
        res.general = [sol]
        return True
    d = sp.separatevars(F, symbols=[x, Y], dict=True)
    if not d:
        return False
    g = sp.simplify(d['coeff'] * d[x])
    h = sp.simplify(d[Y])
    res.method = "Variables separables"
    res.title("Método: Variables separables")
    res.text("La ecuación tiene la forma y' = g(x)·h(y) con:")
    res.math(sp.Eq(sp.Symbol('g(x)'), g))
    res.math(sp.Eq(sp.Symbol('h(y)'), h))
    res.text("Separamos variables e integramos ambos lados:")
    res.math(sp.Eq(sp.Integral(1 / h, Y), sp.Integral(g, x)))
    L = sp.integrate(1 / h, Y)
    R = sp.integrate(g, x)
    res.math(sp.Eq(L, R + C1))
    _solve_implicit(L - R, res)
    eqs = [e for e in sp.solve(h, Y) if not e.has(x)]
    if eqs:
        res.text("Soluciones de equilibrio (h(y) = 0), que se pierden al dividir:")
        for e in eqs:
            res.math(yeq(e))
    return True


def _linear_steps(P, Q, dep, res):
    """Resuelve dep' + P dep = Q con factor integrante. Devuelve dep(x)."""
    res.math(sp.Eq(sp.Symbol(dep.name + "'") + P * dep, Q))
    res.text("Factor integrante μ(x) = e^(∫P dx):")
    mu = sp.simplify(sp.exp(sp.integrate(P, x)))
    res.math(sp.Eq(sp.Symbol('μ(x)'), mu))
    res.text("Multiplicando por μ, el lado izquierdo es la derivada de un producto:")
    muQ = sp.simplify(mu * Q)
    res.math(sp.Eq(sp.Derivative(mu * dep, x), muQ))
    res.text("Integramos ambos lados:")
    I = sp.integrate(muQ, x)
    res.math(sp.Eq(mu * dep, I + C1))
    sol = sp.simplify((I + C1) / mu)
    res.text(f"Despejamos {dep.name}:")
    res.math(sp.Eq(dep, sol))
    return sol


def _linear(F, res):
    a = sp.diff(F, Y)
    if a.has(Y):
        return False
    Q = sp.simplify(F - a * Y)
    if Q.has(Y):
        return False
    P = sp.simplify(-a)
    res.method = "Lineal de primer orden (factor integrante)"
    res.title("Método: Ecuación lineal — factor integrante")
    res.text("Forma estándar y' + P(x)·y = Q(x) con:")
    res.math(sp.Eq(sp.Symbol('P(x)'), P))
    res.math(sp.Eq(sp.Symbol('Q(x)'), Q))
    res.general = [_linear_steps(P, Q, Y, res)]
    return True


def _bernoulli(F, res):
    coef = {}
    for t in sp.Add.make_args(sp.expand(F)):
        c, yp = t.as_independent(Y, as_Add=False)
        if yp == 1:
            n = sp.Integer(0)
        elif yp == Y:
            n = sp.Integer(1)
        elif yp.is_Pow and yp.base == Y and not yp.exp.has(x, Y):
            n = yp.exp
        else:
            return False
        coef[n] = coef.get(n, 0) + c
    others = set(coef) - {1}
    if 1 not in coef or len(others) != 1:
        return False
    n = others.pop()
    if n == 0:
        return False
    P, Q = sp.simplify(-coef[1]), sp.simplify(coef[n])
    res.method = "Bernoulli"
    res.title("Método: Ecuación de Bernoulli")
    res.text(f"Forma y' + P(x)·y = Q(x)·yⁿ con n = {n}:")
    res.math(sp.Eq(sp.Symbol('P(x)'), P))
    res.math(sp.Eq(sp.Symbol('Q(x)'), Q))
    v = sp.Symbol('v')
    res.text("Sustitución v = y^(1−n):")
    res.math(sp.Eq(v, Y ** (1 - n)))
    res.text("Así v' = (1−n)·y^(−n)·y', y la ecuación se vuelve lineal en v:")
    sol_v = _linear_steps(sp.simplify((1 - n) * P), sp.simplify((1 - n) * Q), v, res)
    p = 1 / (1 - n)
    sol = sol_v ** p
    res.text("Regresamos a y = v^(1/(1−n)):")
    res.math(yeq(sol))
    res.general = [sol]
    if p.is_Rational and p.q % 2 == 0:
        res.text("(También la rama con signo negativo:)")
        res.math(yeq(-sol))
        res.general.append(-sol)
    if n > 0:
        res.text("Además y = 0 es solución.")
    return True


def _homogeneous(F, res):
    t = sp.Symbol('t', positive=True)
    if sp.simplify(F.subs({x: t * x, Y: t * Y}, simultaneous=True) - F) != 0:
        return False
    v = sp.Symbol('v')
    G = sp.simplify(F.subs(Y, v * x))
    if G.has(x):
        return False
    res.method = "Homogénea (y = v·x)"
    res.title("Método: Ecuación homogénea")
    res.text("F(tx, ty) = F(x, y): es homogénea. Sustituimos y = v·x, y' = v + x·v':")
    res.math(sp.Eq(v + x * sp.Symbol("v'"), G))
    D = sp.simplify(G - v)
    if D == 0:
        res.text("v' = 0 ⇒ v = C₁:")
        res.math(yeq(C1 * x))
        res.general = [C1 * x]
        return True
    res.text("Separando variables:")
    res.math(sp.Eq(sp.Integral(1 / D, v), sp.Integral(1 / x, x)))
    L = sp.integrate(1 / D, v)
    res.math(sp.Eq(L, sp.log(x) + C1))
    res.text("Regresando a v = y/x:")
    res.math(sp.Eq(L.subs(v, Y / x), sp.log(x) + C1))
    _solve_implicit(L.subs(v, Y / x) - sp.log(x), res)
    return True


def _first_order(expr, res):
    sols = sp.solve(expr.subs(YP[0], Y), YP[1])
    if len(sols) != 1:
        return False
    F = sp.simplify(sols[0])
    res.F = F
    res.title("Forma normal")
    res.math(sp.Eq(dsym(1), F))
    return any(fn(F, res) for fn in (_separable, _linear, _bernoulli, _homogeneous))


def _exact(M, N, res):
    res.F = sp.simplify(-M / N)
    res.title("Forma diferencial M dx + N dy = 0")
    res.math(sp.Eq(sp.Symbol('M'), M))
    res.math(sp.Eq(sp.Symbol('N'), N))
    My, Nx = sp.diff(M, Y), sp.diff(N, x)
    res.text("Criterio de exactitud ∂M/∂y = ∂N/∂x:")
    res.math(sp.Eq(sp.Symbol('∂M/∂y'), My))
    res.math(sp.Eq(sp.Symbol('∂N/∂x'), Nx))
    if sp.simplify(My - Nx) != 0:
        res.warn("No es exacta. Buscamos un factor integrante.")
        r = sp.simplify((My - Nx) / N)
        if not r.has(Y):
            res.text("(M_y − N_x)/N depende solo de x ⇒ μ(x) = e^(∫ (M_y − N_x)/N dx):")
            res.math(sp.Eq(sp.Symbol('(M_y − N_x)/N'), r))
            mu = sp.simplify(sp.exp(sp.integrate(r, x)))
        else:
            r = sp.simplify((Nx - My) / M)
            if r.has(x):
                res.warn("No se encontró un factor integrante μ(x) ni μ(y).")
                return False
            res.text("(N_x − M_y)/M depende solo de y ⇒ μ(y) = e^(∫ (N_x − M_y)/M dy):")
            res.math(sp.Eq(sp.Symbol('(N_x − M_y)/M'), r))
            mu = sp.simplify(sp.exp(sp.integrate(r, Y)))
        res.math(sp.Eq(sp.Symbol('μ'), mu))
        M, N = sp.simplify(mu * M), sp.simplify(mu * N)
        res.text("Multiplicando por μ, la nueva ecuación es exacta:")
        res.math(sp.Eq(sp.Symbol('M'), M))
        res.math(sp.Eq(sp.Symbol('N'), N))
        res.method = "Exacta con factor integrante"
    else:
        res.ok("✔ Es exacta.")
        res.method = "Ecuación exacta"
    res.title("Cálculo de la función potencial ψ(x, y)")
    res.text("Buscamos ψ con ψ_x = M y ψ_y = N. Integrando M respecto de x:")
    psi = sp.integrate(M, x)
    res.math(sp.Eq(sp.Symbol('ψ'), psi + sp.Symbol('g(y)')))
    res.text("Derivando respecto de y e igualando a N:")
    gp = sp.simplify(N - sp.diff(psi, Y))
    res.math(sp.Eq(sp.Symbol("g'(y)"), gp))
    g = sp.integrate(gp, Y)
    res.math(sp.Eq(sp.Symbol('g(y)'), g))
    G = psi + g
    res.text("Solución general ψ(x, y) = C₁:")
    res.math(sp.Eq(G, C1))
    _solve_implicit(G, res)
    return True


# ---------------------------------------------------- orden superior (coef. ctes)

def _linear_cc(expr, n, res):
    a = [sp.expand(sp.diff(expr, YP[k])) for k in range(n + 1)]
    bad = {x, *YP}
    if any(c.free_symbols & bad for c in a):
        return False
    rest = sp.simplify(expr - sum(a[k] * YP[k] for k in range(n + 1)))
    if rest.free_symbols & set(YP):
        return False
    g = sp.simplify(-rest)
    res.method = "Lineal con coeficientes constantes" + (
        "" if g == 0 else " + variación de parámetros")
    res.title("Ecuación lineal con coeficientes constantes")
    res.math(sp.Eq(sum(a[k] * dsym(k) for k in range(n + 1)), g))

    res.title("1) Solución homogénea")
    r = sp.Symbol('r')
    char = sum(a[k] * r ** k for k in range(n + 1))
    res.text("Proponemos y = e^(rx) y obtenemos la ecuación característica:")
    res.math(sp.Eq(char, 0))
    if n == 2:
        disc = sp.simplify(a[1] ** 2 - 4 * a[2] * a[0])
        res.math(sp.Eq(sp.Symbol('Δ = b² − 4ac'), disc))
        if disc.is_number:
            res.text("Δ > 0: dos raíces reales distintas." if disc > 0 else
                     "Δ = 0: raíz real repetida." if disc == 0 else
                     "Δ < 0: raíces complejas conjugadas.")
    rts = sp.roots(sp.Poly(char, r))
    if sum(rts.values()) != n:
        res.warn("No se pudieron hallar todas las raíces de forma exacta.")
        return False
    basis = []
    for rt, m in rts.items():
        re_, im_ = (sp.simplify(p) for p in rt.as_real_imag())
        mult = f" (multiplicidad {m})" if m > 1 else ""
        if im_ == 0:
            res.text(f"Raíz real r = {rt}{mult}")
            basis += [x ** j * sp.exp(rt * x) for j in range(m)]
        elif im_.is_positive:
            res.text(f"Raíces complejas r = α ± βi con α = {re_}, β = {im_}{mult}")
            for j in range(m):
                basis += [x ** j * sp.exp(re_ * x) * sp.cos(im_ * x),
                          x ** j * sp.exp(re_ * x) * sp.sin(im_ * x)]
    res.text("Conjunto fundamental de soluciones:")
    res.math(sp.Tuple(*basis))
    Cs = sp.symbols(f'C1:{n + 1}')
    yh = sum(c * b for c, b in zip(Cs, basis))
    res.math(sp.Eq(sp.Symbol('y_h'), yh))
    if g == 0:
        res.general = [yh]
        return True

    res.title("2) Solución particular — variación de parámetros")
    W = sp.Matrix([[sp.diff(b, x, i) for b in basis] for i in range(n)])
    res.text("Wronskiano:")
    res.math(sp.Eq(sp.Symbol('W'), sp.simplify(W.det())))
    res.text(f"Resolvemos el sistema W·u' = (0, …, g(x)/{a[n]}) y luego integramos:")
    up = W.LUsolve(sp.Matrix([0] * (n - 1) + [g / a[n]]))
    yp = 0
    for i, b in enumerate(basis):
        ui_p = sp.simplify(up[i])
        ui = sp.integrate(ui_p, x).replace(sp.log, lambda a: sp.log(sp.Abs(a)))
        res.math(sp.Eq(sp.Symbol(f"u{i + 1}'"), ui_p))
        res.math(sp.Eq(sp.Symbol(f"u{i + 1}"), ui))
        yp += ui * b
    yp = sp.simplify(yp)
    terms = sp.Add.make_args(sp.expand(yp))
    kept = [t for t in terms if not any(not sp.simplify(t / b).has(x) for b in basis)]
    if len(kept) < len(terms):
        yp = sp.simplify(sp.Add(*kept))
        res.text("(Se omiten los términos que ya forman parte de y_h.)")
    res.math(sp.Eq(sp.Symbol('y_p'), yp))
    res.title("3) Solución general y = y_h + y_p")
    sol = yh + yp
    res.math(yeq(sol))
    res.general = [sol]
    return True


# --------------------------------------------------------------------- general

def _fallback(ode, res):
    hints = sp.classify_ode(ode, yf(x))
    if not hints:
        raise ParseError("No conozco un método para esta ecuación.")
    named = [h for h in hints if h in HINT_NAMES and h != 'factorable']
    hint = named[0] if named else hints[0]
    res.method = HINT_NAMES.get(hint, hint)
    res.title("Resolución con SymPy")
    res.text(f"Método detectado: {res.method}")
    res.warn("(Este tipo aún no tiene explicación paso a paso; se muestra el resultado.)")
    sols = sp.dsolve(ode, yf(x))
    for s in (sols if isinstance(sols, list) else [sols]):
        res.math(s)
        if s.lhs == yf(x) and not s.rhs.has(yf(x)):
            res.general.append(s.rhs)


def _apply_ics(res, ode, ics):
    res.title("Problema de valor inicial")
    for k, x0, v in ics:
        res.math(sp.Eq(sp.Symbol(f"y{chr(39) * k}({x0})"), v))
    if res.implicit is not None and not res.general:
        _, x0, v = ics[0]
        G = res.implicit.lhs
        cval = sp.simplify(G.subs({x: x0, Y: v}))
        res.text("Sustituyendo en la solución implícita:")
        res.math(sp.Eq(C1, cval))
        res.particular_implicit = sp.Eq(G, cval)
        res.math(res.particular_implicit)
        return
    shown = False
    for s in res.general:
        consts = sorted((c for c in s.free_symbols if c.name.startswith('C')), key=str)
        eqs = [sp.Eq(sp.diff(s, x, k).subs(x, x0), v) for k, x0, v in ics]
        if not shown:
            res.text("Sistema para las constantes:")
            for e in eqs:
                res.math(e)
            shown = True
        try:
            sol = sp.solve(eqs, consts, dict=True)
        except Exception:
            sol = []
        if sol:
            for c, val in sol[0].items():
                res.math(sp.Eq(c, val))
            res.particular = sp.simplify(s.subs(sol[0]))
            res.text("Solución particular:")
            res.math(yeq(res.particular))
            return
    try:
        conds = {(yf(x).diff(x, k).subs(x, x0) if k else yf(x0)): v for k, x0, v in ics}
        s = sp.dsolve(ode, yf(x), ics=conds)
        s = s[0] if isinstance(s, list) else s
        res.particular = s.rhs
        res.text("Solución particular (SymPy):")
        res.math(yeq(res.particular))
    except Exception:
        res.warn("No se pudieron determinar las constantes.")


def solve(text, ics_text=''):
    res = Result()
    kind, data = parse_equation(text)
    ics = parse_ics(ics_text)
    res.ics = ics
    if kind == 'exact':
        M, N = data
        ode = sp.Eq(yf(x).diff(x), (-M / N).subs(Y, yf(x)))
        if not _exact(M, N, res):
            _fallback(ode, res)
    else:
        expr, n = data
        res.order = n
        ode = to_ode(expr)
        done = _first_order(expr, res) if n == 1 else _linear_cc(expr, n, res)
        if not done:
            _fallback(ode, res)
    if res.general:
        res.title("Verificación")
        for s in res.general:
            try:
                # ln|u| y ln(u) tienen la misma derivada; se verifica sin el valor absoluto
                ok = sp.checkodesol(ode, sp.Eq(yf(x), s.replace(sp.Abs, lambda a: a)))[0]
            except Exception:
                ok = False
            if ok:
                res.ok("✔ Al sustituir en la ecuación se cumple la igualdad.")
            else:
                res.warn("No se pudo verificar simbólicamente.")
    if ics:
        _apply_ics(res, ode, ics)
    return res


def pretty(e):
    return sp.pretty(e, use_unicode=True, num_columns=110)
