"""Transformada de Laplace paso a paso: directa, inversa, PVI y convolución."""
import sympy as sp

import solvers as S

t = sp.Symbol('t', positive=True)
s = sp.Symbol('s')
tau = sp.Symbol('τ', positive=True)
Ys = sp.Symbol('Y')

_EXTRA = {
    't': t, 'x': t, 's': s,
    'u': sp.Heaviside, 'H': sp.Heaviside, 'heaviside': sp.Heaviside, 'escalon': sp.Heaviside,
    'delta': sp.DiracDelta, 'dirac': sp.DiracDelta, 'senh': sp.sinh,
}


def parse_t(text):
    return S._p(S._normalize(text), _EXTRA).subs(S.x, t)


def parse_s(text):
    return S._p(S._normalize(text), _EXTRA)


def L(f):
    return sp.simplify(sp.laplace_transform(f, t, s, noconds=True))


def Linv(F):
    return sp.simplify(sp.inverse_laplace_transform(F, s, t))


def _lsym(f):
    return sp.Symbol(f"ℒ{{{f}}}")


# ------------------------------------------------------------------ directa

def _rule(g):
    """Describe qué fórmula de la tabla se aplica a g(t) (sin coeficiente)."""
    facs = sp.Mul.make_args(g)
    if any(f.func == sp.DiracDelta for f in facs):
        return "Impulso: ℒ{δ(t−a)} = e^(−as)"
    if any(f.func == sp.Heaviside for f in facs):
        return "Segundo teorema de traslación: ℒ{u(t−a)·g(t)} = e^(−as)·ℒ{g(t+a)}"
    exps = [f for f in facs if f.func == sp.exp and f.has(t)]
    others = [f for f in facs if f.has(t) and f not in exps]
    if exps and others:
        return "Primer teorema de traslación: ℒ{e^(at)·f(t)} = F(s−a)"
    tpows = [f for f in others if f == t or (f.is_Pow and f.base == t)]
    if tpows and len(others) > len(tpows):
        return "Derivada de la transformada: ℒ{tⁿ·f(t)} = (−1)ⁿ·dⁿF/dsⁿ"
    if not g.has(t):
        return "ℒ{1} = 1/s"
    if g == t or (g.is_Pow and g.base == t):
        return "ℒ{tⁿ} = n!/s^(n+1)"
    if g.func == sp.exp:
        return "ℒ{e^(at)} = 1/(s−a)"
    if g.func == sp.sin:
        return "ℒ{sin bt} = b/(s²+b²)"
    if g.func == sp.cos:
        return "ℒ{cos bt} = s/(s²+b²)"
    if g.func == sp.sinh:
        return "ℒ{sinh bt} = b/(s²−b²)"
    if g.func == sp.cosh:
        return "ℒ{cosh bt} = s/(s²−b²)"
    return "Definición ℒ{f} = ∫₀^∞ e^(−st) f(t) dt"


def forward_steps(f, res):
    res.title("Linealidad: se transforma término a término")
    terms = sp.Add.make_args(f)
    if len(terms) > 1:
        res.math(sp.Eq(_lsym(f), sp.Add(*[_lsym(tm) for tm in terms], evaluate=False)))
    parts = []
    for tm in terms:
        c, g = tm.as_independent(t, as_Add=False)
        res.text(f"• {_rule(g)}")
        H = [h for h in sp.Mul.make_args(g) if h.func == sp.Heaviside]
        if H:
            a = sp.solve(H[0].args[0], t)[0]
            inner = sp.simplify((g / H[0]).subs(t, t + a))
            res.math(sp.Eq(_lsym(g), sp.exp(-a * s) * _lsym(inner)))
        Fg = L(c * g)
        res.math(sp.Eq(_lsym(c * g), Fg))
        parts.append(Fg)
    F = sp.Add(*parts)
    res.title("Resultado")
    res.math(sp.Eq(sp.Symbol('F(s)'), F))
    Fs = sp.simplify(sp.together(F))
    if Fs != F:
        res.text("Con denominador común:")
        res.math(sp.Eq(sp.Symbol('F(s)'), sp.factor(Fs)))
    return F


def transform(text):
    res = S.Result()
    res.method = "Transformada de Laplace"
    f = parse_t(text)
    res.title("Función")
    res.math(sp.Eq(sp.Symbol('f(t)'), f))
    F = forward_steps(f, res)
    res.plot = [('f(t)', f)]
    res.F_s = F
    return res


# ------------------------------------------------------------------ inversa

def _inv_rule(term, res):
    num, den = sp.fraction(sp.factor(term))
    facs = [b for b, _ in sp.factor_list(den, s)[1]]
    if len(facs) == 1:
        p = sp.Poly(facs[0], s)
        if p.degree() == 1:
            res.text("• ℒ⁻¹{1/(s−a)ⁿ} = tⁿ⁻¹·e^(at)/(n−1)!")
            return
        if p.degree() == 2:
            A, B, C = p.all_coeffs()
            h = sp.simplify(B / (2 * A))
            k = sp.simplify(C / A - h ** 2)
            res.text("• Completando el cuadrado: ℒ⁻¹{(s−α)/((s−α)²+β²)} = e^(αt)cos βt,"
                     " ℒ⁻¹{β/((s−α)²+β²)} = e^(αt)sin βt")
            res.math(sp.Eq(facs[0], A * ((s + h) ** 2 + k)))
            return
    res.text("• Se aplica la tabla de transformadas inversas")


def inverse_steps(F, res):
    # cada e^(cs) se cambia por un símbolo para que expand no lo mueva al denominador
    marks = {e: (sp.Symbol(f'_E{i}'), sp.simplify(e.args[0] / s))
             for i, e in enumerate(e for e in F.atoms(sp.exp) if e.has(s))}
    G = sp.expand(F.subs({e: m for e, (m, _) in marks.items()}))
    groups = {}
    for tm in sp.Add.make_args(G):
        pw = tm.as_powers_dict()
        a = sp.simplify(-sum(pw.get(m, 0) * c for m, c in marks.values()))
        rest = tm.subs({m: 1 for m, _ in marks.values()})
        groups[a] = groups.get(a, 0) + rest
    total = 0
    for a, R in sorted(groups.items(), key=lambda kv: float(kv[0])):
        R = sp.together(R)
        if a != 0:
            res.title(f"Término con e^(−{a}s)")
            res.math(sp.Eq(sp.Symbol('G(s)'), R))
        Rp = sp.apart(R, s)
        if Rp != R:
            res.text("Fracciones parciales:")
            res.math(Rp)
        g = 0
        for tm in sp.Add.make_args(Rp):
            _inv_rule(tm, res)
            gi = Linv(tm)
            res.math(sp.Eq(sp.Symbol(f"ℒ⁻¹{{{tm}}}"), gi))
            g += gi
        g = sp.simplify(g)
        if a != 0:
            res.text("Segundo teorema de traslación: ℒ⁻¹{e^(−as)G(s)} = g(t−a)·u(t−a)")
            g = g.subs(t, t - a) * sp.Heaviside(t - a)
            res.math(g)
        total += g
    return sp.simplify(total)


def inverse(text):
    res = S.Result()
    res.method = "Transformada inversa de Laplace"
    F = parse_s(text)
    res.title("Función")
    res.math(sp.Eq(sp.Symbol('F(s)'), F))
    f = inverse_steps(F, res)
    res.title("Resultado")
    res.math(sp.Eq(sp.Symbol('f(t)'), f))
    res.plot = [('f(t)', f)]
    return res


# --------------------------------------------------------------- PVI con Laplace

def ivp(text, ics_text):
    res = S.Result()
    res.method = "PVI con transformada de Laplace"
    kind, (expr, n) = S.parse_equation(text, _EXTRA)
    if kind != 'ode':
        raise S.ParseError("Escribe la ecuación con y', y'', …")
    expr = expr.subs(S.x, t)
    a = [sp.diff(expr, S.YP[k]) for k in range(n + 1)]
    if any(c.free_symbols & {t, *S.YP} for c in a):
        raise S.ParseError("Laplace aquí requiere coeficientes constantes.")
    g = sp.simplify(-(expr - sum(a[k] * S.YP[k] for k in range(n + 1))))
    if g.free_symbols & set(S.YP):
        raise S.ParseError("La ecuación debe ser lineal.")

    ics = {k: v for k, x0, v in S.parse_ics(ics_text) if x0 == 0}
    missing = [k for k in range(n) if k not in ics]
    res.title("Ecuación")
    res.math(sp.Eq(sum(a[k] * S.dsym(k) for k in range(n + 1)), g))
    if missing:
        res.warn("Condiciones en t = 0 no dadas se toman como 0: "
                 + ", ".join(f"y{chr(39) * k}(0)" for k in missing))
    y0 = [ics.get(k, sp.Integer(0)) for k in range(n)]

    res.title("1) Transformamos cada término")
    res.text("ℒ{y⁽ᵏ⁾} = sᵏY(s) − sᵏ⁻¹y(0) − … − y⁽ᵏ⁻¹⁾(0):")
    LHS = 0
    for k in range(n + 1):
        if a[k] == 0:
            continue
        Lk = s ** k * Ys - sum(s ** (k - 1 - j) * y0[j] for j in range(k))
        res.math(sp.Eq(_lsym(S.dsym(k)), Lk))
        LHS += a[k] * Lk
    G = 0
    if g != 0:
        res.text("Lado derecho:")
        sub = S.Result()
        G = forward_steps(g, sub)
        res.math(sp.Eq(_lsym(g), G))
    res.title("2) Ecuación algebraica en Y(s)")
    res.math(sp.Eq(sp.expand(LHS), G))
    Y = sp.factor(sp.solve(sp.Eq(LHS, G), Ys)[0])
    res.text("Despejando:")
    res.math(sp.Eq(sp.Symbol('Y(s)'), Y))
    res.title("3) Transformada inversa")
    y = inverse_steps(Y, res)
    res.title("Solución")
    res.math(sp.Eq(sp.Symbol('y(t)'), y))

    res.title("Verificación")
    yfun = sp.Function('y')
    ode = sp.Eq(sum(a[k] * yfun(t).diff(t, k) for k in range(n + 1)), g)
    jumps = y.has(sp.Heaviside) or g.has(sp.Heaviside, sp.DiracDelta)
    if jumps:
        # con u(t−a) o δ(t−a) se comprueba el residuo en puntos lejos de los saltos
        resid = sum(a[k] * sp.diff(y, t, k) for k in range(n + 1)) - g
        pts = [0.37, 1.41, 2.73, 3.89, 5.17, 7.03]
        ok = all(abs(complex(sp.N(resid.subs(t, p)))) < 1e-8 for p in pts)
    else:
        try:
            ok = sp.checkodesol(ode, sp.Eq(yfun(t), y))[0]
        except Exception:
            ok = False
    ok_ic = all(sp.simplify(sp.diff(y, t, k).subs(t, 0) - y0[k]) == 0 for k in range(n))
    if ok and ok_ic:
        res.ok("✔ Cumple la ecuación" + (" (comprobada numéricamente fuera de los saltos)" if jumps else "")
               + " y las condiciones iniciales.")
    elif ok_ic:
        res.ok("✔ Cumple las condiciones iniciales.")
        res.warn("La ecuación no se pudo verificar.")
    else:
        res.warn("No se pudo verificar.")
    res.plot = [('y(t)', y)] + ([('g(t)', g)] if g != 0 and not g.has(sp.DiracDelta) else [])
    return res


# -------------------------------------------------------------------- convolución

def convolution(f_text, g_text):
    res = S.Result()
    res.method = "Convolución (f ∗ g)(t)"
    f, g = parse_t(f_text), parse_t(g_text)
    res.title("Definición")
    res.text("(f ∗ g)(t) = ∫₀ᵗ f(τ)·g(t−τ) dτ")
    integrand = f.subs(t, tau) * g.subs(t, t - tau)
    res.math(sp.Integral(integrand, (tau, 0, t)))
    conv = sp.simplify(sp.integrate(integrand, (tau, 0, t)))
    res.math(sp.Eq(sp.Symbol('(f ∗ g)(t)'), conv))
    res.title("Teorema de convolución: ℒ{f ∗ g} = F(s)·G(s)")
    Ff, Gg = L(f), L(g)
    res.math(sp.Eq(sp.Symbol('F(s)'), Ff))
    res.math(sp.Eq(sp.Symbol('G(s)'), Gg))
    prod = sp.simplify(Ff * Gg)
    res.math(sp.Eq(sp.Symbol('F(s)·G(s)'), prod))
    if sp.simplify(L(conv) - prod) == 0:
        res.ok("✔ ℒ{f ∗ g} coincide con F(s)·G(s).")
    else:
        res.warn("No se pudo comprobar el teorema simbólicamente.")
    res.plot = [('f(t)', f), ('g(t)', g), ('(f∗g)(t)', conv)]
    return res
