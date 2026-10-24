"""Sistemas de EDO: valores/vectores propios, variación de parámetros, linealización."""
import re

import sympy as sp

import solvers as S

t = sp.Symbol('t')


def _eq(lhs, rhs):
    return sp.Eq(lhs, rhs, evaluate=False)


_LINE = re.compile(r"^\s*(?:d([a-zA-Z])\s*/\s*dt|([a-zA-Z])\s*['′’])\s*=\s*(.+)$")
_IC = re.compile(r"([a-zA-Z])\s*\(\s*([^)]+?)\s*\)\s*=\s*([^,;\n]+)")


def parse_system(text):
    lines = [l for l in re.split(r"[\n;]+", text) if l.strip()]
    if not lines:
        raise S.ParseError("Escribe el sistema, una ecuación por línea: x' = …")
    names, rhs_txt = [], []
    for l in lines:
        m = _LINE.match(l)
        if not m:
            raise S.ParseError(f"No entiendo la línea '{l.strip()}'. Usa el formato x' = …")
        names.append(m.group(1) or m.group(2))
        rhs_txt.append(m.group(3))
    if len(set(names)) != len(names) or 't' in names or 'e' in names:
        raise S.ParseError("Usa letras distintas para las variables (no t ni e).")
    vars_ = [sp.Symbol(n) for n in names]
    extra = dict(zip(names, vars_), t=t)
    rhs = [sp.nsimplify(S._p(r, extra), rational=True) for r in rhs_txt]
    unknown = set().union(*(r.free_symbols for r in rhs)) - set(vars_) - {t}
    if unknown:
        raise S.ParseError(f"Variables desconocidas: {', '.join(map(str, unknown))}")
    return names, vars_, rhs


def parse_sys_ics(text, names):
    out = {}
    for m in _IC.finditer(text or ''):
        if m.group(1) in names:
            out[m.group(1)] = (sp.nsimplify(S.parse_value(m.group(2)), rational=True),
                               sp.nsimplify(S.parse_value(m.group(3).strip()), rational=True))
    if (text or '').strip() and not out:
        raise S.ParseError("Condiciones iniciales: usa el formato x(0)=1, y(0)=0")
    return out


def _clean(v):
    """Escala un vector para quitar denominadores enteros."""
    v = sp.simplify(v)
    dens = [sp.fraction(sp.together(c))[1] for c in v]
    if all(d.is_Integer for d in dens):
        v = v * sp.ilcm(*[int(d) for d in dens], 1)
    return v


def classify(eigs, A=None):
    """Clasifica el punto crítico de un sistema 2×2 a partir de sus valores propios."""
    l1, l2 = (complex(sp.N(e)) for e in eigs)
    eps = 1e-12
    if abs(l1.imag) > eps:
        a = l1.real
        if abs(a) < eps:
            return "Centro", "estable (no asintóticamente)", "Trayectorias cerradas (órbitas periódicas)."
        if a < 0:
            return "Foco (espiral) estable", "asintóticamente estable", "Re(λ) < 0: las espirales se acercan al origen."
        return "Foco (espiral) inestable", "inestable", "Re(λ) > 0: las espirales se alejan."
    a, b = l1.real, l2.real
    if abs(a) < eps or abs(b) < eps:
        return "Caso degenerado (λ = 0)", "no aislado", "Hay una recta de puntos de equilibrio."
    if a * b < 0:
        return "Punto silla", "inestable", "λ reales de signo opuesto."
    stab = "asintóticamente estable" if a < 0 else "inestable"
    kind = "estable" if a < 0 else "inestable"
    if abs(a - b) < eps:
        if A is not None and sp.simplify(A - A[0, 0] * sp.eye(2)) == sp.zeros(2):
            return f"Nodo propio (estrella) {kind}", stab, "λ repetido con dos vectores propios."
        return f"Nodo impropio {kind}", stab, "λ repetido con un solo vector propio."
    return f"Nodo {kind}", stab, "λ reales distintos del mismo signo."


def _eig_list(M):
    out = []
    for v, m in M.eigenvals().items():
        out += [v] * m
    return out


def _linear(A, Fv, vars_, ics, res):
    n = A.shape[0]
    nonhom = Fv != sp.zeros(n, 1)
    res.method = "Sistema lineal — valores propios" + (" + variación de parámetros" if nonhom else "")
    res.title("Forma matricial X' = A·X" + (" + F(t)" if nonhom else ""))
    res.math(_eq(sp.Symbol('X'), sp.Matrix(vars_)))
    res.math(_eq(sp.Symbol('A'), A))
    if nonhom:
        res.math(_eq(sp.Symbol('F(t)'), Fv))

    res.title("1) Valores propios")
    lam = sp.Symbol('λ')
    cp = sp.expand((A - lam * sp.eye(n)).det())
    res.text("Ecuación característica det(A − λI) = 0:")
    res.math(sp.Eq(cp, 0))
    if n == 2:
        tr, dt = A.trace(), A.det()
        res.text(f"traza(A) = {tr},  det(A) = {dt},  tr² − 4·det = {sp.simplify(tr ** 2 - 4 * dt)}")

    res.title("2) Vectores propios y soluciones fundamentales")
    sols, real_dirs = [], []
    for val, mult, vecs in A.eigenvects():
        val = sp.simplify(val)
        re_, im_ = (sp.simplify(p) for p in val.as_real_imag())
        if im_ != 0:
            if im_.is_negative:
                continue
            v = _clean(vecs[0])
            a = v.applyfunc(lambda c: sp.simplify(sp.re(sp.expand_complex(c))))
            b = v.applyfunc(lambda c: sp.simplify(sp.im(sp.expand_complex(c))))
            res.text(f"λ = {val}   (α = {re_}, β = {im_})")
            res.math(_eq(sp.Symbol('v'), v))
            res.text("Escribimos v = a + b·i:")
            res.math(_eq(sp.Symbol('a'), a))
            res.math(_eq(sp.Symbol('b'), b))
            res.text("Las partes real e imaginaria de v·e^(λt) son soluciones reales:")
            X1 = sp.exp(re_ * t) * (a * sp.cos(im_ * t) - b * sp.sin(im_ * t))
            X2 = sp.exp(re_ * t) * (a * sp.sin(im_ * t) + b * sp.cos(im_ * t))
            for X in (X1, X2):
                sols.append(sp.simplify(X))
                res.math(_eq(sp.Symbol(f"X{len(sols)}"), sols[-1]))
            continue
        res.text(f"λ = {val}" + (f"   (multiplicidad {mult})" if mult > 1 else ""))
        for v in vecs:
            v = _clean(v)
            real_dirs.append((val, v))
            res.math(_eq(sp.Symbol('v'), v))
            sols.append(v * sp.exp(val * t))
            res.math(_eq(sp.Symbol(f"X{len(sols)}"), sols[-1]))
        if len(vecs) < mult:
            if mult - len(vecs) != 1:
                return _expm(A, Fv, vars_, ics, res)
            v = _clean(vecs[0])
            res.text("Solo hay un vector propio: buscamos un vector generalizado w con (A − λI)·w = v")
            sol, params = (A - val * sp.eye(n)).gauss_jordan_solve(v)
            w = sol.subs({p: 0 for p in params})
            res.math(_eq(sp.Symbol('w'), w))
            sols.append(sp.exp(val * t) * (t * v + w))
            res.math(_eq(sp.Symbol(f"X{len(sols)}"), sols[-1]))
    return _assemble(A, Fv, vars_, ics, res, sols, real_dirs)


def _expm(A, Fv, vars_, ics, res):
    res.text("Valor propio con multiplicidad alta: usamos la matriz exponencial e^(At).")
    Phi = sp.simplify((A * t).exp())
    res.math(_eq(sp.Symbol('e^(At)'), Phi))
    return _assemble(A, Fv, vars_, ics, res, [Phi[:, j] for j in range(A.shape[0])], [])


def _assemble(A, Fv, vars_, ics, res, sols, real_dirs):
    n = A.shape[0]
    Cs = sp.symbols(f'C1:{n + 1}')
    Xh = sp.zeros(n, 1)
    for c, X in zip(Cs, sols):
        Xh += c * X
    res.title("Solución " + ("homogénea" if Fv != sp.zeros(n, 1) else "general"))
    res.math(_eq(sp.Symbol('X_h'), sp.Add(*[sp.MatMul(c, X, evaluate=False) for c, X in zip(Cs, sols)],
                                             evaluate=False)))
    Xp = sp.zeros(n, 1)
    if Fv != sp.zeros(n, 1):
        res.title("3) Solución particular — variación de parámetros")
        Phi = sp.Matrix.hstack(*sols)
        res.text("Matriz fundamental Φ(t) = [X₁ X₂ …]:")
        res.math(_eq(sp.Symbol('Φ(t)'), Phi))
        Pinv = sp.simplify(Phi.inv())
        res.math(_eq(sp.Symbol('Φ⁻¹(t)'), Pinv))
        up = sp.simplify(Pinv * Fv)
        res.text("u'(t) = Φ⁻¹(t)·F(t):")
        res.math(_eq(sp.Symbol("u'"), up))
        u = up.applyfunc(lambda e: sp.integrate(e, t))
        res.math(_eq(sp.Symbol('u'), u))
        Xp = sp.simplify(Phi * u)
        res.text("X_p = Φ(t)·u(t):")
        res.math(_eq(sp.Symbol('X_p'), Xp))
        res.title("Solución general X = X_h + X_p")
    X = Xh + Xp
    for v, comp in zip(vars_, X):
        res.math(_eq(sp.Symbol(f"{v}(t)"), sp.collect(sp.expand(comp), Cs)))

    res.title("Verificación")
    if sp.simplify(X.diff(t) - A * X - Fv) == sp.zeros(n, 1):
        res.ok("✔ X' = A·X" + (" + F(t)" if Fv != sp.zeros(n, 1) else "") + " se cumple.")
    else:
        res.warn("No se pudo verificar simbólicamente.")

    part = None
    if ics:
        res.title("Problema de valor inicial")
        eqs = [sp.Eq(comp.subs(t, ics[str(v)][0]), ics[str(v)][1])
               for v, comp in zip(vars_, X) if str(v) in ics]
        for e in eqs:
            res.math(e)
        sol = sp.solve(eqs, Cs, dict=True)
        if sol:
            for c, val in sol[0].items():
                res.math(sp.Eq(c, val))
            part = sp.expand(X.subs(sol[0]))
            res.text("Solución particular:")
            for v, comp in zip(vars_, part):
                res.math(_eq(sp.Symbol(f"{v}(t)"), comp))
        else:
            res.warn("No se pudieron determinar las constantes.")

    info = {'general': X, 'particular': part, 'dirs': real_dirs, 'equilibria': [], 'class': None}
    if n == 2 and not Fv.has(t):
        name, stab, why = classify(_eig_list(A), A)
        eq = sp.zeros(2, 1) if Fv == sp.zeros(2, 1) else -A.inv() * Fv if A.det() != 0 else None
        res.title("Clasificación del punto crítico")
        if eq is not None:
            res.text(f"Punto de equilibrio: ({eq[0]}, {eq[1]})")
            info['equilibria'] = [tuple(float(c) for c in eq)]
        res.ok(f"{name} — {stab}")
        res.text(why)
        info['class'] = name
    return info


def _nonlinear(rhs, vars_, res):
    n = len(vars_)
    res.method = "Sistema no lineal — linealización"
    res.title("Sistema")
    for v, r in zip(vars_, rhs):
        res.math(_eq(sp.Symbol(f"{v}'"), r))
    info = {'general': None, 'particular': None, 'dirs': [], 'equilibria': [], 'class': None}
    if any(r.has(t) for r in rhs):
        res.warn("El sistema no es autónomo: se resuelve solo numéricamente (RK4).")
        return info
    res.title("1) Puntos de equilibrio (todas las derivadas = 0)")
    try:
        eqs = sp.solve(rhs, vars_, dict=True)
    except Exception:
        eqs = []
    eqs = [p for p in eqs if len(p) == n and all(val.is_real for val in p.values())][:12]
    if not eqs:
        res.warn("No se encontraron equilibrios de forma exacta.")
    for p in eqs:
        res.math(sp.Tuple(*[p[v] for v in vars_]))
    res.title("2) Matriz jacobiana")
    J = sp.Matrix(rhs).jacobian(vars_)
    res.math(_eq(sp.Symbol('J'), J))
    res.title("3) Clasificación por linealización")
    for p in eqs:
        Jp = sp.simplify(J.subs(p))
        pt = tuple(p[v] for v in vars_)
        res.text(f"En {pt}:")
        res.math(_eq(sp.Symbol('J'), Jp))
        ev = _eig_list(Jp)
        res.text("λ = " + ", ".join(str(sp.simplify(e)) for e in ev))
        if n == 2:
            name, stab, _ = classify(ev, Jp)
            res.ok(f"{name} — {stab}")
        else:
            reals = [complex(sp.N(e)).real for e in ev]
            res.ok("asintóticamente estable" if max(reals) < 0 else
                   "inestable" if max(reals) > 0 else "no concluyente")
        info['equilibria'].append(tuple(float(c) for c in pt))
    res.text("Nota: si la linealización da un centro o λ = 0, el criterio no es concluyente"
             " para el sistema no lineal.")
    return info


def solve(text, ics_text=''):
    res = S.Result()
    names, vars_, rhs = parse_system(text)
    ics = parse_sys_ics(ics_text, names)
    n = len(vars_)
    A = sp.Matrix([[sp.diff(r, v) for v in vars_] for r in rhs])
    linear = all(not (a.free_symbols & (set(vars_) | {t})) for a in A)
    if linear:
        Fv = sp.Matrix([sp.simplify(r - sum(A[i, j] * vars_[j] for j in range(n)))
                        for i, r in enumerate(rhs)])
        info = _linear(A, Fv, vars_, ics, res)
    else:
        info = _nonlinear(rhs, vars_, res)
        if ics:
            res.text("La trayectoria con las condiciones iniciales se calcula numéricamente (RK4) en la gráfica.")
    res.sys = dict(info, vars=vars_, names=names, rhs=rhs, ics=ics, linear=linear)
    return res
