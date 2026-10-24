"""Prueba del módulo de Laplace y de sistemas con todos los ejemplos."""
import sys
import time

import laplace as LP
import solvers as S
import systems as SY
from examples import CONVOLUTION, INVERSE, IVP_LAPLACE, SYSTEMS, TRANSFORM

verbose = '-v' in sys.argv
cases = ([(LP.transform, (f,)) for f in TRANSFORM] + [(LP.inverse, (F,)) for F in INVERSE]
         + [(LP.ivp, c) for c in IVP_LAPLACE] + [(LP.convolution, c) for c in CONVOLUTION]
         + [(SY.solve, c) for c in SYSTEMS])
fails = 0
for fn, args in cases:
    t0 = time.time()
    try:
        r = fn(*args)
        warn = any(k == 'warn' for k, _ in r.steps)
        status = 'WRN' if warn else 'OK '
        last = next((c for k, c in reversed(r.steps) if k == 'math'), '')
        detail = f"{r.method} | {last}"
    except Exception as e:
        status, detail = 'ERR', repr(e)
        fails += 1
    print(f"{status} {time.time() - t0:5.1f}s  {fn.__module__}.{fn.__name__}{args}  -> {detail}")
    if verbose and status != 'ERR':
        for k, c in r.steps:
            print('   ', S.pretty(c) if k == 'math' else c)
sys.exit(1 if fails else 0)
