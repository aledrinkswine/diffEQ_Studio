"""Prueba rápida del motor con los ejemplos de la app."""
import sys
import time

import solvers as S
from examples import FIRST_ORDER, HIGHER_ORDER

verbose = '-v' in sys.argv
fails = 0
for eq, ics in FIRST_ORDER + HIGHER_ORDER:
    t = time.time()
    try:
        r = S.solve(eq, ics)
        ok = any(k == 'ok' for k, _ in r.steps) or r.implicit is not None
        status = 'OK ' if ok else '?? '
        detail = f"{r.method} | part={r.particular if r.particular is not None else r.particular_implicit}"
    except Exception as e:
        status, detail = 'ERR', repr(e)
        fails += 1
    print(f"{status} {time.time() - t:5.1f}s  {eq}  [{ics}]  -> {detail}")
    if verbose and status != 'ERR':
        for k, c in r.steps:
            print('   ', S.pretty(c) if k == 'math' else c)
sys.exit(1 if fails else 0)
