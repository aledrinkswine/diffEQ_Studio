"""Imprime los pasos de una ecuación: python show_steps.py "y' = x*y" "y(0)=1" """
import sys

import solvers as S

r = S.solve(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else '')
print('Método:', r.method)
for k, c in r.steps:
    print(S.pretty(c) if k == 'math' else c)
