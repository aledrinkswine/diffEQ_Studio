"""Prueba de humo de la interfaz: resuelve todos los ejemplos y guarda las gráficas."""
import os
import sys

import app as A

out = sys.argv[1] if len(sys.argv) > 1 else 'gui_test_out'
os.makedirs(out, exist_ok=True)
root = A.App()
nb = root.nametowidget(root.winfo_children()[-1])
first, higher, lap, sysm, numeric, field = [root.nametowidget(t) for t in nb.tabs()[:6]]

jobs = [(first, None, i) for i in range(len(first.examples))]
jobs += [(higher, None, i) for i in range(len(higher.examples))]
jobs += [(lap, key, i) for key, *_, ex in A.LAPLACE_MODES for i in range(len(ex))]
jobs += [(sysm, None, i) for i in range(len(A.SYSTEMS))]
errors = []


def step(k=0):
    if k < len(jobs):
        tab, mode, i = jobs[k]
        if mode and tab.mode.get() != mode:
            tab.mode.set(mode)
            tab.set_mode()
        tab.combo.current(i)
        tab.load_example()
        wait(tab, k)
    else:
        numeric.compute()
        root.after(4000, finish)


def wait(tab, k):
    if 'disabled' in tab.btn.state():
        root.after(100, wait, tab, k)
        return
    body = tab.text.get('1.0', 'end')
    graph = ' '.join(t.get_text() for ax in tab.plot.fig.axes for t in ax.texts)
    if 'No se pudo' in body or 'No pude' in body or 'No se pudo graficar' in graph:
        errors.append((k, body.strip()[:150], graph[:150]))
    tab.plot.fig.savefig(os.path.join(out, f'tab_{k:02d}.png'))
    step(k + 1)


def finish():
    numeric.plot.fig.savefig(os.path.join(out, 'numeric.png'))
    print('numeric:', numeric.status.cget('text'), '| rows:', len(numeric.tree.get_children()))
    field.on_click(type('E', (), {'inaxes': field.plot.ax, 'xdata': 0.0, 'ydata': 0.5})())
    field.plot.fig.savefig(os.path.join(out, 'field.png'))
    print(f'{len(jobs)} ejemplos, errores:', errors)
    root.destroy()


root.after(500, step)
root.mainloop()
