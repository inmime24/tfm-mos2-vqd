"""
Ajusta los 12 parametros del modelo tight-binding (sk_model.py) para que
sus bandas reproduzcan las de DFT mediante minimos cuadrados no lineales
(least_squares, metodo TRF) con varios arranques aleatorios, quedandose
con el mejor resultado.

Entrada: kpath.npz (puntos k y bandas de referencia).
Salida: fit_result.npz (los 12 parametros ajustados, con sus nombres).
"""


import numpy as np
from scipy.optimize import least_squares
import sk_model as m

#CARGA DE DATOS
d = np.load('kpath.npz')
kcart = d['kcart'] # (106) valores de k
target = d['target_bands']  # (106, 11), energias relativas a Ef

#FUNCIÓN DE COSTE
def residuals(params):
    bands = m.bands_at_k(kcart, params)
    return (bands - target).ravel()

#COTAS INFERIORES Y SUPERIORES PARA CADA PARÁMETRO
lb = np.array([-8,-8,-8, -10,-10,  -6,-6,  -3,-3,-3,  -3,-3])
ub = np.array([ 4, 4, 4,  -1,-1,   0, 6,   3, 3, 3,   3, 3])

#OPTIMIZACIÓN
rng = np.random.default_rng(42)
best = None
for trial in range(6):
    p0 = rng.uniform(lb, ub)
    res = least_squares(residuals, p0, bounds=(lb, ub), method='trf',
                         xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=4000)
    cost = res.cost
    print(f'trial {trial}: cost={cost:.4f}  rmse={np.sqrt(2*cost/target.size):.4f} eV')
    if best is None or cost < best.cost:
        best = res

#RESULTADOS
print('\nBest fit RMSE:', np.sqrt(2*best.cost/target.size), 'eV')
for name, val in zip(m.PARAM_NAMES, best.x):
    print(f'  {name:12s} = {val:+.4f} eV')

np.savez('fit_result.npz', params=best.x, param_names=m.PARAM_NAMES)
