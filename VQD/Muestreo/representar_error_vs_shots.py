"""
Grafica el RMSE (y el error maximo) de VQD con muestreo frente a la
diagonalizacion exacta, en el punto K, en funcion del numero de disparos
(shots) por termino de Pauli. Eje X en escala logaritmica.

Uso:
    python3 representar_error_vs_shots.py --salida error_vs_shots.png
"""
import json
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

_cmap_verde = LinearSegmentedColormap.from_list('grad', ['darkgreen', 'chartreuse'])
COLOR_RMSE = _cmap_verde(0.0)      # verde oscuro
COLOR_MEDIO = _cmap_verde(0.5)     # verde medio
COLOR_MAXIMO = _cmap_verde(1.0)    # verde claro
COLOR_TIEMPO = 'darkorchid'

ARCHIVOS = {
    4096: 'resultados/ruido_5_shots4096.json',
    8192: 'resultados/ruido_8_shots8192_K.json',
    16384: 'resultados/ruido_9_shots16384_K.json',
    32768: 'resultados/ruido_6_shots32768.json',
    65536: 'resultados/ruido_7_shots65536.json',
}
SALIDA = 'error_vs_shots.png'


def representar(archivos=ARCHIVOS, salida=SALIDA):
    shots_list = sorted(archivos.keys())
    rmses, medios, maximos, tiempos = [], [], [], []

    for shots in shots_list:
        with open(archivos[shots]) as f:
            d = json.load(f)
        punto = d['resultados_por_punto'][0]
        errores = np.array([n['error_vs_exacto'] for n in punto['niveles']])
        rmses.append(np.sqrt(np.mean(errores**2)))
        medios.append(errores.mean())
        maximos.append(errores.max())
        tiempos.append(punto['tiempo_total_punto_s'] / 60)

    fig, ax1 = plt.subplots(figsize=(7.5, 5.5))
    ax2 = ax1.twinx()

    l1, = ax1.plot(shots_list, rmses, 'o-', color=COLOR_RMSE, markersize=6, lw=1.3, label='RMSE')
    l2, = ax1.plot(shots_list, medios, '^-', color=COLOR_MEDIO, markersize=6, lw=1.2, label='Error medio')
    l3, = ax1.plot(shots_list, maximos, 's-', color=COLOR_MAXIMO, markersize=6, lw=1.1, label='Error máximo')
    l4, = ax2.plot(shots_list, tiempos, 'D-', color=COLOR_TIEMPO, markersize=5, lw=1.2, label='Tiempo total')

    ax1.set_xscale('log', base=2)
    ax1.set_xticks(shots_list)
    ax1.set_xticklabels([str(s) for s in shots_list])
    ax1.set_xlabel('Shots por término de Pauli')
    ax1.set_ylabel('Error frente a diagonalización exacta (eV)')
    ax2.set_ylabel('Tiempo total del punto (min)')
    ax2.tick_params(axis='y')
    ax1.set_title('Error y tiempo de VQD vs. número de shots')
    ax1.grid(True, alpha=0.3)
    ax1.legend(handles=[l1, l2, l3, l4], fontsize=9, loc='upper left')

    plt.tight_layout()
    plt.savefig(salida, dpi=150, bbox_inches='tight')
    print(f'Guardado: {salida}')
    for shots, rmse, medio, maximo, tiempo in zip(shots_list, rmses, medios, maximos, tiempos):
        print(f'  {shots:>6} shots: RMSE={rmse:.4f}  medio={medio:.4f}  max={maximo:.4f} eV  tiempo={tiempo:.1f} min')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--salida', default=SALIDA)
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    representar(salida=args.salida)
