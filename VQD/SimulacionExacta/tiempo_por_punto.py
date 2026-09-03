"""
tiempo_por_punto.py

Grafica el tiempo de simulacion (clasica) por punto k, en orden de
ejecucion, para mostrar el efecto del warm-start: el primer punto
(arranque en frio, multi-semilla) tarda notablemente mas que el resto
(warm-start, un unico intento por nivel).

Nota metodologica importante: este es tiempo de SIMULACION CLASICA del
circuito, no tiempo de ejecucion en hardware cuantico real -- no debe
interpretarse como una estimacion de coste en un dispositivo cuantico.
Sirve para cuantificar el efecto relativo del warm-start dentro de la
misma metodologia de simulacion.

Uso:
    python3 tiempo_por_punto.py
"""
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

VQD_JSON = 'resultados/bandas_vqd11.json'
SALIDA = 'tiempo_por_punto.png'


def representar(vqd_json=VQD_JSON, salida=SALIDA):
    with open(vqd_json) as f:
        resultados = json.load(f)['resultados_por_punto']

    resultados_ordenados = sorted(resultados, key=lambda p: p['indice_punto'])
    indices = [p['indice_punto'] for p in resultados_ordenados]
    tiempos_min = [p['tiempo_total_punto_s'] / 60 for p in resultados_ordenados]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(indices, tiempos_min, 'o-', color='black', markersize=4, linewidth=1)
    ax.axvline(0.5, color='gray', lw=0.8, ls=':')
    ax.annotate('arranque en frio', xy=(0, tiempos_min[0]), xytext=(3, tiempos_min[0] - 2),
                fontsize=9, arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))
    ax.set_xlabel('Punto (orden de ejecucion)')
    ax.set_ylabel('Tiempo (minutos)')
    ax.set_title('Tiempo de simulacion por punto: arranque en frio vs. warm-start')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(salida, dpi=150, bbox_inches='tight')

    print(f'Punto 0 (frio): {tiempos_min[0]:.1f} min')
    print(f'Media resto (warm-start): {np.mean(tiempos_min[1:]):.1f} min')
    print(f'Factor de aceleracion: {tiempos_min[0] / np.mean(tiempos_min[1:]):.2f}x')
    print(f'Guardado: {salida}')


if __name__ == '__main__':
    representar()
