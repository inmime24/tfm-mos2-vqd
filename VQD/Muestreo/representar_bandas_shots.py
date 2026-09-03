"""
Figura combinada: TB exacto (gris, referencia continua) + VQD por muestreo en
la cadena de warm-start (verde) + VQD por muestreo en los puntos especiales
M/Gamma/K (morado).

Uso:
    python3 representar_figura2_ruido.py \
        --warmstart resultados/ruido_1_warmstart.json \
        --especiales resultados/ruido_2_M.json resultados/ruido_3_Gamma.json resultados/ruido_4_K.json \
        --salida figura2_ruido.png
"""
import json
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

import sk_model as m
from matplotlib.colors import LinearSegmentedColormap

WARMSTART_JSON = 'resultados/ruido_1_warmstart.json'
ESPECIALES_JSONS = ['resultados/ruido_2_M.json', 'resultados/ruido_3_Gamma.json', 'resultados/ruido_4_K.json']
FIT_RESULT_NPZ = 'datos/fit_result.npz'
KPATH_NPZ = 'datos/kpath.npz'
SALIDA = 'bandas_ruido.png'
UMBRAL_SALTO = 0.3  # fraccion del k-path total; saltos mayores no se conectan con linea
N_CICLOS_GRADIENTE = 2  # repeticiones oscuro->claro


def generar_colores_gradiente(color_oscuro, color_claro, n, n_ciclos=N_CICLOS_GRADIENTE):

    cmap = LinearSegmentedColormap.from_list('grad', [color_oscuro, color_claro])
    n_por_ciclo = int(np.ceil(n / n_ciclos))
    colores = []
    for i in range(n):
        ciclo_pos = i % n_por_ciclo
        frac = ciclo_pos / max(n_por_ciclo - 1, 1)
        colores.append(cmap(frac))
    return colores


def cargar_puntos_ordenados(json_path):
    with open(json_path) as f:
        puntos = json.load(f)['resultados_por_punto']
    return sorted(puntos, key=lambda p: p['indice_punto'])


def dibujar_dataset(ax, puntos, kdist_completo, colores_por_nivel, etiqueta_leyenda, umbral_salto):
    n_niveles = max((len(p['niveles']) for p in puntos), default=0)
    ya_etiquetado = False
    for nivel in range(n_niveles):
        xs, ys = [], []
        for p in puntos:
            if nivel < len(p['niveles']):
                idx = int(p['etiqueta'])
                xs.append(kdist_completo[idx])
                ys.append(p['niveles'][nivel]['energia_exacta_final'])
        if not xs:
            continue
        color = colores_por_nivel[nivel]
        # cortar la linea donde haya un salto grande en x (p.ej. el cierre en M)
        rango_total = kdist_completo.max() - kdist_completo.min()
        segmentos_x, segmentos_y = [[xs[0]]], [[ys[0]]]
        for i in range(1, len(xs)):
            if abs(xs[i] - xs[i-1]) > umbral_salto * rango_total:
                segmentos_x.append([])
                segmentos_y.append([])
            segmentos_x[-1].append(xs[i])
            segmentos_y[-1].append(ys[i])
        for seg_x, seg_y in zip(segmentos_x, segmentos_y):
            ax.plot(seg_x, seg_y, 'o-', color=color, markersize=4, lw=1.2, zorder=5,
                     label=etiqueta_leyenda if not ya_etiquetado else None)
            ya_etiquetado = True


def representar(warmstart_json=WARMSTART_JSON, especiales_jsons=ESPECIALES_JSONS,
                 fit_result_npz=FIT_RESULT_NPZ, kpath_npz=KPATH_NPZ, salida=SALIDA,
                 umbral_salto=UMBRAL_SALTO):

    d = np.load(kpath_npz)
    kcart, kdist_completo = d['kcart'], d['kdist']
    params = np.load(fit_result_npz)['params']
    bands_tb = m.bands_at_k(kcart, params)
    n_bandas = bands_tb.shape[1]

    colores_verdes = generar_colores_gradiente('darkgreen', 'chartreuse', n_bandas)
    colores_morados = generar_colores_gradiente('indigo', 'darkorchid', n_bandas)

    fig, ax = plt.subplots(figsize=(9, 6))

    for b in range(n_bandas):
        ax.plot(kdist_completo, bands_tb[:, b], color='#999999', lw=1.2, zorder=1,
                 label='TB exacto (diagonalización)' if b == 0 else None)

    puntos_warmstart = cargar_puntos_ordenados(warmstart_json)
    dibujar_dataset(ax, puntos_warmstart, kdist_completo, colores_verdes,
                     'VQD con ruido (warm-start)', umbral_salto)

   
    ya_etiquetado_morado = False
    for archivo in especiales_jsons:
        puntos = cargar_puntos_ordenados(archivo)
        dibujar_dataset(ax, puntos, kdist_completo, colores_morados,
                         'VQD con ruido (puntos especiales)' if not ya_etiquetado_morado else '_nolegend_',
                         umbral_salto)
        ya_etiquetado_morado = True

    xticks = [0.0, 0.605566, 1.304813, 1.654436]
    ax.set_xticks(xticks)
    ax.set_xticklabels(['M', r'$\Gamma$', 'K', 'M'])
    for x in xticks:
        ax.axvline(x, color='black', lw=0.7)
    ax.axhline(0, color='gray', lw=0.6, ls=':')
    ax.set_xlim(kdist_completo.min(), kdist_completo.max())
    ax.set_ylabel('E - E$_F$ (eV)')
    ax.set_xlabel('k-path')
    ax.set_title('TB exacto vs. VQD con muestreo')
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(salida, dpi=150, bbox_inches='tight')
    print(f'Guardado: {salida}')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--warmstart', default=WARMSTART_JSON)
    p.add_argument('--especiales', nargs='+', default=ESPECIALES_JSONS)
    p.add_argument('--fit-result', default=FIT_RESULT_NPZ)
    p.add_argument('--kpath', default=KPATH_NPZ)
    p.add_argument('--salida', default=SALIDA)
    p.add_argument('--umbral-salto', type=float, default=UMBRAL_SALTO)
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    representar(warmstart_json=args.warmstart, especiales_jsons=args.especiales,
                fit_result_npz=args.fit_result, kpath_npz=args.kpath,
                salida=args.salida, umbral_salto=args.umbral_salto)
