"""
Figura combinada: DFT (SIESTA, lineas grises) + TB exacto (diagonalizacion
directa del modelo, lineas azules) + VQD en simulación exacta (puntos de colores, uno
por nivel), todo en una sola grafica.
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

VQD_JSON = 'resultados/bandas_vqd11.json'
FIT_RESULT_NPZ = 'datos/fit_result.npz'
KPATH_NPZ = 'datos/kpath.npz'
SALIDA = 'bandas_completa.png'
ERROR_MAX = None
SOLAP_MAX = None


def representar(vqd_json=VQD_JSON, fit_result_npz=FIT_RESULT_NPZ, kpath_npz=KPATH_NPZ,
                 salida=SALIDA, error_max=ERROR_MAX, solap_max=SOLAP_MAX):

    d = np.load(kpath_npz)
    kcart, kdist_completo, target_bands = d['kcart'], d['kdist'], d['target_bands']
    params = np.load(fit_result_npz)['params']
    bands_tb = m.bands_at_k(kcart, params)

    with open(vqd_json) as f:
        resultados = json.load(f)['resultados_por_punto']

    n_niveles = max((len(p['niveles']) for p in resultados), default=0)
    n_incompletos = sum(1 for p in resultados if len(p['niveles']) < n_niveles)
    if n_incompletos:
        print(f"Aviso: {n_incompletos} punto(s) con menos de {n_niveles} niveles.")

    n_excluidos = 0

    fig, ax = plt.subplots(figsize=(9, 6))

    # DFT: 
    for b in range(target_bands.shape[1]):
        ax.plot(kdist_completo, target_bands[:, b], color='#999999', lw=1.2,
                zorder=1, label='DFT (SIESTA)' if b == 0 else None)

    # TB exacto:
    for b in range(bands_tb.shape[1]):
        ax.plot(kdist_completo, bands_tb[:, b], color='#1874CD', lw=1.1,
                zorder=2, label='TB exacto (diagonalización)' if b == 0 else None)

    # VQD:
    cmap = plt.get_cmap('tab20')
    for nivel in range(n_niveles):
        xs, ys = [], []
        for p in resultados:
            if nivel >= len(p['niveles']):
                continue
            n = p['niveles'][nivel]
            if error_max is not None and n['error'] > error_max:
                n_excluidos += 1
                continue
            if solap_max is not None and n['solapamiento_maximo'] > solap_max:
                n_excluidos += 1
                continue
            idx = int(p['etiqueta'])
            xs.append(kdist_completo[idx])
            ys.append(n['energia_vqd'])
        if not xs:
            continue
        ax.scatter(xs, ys, color=cmap(nivel % 20), s=20, zorder=5,
                   edgecolor='none', label=f'VQD nivel {nivel}')

    if n_excluidos:
        print(f"Excluidos {n_excluidos} puntos por superar error_max/solap_max.")

    xticks = [0.0, 0.605566, 1.304813, 1.654436]
    ax.set_xticks(xticks)
    ax.set_xticklabels(['M', r'$\Gamma$', 'K', 'M'])
    for x in xticks:
        ax.axvline(x, color='black', lw=0.7)
    ax.axhline(0, color='gray', lw=0.6, ls=':')
    ax.set_xlim(kdist_completo.min(), kdist_completo.max())
    ax.set_ylabel('E - E$_F$ (eV)')
    ax.set_xlabel('k-path')
    ax.set_title(f'DFT, TB exacto y VQD')
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=8,
              framealpha=0.9, borderaxespad=0)

    plt.tight_layout()
    plt.savefig(salida, dpi=150, bbox_inches='tight')
    print(f'Guardado: {salida}  ({len(resultados)} puntos, {n_niveles} niveles/punto)')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--vqd', default=VQD_JSON)
    p.add_argument('--fit-result', default=FIT_RESULT_NPZ)
    p.add_argument('--kpath', default=KPATH_NPZ)
    p.add_argument('--salida', default=SALIDA)
    p.add_argument('--error-max', type=float, default=ERROR_MAX)
    p.add_argument('--solap-max', type=float, default=SOLAP_MAX)
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    representar(vqd_json=args.vqd, fit_result_npz=args.fit_result, kpath_npz=args.kpath,
                salida=args.salida, error_max=args.error_max, solap_max=args.solap_max)
