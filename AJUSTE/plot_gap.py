"""
Compara graficamente las bandas del modelo tight-binding ajustado frente a
las de DFT (SIESTA), a lo largo del camino M-Gamma-K-M, y calcula el gap
directo en K para las dos.

Entrada: kpath.npz (bandas de referencia), fit_result.npz (parametros ajustados).
Salida: gap_comparison.png.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sk_model as m

#CARGA DE DATOS
d = np.load('kpath.npz')
kdist, kcart, target = d['kdist'], d['kcart'], d['target_bands']
r = np.load('fit_result.npz')
params = r['params']
bands = m.bands_at_k(kcart, params)

n_bands = target.shape[1]
COLOR_DFT = '#BBBBBB'      
COLOR_TB = '#B5342A'       

#FIGURA
fig, ax = plt.subplots(figsize=(7.5, 6))

for b in range(n_bands):
    ax.plot(kdist, target[:, b], color=COLOR_DFT, lw=1.8, zorder=1,
             label='DFT (SIESTA)' if b == 0 else None)
    ax.plot(kdist, bands[:, b], color=COLOR_TB, lw=1.3, zorder=2,
             label='TB ajustado' if b == 0 else None)

# --- nivel de fermi ---
ax.axhline(0.0, color='gray', lw=0.9, ls='--', label='$E_F$', zorder=0) 

# --- gap directo en K: DFT vs TB ---
idx_K = int(np.argmin(np.abs(kdist - 1.304813)))
xK = kdist[idx_K]

E_dft_K, E_tb_K = target[idx_K], bands[idx_K]
vbm_dft, cbm_dft = E_dft_K[E_dft_K < 0].max(), E_dft_K[E_dft_K > 0].min()
vbm_tb, cbm_tb = E_tb_K[E_tb_K < 0].max(), E_tb_K[E_tb_K > 0].min()
gap_dft, gap_tb = cbm_dft - vbm_dft, cbm_tb - vbm_tb

off = 0.035
# DFT:
ax.annotate('', xy=(xK - off, cbm_dft), xytext=(xK - off, vbm_dft),
            arrowprops=dict(arrowstyle='<->', color='#666666', lw=1.4))
ax.plot(xK, vbm_dft, 'o', color='#666666', ms=6, zorder=5)
ax.plot(xK, cbm_dft, 'o', color='#666666', ms=6, zorder=5)
ax.text(xK - off - 0.05, (vbm_dft + cbm_dft) / 2, f'{gap_dft:.3f} eV\n(DFT)',
        ha='right', va='center', fontsize=9.5, color='#555555')

# TB: 
ax.annotate('', xy=(xK + off, cbm_tb), xytext=(xK + off, vbm_tb),
            arrowprops=dict(arrowstyle='<->', color=COLOR_TB, lw=1.4))
ax.plot(xK, vbm_tb, 'o', color=COLOR_TB, ms=6, zorder=5)
ax.plot(xK, cbm_tb, 'o', color=COLOR_TB, ms=6, zorder=5)
ax.text(xK + off + 0.05, (vbm_tb + cbm_tb) / 2, f'{gap_tb:.3f} eV\n(TB)',
        ha='left', va='center', fontsize=9.5, color=COLOR_TB)

# --- ejes, etiquetas de alta simetria ---
xticks = [0.0, 0.605566, 1.304813, 1.654436]
ax.set_xticks(xticks)
ax.set_xticklabels(['M', r'$\Gamma$', 'K', 'M'])
for x in xticks:
    ax.axvline(x, color='black', lw=0.8)

ax.set_xlim(kdist.min(), kdist.max())
ax.set_ylim(-3, 3)
ax.set_ylabel('Energy (eV)')
ax.set_xlabel('k-path')
ax.set_title('MoS$_2$ monocapa: DFT vs. TB ajustado')
ax.legend(loc='upper right', fontsize=9.5)

plt.tight_layout()
plt.savefig('gap_comparison.png', dpi=150)
print('guardado. DFT gap(K)=%.4f eV   TB gap(K)=%.4f eV' % (gap_dft, gap_tb))
