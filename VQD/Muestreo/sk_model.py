"""
Modelo tight-binding Slater-Koster para MoS2 monocapa (11 orbitales: 5 d de Mo + 3 p de S_top + 3 p de S_bottom).

Orden de orbitales d (Mo): xy, yz, zx, x2-y2, 3z2-r2 (son los utilizados en las tablas Slater-Koster)
Orden de orbitales p (S):  px, py, pz
Base total (11): [Mo: xy,yz,zx,x2-y2,3z2-r2, S1: px,py,pz, S2: px,py,pz]
"""
import numpy as np

# ---------------- Geometria real (desde mos2_bands.fdf) ----------------
ALAT = 3.170
A1 = ALAT*np.array([1.0, 0.0, 0.0])
A2 = ALAT*np.array([-0.5, 0.866025, 0.0])
A3 = ALAT*np.array([0.0, 0.0, 3.890823])
LATVECS = np.array([A1, A2, A3])

FRAC = {
    'Mo': np.array([0.0, 0.0, 0.26800061]),
    'S1': np.array([2/3, 1/3, 0.13886276]),
    'S2': np.array([2/3, 1/3, 0.39713732]),
}
CART = {k: FRAC[k] @ LATVECS for k in FRAC}


ORBS_PER_ATOM = {'Mo': 5, 'S1': 3, 'S2': 3}
ATOMS = ['Mo', 'S1', 'S2']
NORB = sum(ORBS_PER_ATOM.values())  # 11

OFFSET = {}
o = 0
for at in ATOMS:
    OFFSET[at] = o
    o += ORBS_PER_ATOM[at]


# ---------------- Slater-Koster two-center integrals ----------------
SQ3 = np.sqrt(3.0)

def pp_block(l, m, n, Vppsig, Vpppi):
    E = np.zeros((3, 3))
    dc = np.array([l, m, n])
    for i in range(3):
        for j in range(3):
            if i == j:
                E[i, j] = dc[i]**2*Vppsig + (1-dc[i]**2)*Vpppi
            else:
                E[i, j] = dc[i]*dc[j]*(Vppsig-Vpppi)
    return E


def pd_block(l, m, n, Vpdsig, Vpdpi):
    E = np.zeros((3, 5))
    E[0, 0] = SQ3*l**2*m*Vpdsig + m*(1-2*l**2)*Vpdpi
    E[0, 1] = SQ3*l*m*n*Vpdsig - 2*l*m*n*Vpdpi
    E[0, 2] = SQ3*l**2*n*Vpdsig + n*(1-2*l**2)*Vpdpi
    E[0, 3] = 0.5*SQ3*l*(l**2-m**2)*Vpdsig + l*(1-l**2+m**2)*Vpdpi
    E[0, 4] = l*(n**2-0.5*(l**2+m**2))*Vpdsig - SQ3*l*n**2*Vpdpi
    E[1, 0] = SQ3*m**2*l*Vpdsig + l*(1-2*m**2)*Vpdpi
    E[1, 1] = SQ3*m**2*n*Vpdsig + n*(1-2*m**2)*Vpdpi
    E[1, 2] = SQ3*l*m*n*Vpdsig - 2*l*m*n*Vpdpi
    E[1, 3] = 0.5*SQ3*m*(l**2-m**2)*Vpdsig - m*(1+l**2-m**2)*Vpdpi
    E[1, 4] = m*(n**2-0.5*(l**2+m**2))*Vpdsig - SQ3*m*n**2*Vpdpi
    E[2, 0] = SQ3*l*m*n*Vpdsig - 2*l*m*n*Vpdpi
    E[2, 1] = SQ3*n**2*m*Vpdsig + m*(1-2*n**2)*Vpdpi
    E[2, 2] = SQ3*n**2*l*Vpdsig + l*(1-2*n**2)*Vpdpi
    E[2, 3] = 0.5*SQ3*n*(l**2-m**2)*Vpdsig - n*(l**2-m**2)*Vpdpi
    E[2, 4] = n*(n**2-0.5*(l**2+m**2))*Vpdsig + SQ3*n*(l**2+m**2)*Vpdpi
    return E


def dd_block(l, m, n, Vddsig, Vddpi, Vdddelta):
    E = np.zeros((5, 5))
    l2, m2, n2 = l*l, m*m, n*n
    E[0,0] = 3*l2*m2*Vddsig + (l2+m2-4*l2*m2)*Vddpi + (n2+l2*m2)*Vdddelta
    E[1,1] = 3*m2*n2*Vddsig + (m2+n2-4*m2*n2)*Vddpi + (l2+m2*n2)*Vdddelta
    E[2,2] = 3*n2*l2*Vddsig + (n2+l2-4*n2*l2)*Vddpi + (m2+n2*l2)*Vdddelta

    E[0,1] = 3*l*m2*n*Vddsig + l*n*(1-4*m2)*Vddpi + l*n*(m2-1)*Vdddelta
    E[1,0] = E[0,1]
    E[1,2] = 3*m*n2*l*Vddsig + m*l*(1-4*n2)*Vddpi + m*l*(n2-1)*Vdddelta
    E[2,1] = E[1,2]
    E[2,0] = 3*n*l2*m*Vddsig + n*m*(1-4*l2)*Vddpi + n*m*(l2-1)*Vdddelta
    E[0,2] = E[2,0]

    E[0,3] = 1.5*l*m*(l2-m2)*Vddsig + 2*l*m*(m2-l2)*Vddpi + 0.5*l*m*(l2-m2)*Vdddelta
    E[3,0] = E[0,3]
    E[1,3] = 1.5*m*n*(l2-m2)*Vddsig - m*n*(1+2*(l2-m2))*Vddpi + m*n*(1+0.5*(l2-m2))*Vdddelta
    E[3,1] = E[1,3]
    E[2,3] = 1.5*n*l*(l2-m2)*Vddsig + n*l*(1-2*(l2-m2))*Vddpi - n*l*(1-0.5*(l2-m2))*Vdddelta
    E[3,2] = E[2,3]

    E[0,4] = SQ3*l*m*(n2-0.5*(l2+m2))*Vddsig - 2*SQ3*l*m*n2*Vddpi + 0.5*SQ3*l*m*(1+n2)*Vdddelta
    E[4,0] = E[0,4]
    E[1,4] = SQ3*m*n*(n2-0.5*(l2+m2))*Vddsig + SQ3*m*n*(l2+m2-n2)*Vddpi - 0.5*SQ3*m*n*(l2+m2)*Vdddelta
    E[4,1] = E[1,4]
    E[2,4] = SQ3*n*l*(n2-0.5*(l2+m2))*Vddsig + SQ3*n*l*(l2+m2-n2)*Vddpi - 0.5*SQ3*n*l*(l2+m2)*Vdddelta
    E[4,2] = E[2,4]

    E[3,3] = 0.75*(l2-m2)**2*Vddsig + (l2+m2-(l2-m2)**2)*Vddpi + (n2+0.25*(l2-m2)**2)*Vdddelta
    E[3,4] = 0.5*SQ3*(l2-m2)*(n2-0.5*(l2+m2))*Vddsig + SQ3*n2*(m2-l2)*Vddpi + 0.25*SQ3*(1+n2)*(l2-m2)*Vdddelta
    E[4,3] = E[3,4]
    E[4,4] = (n2-0.5*(l2+m2))**2*Vddsig + 3*n2*(l2+m2)*Vddpi + 0.75*(l2+m2)**2*Vdddelta
    return E


PARAM_NAMES = [
    'e_Mo_a1', 'e_Mo_e', 'e_Mo_e2',
    'e_S_p_para', 'e_S_p_perp',
    'Vpd_sigma', 'Vpd_pi',
    'Vdd_sigma', 'Vdd_pi', 'Vdd_delta',
    'Vpp_sigma', 'Vpp_pi',
]
N_PARAMS = len(PARAM_NAMES)


def onsite_matrix(p): #contruye la parte digonal del modelo
    H0 = np.zeros((NORB, NORB))
    o = OFFSET['Mo']
    H0[o+0, o+0] = p['e_Mo_e']
    H0[o+1, o+1] = p['e_Mo_e2']
    H0[o+2, o+2] = p['e_Mo_e2']
    H0[o+3, o+3] = p['e_Mo_e']
    H0[o+4, o+4] = p['e_Mo_a1']
    for s in ('S1', 'S2'):
        o = OFFSET[s]
        H0[o+0, o+0] = p['e_S_p_para']
        H0[o+1, o+1] = p['e_S_p_para']
        H0[o+2, o+2] = p['e_S_p_perp']
    return H0


def build_neighbor_shells(cutoff=3.3):#contruye la lista de primeros vecinos
    shells = []
    rng = range(-2, 3)
    for ai in ATOMS:
        for aj in ATOMS:
            for n1 in rng:
                for n2 in rng:
                    Rcell = n1*A1 + n2*A2
                    pos_j = CART[aj] + Rcell
                    d = pos_j - CART[ai]
                    dist = np.linalg.norm(d)
                    if 1e-3 < dist < cutoff:
                        shells.append({'ai': ai, 'aj': aj, 'n1': n1, 'n2': n2, 'd': d, 'dist': dist})
    return shells


SHELLS = build_neighbor_shells(cutoff=3.3)


def hop_block(ai, aj, l, m, n, p):
    is_d_i = (ai == 'Mo')
    is_d_j = (aj == 'Mo')
    if is_d_i and is_d_j:
        return dd_block(l, m, n, p['Vdd_sigma'], p['Vdd_pi'], p['Vdd_delta'])
    if is_d_i and not is_d_j:
        # <d|H|p> = <p|H|d>^T evaluated with the p->d direction (pd_block is odd under
        # bond reversal, so we must flip (l,m,n) which here is given as d->p).
        return pd_block(-l, -m, -n, p['Vpd_sigma'], p['Vpd_pi']).T
    if (not is_d_i) and is_d_j:
        return pd_block(l, m, n, p['Vpd_sigma'], p['Vpd_pi'])
    return pp_block(l, m, n, p['Vpp_sigma'], p['Vpp_pi'])


def build_Hk(kcart, params_vec):
    p = dict(zip(PARAM_NAMES, params_vec))
    H = onsite_matrix(p).astype(complex)
    for sh in SHELLS: #vecinos
        ai, aj, d = sh['ai'], sh['aj'], sh['d'] 
        Rc = sh['n1']*A1 + sh['n2']*A2
        dist = sh['dist']
        l, m, n = d / dist
        block = hop_block(ai, aj, l, m, n, p) #término de la tabla de slater
        oi, oj = OFFSET[ai], OFFSET[aj]
        ni, nj = ORBS_PER_ATOM[ai], ORBS_PER_ATOM[aj]
        phase = np.exp(1j*np.dot(kcart, Rc))
        H[oi:oi+ni, oj:oj+nj] += block*phase
    return H


def bands_at_k(kcart_array, params_vec):
    nk = kcart_array.shape[0]
    out = np.zeros((nk, NORB))
    for ik in range(nk):
        H = build_Hk(kcart_array[ik], params_vec)
        H = 0.5*(H + H.conj().T)
        ev = np.linalg.eigvalsh(H)
        out[ik] = ev
    return out
