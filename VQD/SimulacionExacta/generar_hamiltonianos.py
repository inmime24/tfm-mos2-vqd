"""
Dado el modelo tight-binding ya ajustado (fit_result.npz) y el camino k
(kpath.npz), genera los Hamiltonianos H(k) en los puntos solicitados y los
convierte a SparsePauliOp (Jordan-Wigner), guardando todo en un unico
archivo.

Modos:
    "puntos_especiales": argumento = lista de nombres, p.ej. ["K"] o ["K","M","Gamma"]
    "equiespaciado":      argumento = N (numero de puntos)

Uso:
    generar_hamiltonianos(
        fit_result_path='datos/fit_result.npz',
        kpath_path='datos/kpath.npz',
        modo='puntos_especiales',
        argumento=['K'],
        salida='resultados/hamiltonianos_precalculados.pkl',
    )
"""
import numpy as np
import pickle
from qiskit.quantum_info import SparsePauliOp
import sk_model as m

# Funciones auxiliares

def coordenadas_alta_simetria(nombre):
    """Devuelve el vector k cartesiano de un punto de alta simetria
    (M, Gamma, K), usando la misma geometria que sk_model.py."""
    A1, A2, A3 = m.A1, m.A2, m.A3
    vol = np.dot(A1, np.cross(A2, A3))
    B1 = 2 * np.pi * np.cross(A2, A3) / vol
    B2 = 2 * np.pi * np.cross(A3, A1) / vol
    B3 = 2 * np.pi * np.cross(A1, A2) / vol

    tabla_fraccional = {
        'M': np.array([0.5, 0.0, 0.0]),
        'Gamma': np.array([0.0, 0.0, 0.0]),
        'K': np.array([1 / 3, 1 / 3, 0.0]),
    }
    if nombre not in tabla_fraccional:
        raise ValueError(f"Punto de alta simetria desconocido: '{nombre}'. "
                          f"Opciones: {list(tabla_fraccional.keys())}")

    kfrac = tabla_fraccional[nombre]
    kcart = kfrac[0] * B1 + kfrac[1] * B2 + kfrac[2] * B3
    return kcart


def seleccionar_indices_equiespaciados(N, total_puntos):
    """Selecciona N indices igualmente repartidos entre 0 y total_puntos-1."""
    if N >= total_puntos:
        print(f"Aviso: se pidieron {N} puntos pero el camino solo tiene "
              f"{total_puntos}; se devuelven todos.")
        return list(range(total_puntos))

    if N == 1:
        return [0]

    indices = []
    for i in range(N):
        idx = round(i * (total_puntos - 1) / (N - 1))
        indices.append(idx)

    vistos = set()
    indices_unicos = []
    for idx in indices:
        if idx not in vistos:
            vistos.add(idx)
            indices_unicos.append(idx)

    if len(indices_unicos) < N:
        print(f"Aviso: se pidieron {N} puntos equiespaciados, pero tras "
              f"quitar duplicados solo quedaron {len(indices_unicos)}.")

    return indices_unicos


def construir_sparse_pauli_op(H, tol=1e-12):
    """Descompone la matriz H (NxN, hermitica) en cadenas de Pauli via
    Jordan-Wigner. Identidades usadas:
        n_i               = (I - Z_i)/2
        Re(H_ij), i<j:     coef 1/2 en X_i X_j y en Y_i Y_j (con Z-string)
        Im(H_ij), i<j:     coef 1/2 en Y_i X_j y -1/2 en X_i Y_j (con Z-string)
    """
    N = H.shape[0]
    terminos = {}

    def anadir_termino(soporte, coef):
        if abs(coef) < tol:
            return
        clave = tuple(sorted(soporte.items()))
        terminos[clave] = terminos.get(clave, 0.0) + coef

    id_coef = 0.0
    for i in range(N):
        hii = H[i, i].real
        id_coef += hii / 2
        anadir_termino({i: 'Z'}, -hii / 2)
    anadir_termino({}, id_coef)

    for i in range(N):
        for j in range(i + 1, N):
            hij = H[i, j]
            a, b = hij.real, hij.imag
            zstring = {q: 'Z' for q in range(i + 1, j)}
            if abs(a) > tol:
                anadir_termino({**zstring, i: 'X', j: 'X'}, a / 2)
                anadir_termino({**zstring, i: 'Y', j: 'Y'}, a / 2)
            if abs(b) > tol:
                anadir_termino({**zstring, i: 'Y', j: 'X'}, b / 2)
                anadir_termino({**zstring, i: 'X', j: 'Y'}, -b / 2)

    pauli_strings = []
    coeffs = []
    for soporte_items, coef in terminos.items():
        soporte = dict(soporte_items)
        # SparsePauliOp: caracter mas a la izquierda = qubit de indice mas alto
        cadena = ''.join(soporte.get(q, 'I') for q in range(N - 1, -1, -1))
        pauli_strings.append(cadena)
        coeffs.append(coef.real if hasattr(coef, 'real') else coef)

    return SparsePauliOp(pauli_strings, coeffs=np.array(coeffs, dtype=complex))


# Funcion principal

def generar_hamiltonianos(fit_result_path, kpath_path, modo, argumento, salida):
    r = np.load(fit_result_path)
    params = r['params']

    d = np.load(kpath_path)
    kcart_completo = d['kcart']
    total_puntos = kcart_completo.shape[0]

    if modo == 'puntos_especiales':
        nombres = argumento
        lista_k = [(nombre, coordenadas_alta_simetria(nombre)) for nombre in nombres]

    elif modo == 'equiespaciado':
        N = argumento
        indices = seleccionar_indices_equiespaciados(N, total_puntos)
        lista_k = [(idx, kcart_completo[idx]) for idx in indices]

    else:
        raise ValueError(f"Modo desconocido: '{modo}'. Usa 'puntos_especiales' o 'equiespaciado'.")

    resultados = []
    for etiqueta, kcart in lista_k:
        H = m.build_Hk(kcart, params)
        H = 0.5 * (H + H.conj().T)
        pauli_op = construir_sparse_pauli_op(H)
        resultados.append({'etiqueta': etiqueta, 'kcart': kcart, 'pauli_op': pauli_op, 'H': H})
        print(f'{etiqueta}: {len(pauli_op.paulis)} terminos de Pauli')

    with open(salida, 'wb') as f:
        pickle.dump(resultados, f)

    print(f'\nGuardado {len(resultados)} Hamiltonianos en {salida}')
    return resultados


if __name__ == '__main__':
    generar_hamiltonianos(
        fit_result_path='datos/fit_result.npz',
        kpath_path='datos/kpath.npz',
        modo='equiespaciado',
        argumento=40,
        salida='resultados/hamiltonianos_40.pkl',
    )