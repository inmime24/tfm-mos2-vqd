"""
calcular_bandas_vqd_ruido.py

Version CON RUIDO (muestreo via Qulacs) de calcular_bandas_vqd.py. Calcula,
mediante VQD, los autovalores de un conjunto de puntos k (definido por el
.pkl de entrada -- ese archivo decide "que puntos"), usando una funcion de
coste con muestreo real en vez de valor esperado exacto. La penalizacion de
solapamiento entre estados se mantiene EXACTA (decision ya tomada: solo la
energia se mide con ruido).

scipy.optimize.minimize NO respeta maxfun de forma fiable con L-BFGS-B
cuando el gradiente se estima por diferencias finitas (comprobado
empiricamente: con maxfun=30 se llegaron a hacer mas de 390 llamadas sin
detenerse). Por eso se usa un mecanismo propio de parada forzada por
excepcion, que garantiza el presupuesto de evaluaciones exacto para ambos
optimizadores (L-BFGS-B y SPSA), permitiendo una comparacion justa.

Guardado incremental identico a calcular_bandas_vqd.py: permite reanudar
si el trabajo se corta por el limite de tiempo de la cola.

Uso:
    python3 calcular_bandas_vqd_ruido.py --pkl datos/hamiltonianos_16puntos.pkl \
        --n-estados 11 --shots 4096 --metodo SPSA --presupuesto 3000 \
        --beta 10 --salida resultados/bandas_ruido_warmstart.json
"""
import os
import time
import json
import pickle
import argparse

import numpy as np
from scipy.optimize import minimize
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qulacs import QuantumState
from qulacs.gate import H, Sdag
from qiskit_algorithms.optimizers import SPSA

from ansatz_particula import crear_ansatz

PKL_HAMILTONIANOS = 'datos/hamiltonianos_precalculados.pkl'
N_ESTADOS = 11
SHOTS = 4096
METODO = 'SPSA'   # pendiente de confirmar/ajustar con la prueba comparativa
PRESUPUESTO_EVALS = 3000
BETA = 10.0
SALIDA_JSON = 'resultados/bandas_vqd_ruido.json'


# ============================================================
# Utilidades de circuito (identicas a calcular_bandas_vqd.py)
# ============================================================

def construir_circuito(ansatz, n_qubits, parametros):
    estado_inicial = QuantumCircuit(n_qubits)
    estado_inicial.x(0)
    return estado_inicial.compose(ansatz.assign_parameters(parametros))


def guardar_progreso(salida_json, resultados_por_punto):
    os.makedirs(os.path.dirname(salida_json), exist_ok=True)
    tmp = salida_json + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({'resultados_por_punto': resultados_por_punto}, f, indent=2)
    os.replace(tmp, salida_json)


# ============================================================
# Muestreo via Qulacs (verificado con un estado de Bell)
# ============================================================

def medir_termino_qulacs(qs_base, soporte, shots):
    qs = qs_base.copy()
    qubits_medidos = []
    for q, p in soporte.items():
        if p == 'X':
            H(q).update_quantum_state(qs)
        elif p == 'Y':
            Sdag(q).update_quantum_state(qs)
            H(q).update_quantum_state(qs)
        qubits_medidos.append(q)
    resultados = qs.sampling(shots)
    valor = 0.0
    for r in resultados:
        paridad = 1
        for qi in qubits_medidos:
            if (r >> qi) & 1:
                paridad *= -1
        valor += paridad
    return valor / shots


def evaluar_pauli_con_muestreo(pauli_op, psi_qiskit, shots):
    N = psi_qiskit.num_qubits
    qs_base = QuantumState(N)
    qs_base.load(psi_qiskit.data)

    labels = pauli_op.paulis.to_labels()
    coeffs = pauli_op.coeffs.real

    valor = 0.0
    for label, coef in zip(labels, coeffs):
        soporte = {i: ch for i, ch in enumerate(reversed(label)) if ch != 'I'}
        if not soporte:
            valor += coef.real
            continue
        valor += coef.real * medir_termino_qulacs(qs_base, soporte, shots)
    return valor


# ============================================================
# Optimizacion con presupuesto de evaluaciones garantizado
# ============================================================

class PresupuestoAgotado(Exception):
    pass


def crear_funcion_costo_con_limite(funcion_costo_base, maxfun):
    estado = {'contador': 0, 'mejor_valor': np.inf, 'mejor_theta': None}

    def funcion_costo(theta):
        if estado['contador'] >= maxfun:
            raise PresupuestoAgotado()
        estado['contador'] += 1
        valor = funcion_costo_base(theta)
        if valor < estado['mejor_valor']:
            estado['mejor_valor'] = valor
            estado['mejor_theta'] = np.array(theta)
        return valor

    return funcion_costo, estado


def optimizar_con_presupuesto(funcion_costo_base, p0, metodo, maxfun):
    funcion_costo, estado = crear_funcion_costo_con_limite(funcion_costo_base, maxfun)

    t0 = time.time()
    try:
        if metodo == 'L-BFGS-B':
            minimize(funcion_costo, p0, method='L-BFGS-B',
                     options={'maxfun': maxfun, 'maxiter': maxfun})
        elif metodo == 'SPSA':
            SPSA(maxiter=maxfun // 2).minimize(funcion_costo, p0)
        else:
            raise ValueError(f"Metodo desconocido: {metodo}")
    except PresupuestoAgotado:
        pass
    tiempo = time.time() - t0

    return estado['mejor_theta'], estado['contador'], tiempo


# ============================================================
# Orquestador principal (mismo esqueleto que calcular_bandas_vqd.py)
# ============================================================

def calcular_bandas_vqd_ruido(pkl_hamiltonianos=PKL_HAMILTONIANOS, n_estados=N_ESTADOS,
                               shots=SHOTS, metodo=METODO, presupuesto_evals=PRESUPUESTO_EVALS,
                               beta=BETA, salida_json=SALIDA_JSON):

    with open(pkl_hamiltonianos, 'rb') as f:
        puntos = pickle.load(f)

    N = puntos[0]['H'].shape[0]
    ansatz = crear_ansatz(N, reps=2, entanglement='full')

    print(f"{len(puntos)} puntos k  |  N={N} qubits  |  shots={shots}  |  "
          f"metodo={metodo}  |  presupuesto={presupuesto_evals} evals/nivel  |  beta={beta}\n")

    resultados_por_punto = []
    punto_parcial = None
    if os.path.exists(salida_json):
        with open(salida_json) as f:
            resultados_por_punto = json.load(f)['resultados_por_punto']
        if resultados_por_punto and len(resultados_por_punto[-1]['niveles']) < n_estados:
            punto_parcial = resultados_por_punto.pop()
            print(f"Reanudando: punto {punto_parcial['indice_punto']} "
                  f"({punto_parcial['etiqueta']}) tenia {len(punto_parcial['niveles'])}/"
                  f"{n_estados} niveles -- se completa desde ahi.")
        else:
            print(f"Reanudando: {len(resultados_por_punto)} puntos ya completos.")

    indice_inicio = len(resultados_por_punto) if punto_parcial is None else punto_parcial['indice_punto']

    for idx_punto in range(indice_inicio, len(puntos)):
        punto = puntos[idx_punto]
        pauli_op, H_matrix, etiqueta = punto['pauli_op'], punto['H'], punto['etiqueta']
        ev_directo = np.linalg.eigvalsh(0.5 * (H_matrix + H_matrix.conj().T))

        if idx_punto > 0:
            theta_previo_por_nivel = [np.array(r['theta']) for r in resultados_por_punto[-1]['niveles']]
            etiqueta_previa = resultados_por_punto[-1]['etiqueta']
        else:
            theta_previo_por_nivel = None
            etiqueta_previa = None

        if punto_parcial is not None and idx_punto == punto_parcial['indice_punto']:
            niveles_resultado = punto_parcial['niveles']
            estados_encontrados = [
                Statevector(construir_circuito(ansatz, N, np.array(r['theta'])))
                for r in niveles_resultado
            ]
            t_inicio_punto = time.time() - punto_parcial.get('tiempo_total_punto_s', 0.0)
            nivel_inicio = len(niveles_resultado)
        else:
            niveles_resultado = []
            estados_encontrados = []
            t_inicio_punto = time.time()
            nivel_inicio = 0

        print(f"--- Punto {idx_punto} ({etiqueta}) ---")

        for nivel in range(nivel_inicio, n_estados):

            def funcion_costo_ruidosa(theta, _estados=estados_encontrados):
                psi = Statevector(construir_circuito(ansatz, N, theta))
                energia_ruidosa = evaluar_pauli_con_muestreo(pauli_op, psi, shots)
                # la penalizacion de solapamiento se mantiene EXACTA (no ruidosa)
                penal = sum(beta * abs(prev.inner(psi)) ** 2 for prev in _estados)
                return energia_ruidosa + penal

            if theta_previo_por_nivel is not None:
                p0 = theta_previo_por_nivel[nivel]
                origen = f'warm_start:{etiqueta_previa}'
            else:
                rng = np.random.default_rng(1)
                p0 = rng.uniform(0, 2 * np.pi, ansatz.num_parameters)
                origen = 'frio:semilla=1'

            theta_final, n_evals, tiempo_s = optimizar_con_presupuesto(
                funcion_costo_ruidosa, p0, metodo, presupuesto_evals)

            psi_final = Statevector(construir_circuito(ansatz, N, theta_final))
            energia_exacta_final = float(psi_final.expectation_value(pauli_op).real)
            energia_ruidosa_final = evaluar_pauli_con_muestreo(pauli_op, psi_final, shots)

            energia_exacta_ref = float(ev_directo[nivel])
            error_vs_exacto = abs(energia_exacta_final - energia_exacta_ref)
            solapamientos_previos = [
                float(abs(estados_encontrados[j].inner(psi_final)) ** 2)
                for j in range(len(estados_encontrados))
            ]
            solap_max = max(solapamientos_previos) if solapamientos_previos else 0.0

            estados_encontrados.append(psi_final)

            print(f"  nivel {nivel}: energia_exacta={energia_exacta_final:.6f} "
                  f"(ref={energia_exacta_ref:.6f}) error={error_vs_exacto:.4e} "
                  f"energia_ruidosa={energia_ruidosa_final:.6f} "
                  f"solap_max={solap_max:.2e} tiempo={tiempo_s:.1f}s "
                  f"evals={n_evals} origen={origen}")

            niveles_resultado.append({
                'nivel': nivel,
                'energia_exacta_final': energia_exacta_final,
                'energia_ruidosa_final': energia_ruidosa_final,
                'energia_exacta_referencia': energia_exacta_ref,
                'error_vs_exacto': error_vs_exacto,
                'solapamiento_maximo': solap_max,
                'theta': theta_final.tolist(),
                'tiempo_s': tiempo_s,
                'n_evaluaciones': n_evals,
                'origen_semilla': origen,
                'shots': shots,
                'metodo': metodo,
            })

            tiempo_total_punto_parcial = time.time() - t_inicio_punto
            guardar_progreso(salida_json, resultados_por_punto + [{
                'indice_punto': idx_punto,
                'etiqueta': etiqueta,
                'kcart': list(punto['kcart']),
                'tiempo_total_punto_s': tiempo_total_punto_parcial,
                'niveles': niveles_resultado,
            }])

        tiempo_total_punto = time.time() - t_inicio_punto
        resultados_por_punto.append({
            'indice_punto': idx_punto,
            'etiqueta': etiqueta,
            'kcart': list(punto['kcart']),
            'tiempo_total_punto_s': tiempo_total_punto,
            'niveles': niveles_resultado,
        })
        guardar_progreso(salida_json, resultados_por_punto)
        punto_parcial = None
        print(f"  Punto {idx_punto} completo en {tiempo_total_punto:.1f}s\n")

    print(f"Barrido completo: {len(resultados_por_punto)} puntos. Guardado en {salida_json}")
    return resultados_por_punto


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--pkl', default=PKL_HAMILTONIANOS)
    p.add_argument('--n-estados', type=int, default=N_ESTADOS)
    p.add_argument('--shots', type=int, default=SHOTS)
    p.add_argument('--metodo', default=METODO, choices=['L-BFGS-B', 'SPSA'])
    p.add_argument('--presupuesto', type=int, default=PRESUPUESTO_EVALS)
    p.add_argument('--beta', type=float, default=BETA)
    p.add_argument('--salida', default=SALIDA_JSON)
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    calcular_bandas_vqd_ruido(pkl_hamiltonianos=args.pkl, n_estados=args.n_estados,
                               shots=args.shots, metodo=args.metodo,
                               presupuesto_evals=args.presupuesto, beta=args.beta,
                               salida_json=args.salida)
