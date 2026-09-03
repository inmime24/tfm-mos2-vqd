"""
Calcula, mediante VQD, los 11 autovalores mas bajos en ~40 puntos k del
camino M-Gamma-K-M, encadenando warm-start entre puntos vecinos.

Guardado incremental tras CADA (punto, nivel) -- si el trabajo se corta
por el limite de tiempo de la cola, al relanzarlo se reanuda exactamente
donde se quedo, sin repetir nada ya hecho.

Por cada (punto, nivel) se registra: energia VQD, energia exacta, error,
solapamiento maximo con los estados ya encontrados en ese mismo punto,
tiempo de ejecucion, numero de evaluaciones, si convergió segun scipy, y
el origen de la semilla usada (warm-start del punto anterior, o arranque
en frio con su propia semilla).

El primer punto de la lista usa arranque en frio multi-semilla.

Uso:
    python3 calcular_bandas_vqd.py [--pkl datos/hamiltonianos_40.pkl]
                                    [--n-estados 11] [--reps 2]
                                    [--metodo L-BFGS-B] [--maxiter 5000]
                                    [--beta 10] [--salida ruta.json]
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

from ansatz_particula import crear_ansatz

PKL_HAMILTONIANOS = 'datos/hamiltonianos_40.pkl'
N_ESTADOS = 11
SEMILLA_ESTADO0_FRIO = 2         
N_SEMILLAS_FRIO = 2              
REPS = 2
METODO = 'L-BFGS-B'
MAXITER = 5000
BETA = 10.0
SALIDA_JSON = 'resultados/bandas_vqd.json'


# funciones necesarias

def construir_circuito(ansatz, n_qubits, parametros):
    estado_inicial = QuantumCircuit(n_qubits)
    estado_inicial.x(0)
    return estado_inicial.compose(ansatz.assign_parameters(parametros))


def opciones_optimizador(metodo, maxiter):
    if metodo.upper() == 'COBYLA':
        return {'maxiter': maxiter}
    return {'maxiter': maxiter, 'maxfun': maxiter}


def optimizar_desde(funcion_costo, p0, metodo, maxiter):
    opciones = opciones_optimizador(metodo, maxiter)
    t_inicio = time.time()
    resultado_opt = minimize(funcion_costo, p0, method=metodo, options=opciones)
    tiempo_s = time.time() - t_inicio
    return resultado_opt, tiempo_s


def optimizar_con_semilla(funcion_costo, n_params, semilla, metodo, maxiter):
    rng = np.random.default_rng(semilla)
    p0 = rng.uniform(0, 2 * np.pi, n_params)
    return optimizar_desde(funcion_costo, p0, metodo, maxiter)


def guardar_progreso(salida_json, resultados_por_punto):
    os.makedirs(os.path.dirname(salida_json), exist_ok=True)
    tmp = salida_json + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({'resultados_por_punto': resultados_por_punto}, f, indent=2)
    os.replace(tmp, salida_json)  # escritura atomica: nunca deja el .json a medias



# Resolver un unico nivel en un punto dado

def resolver_nivel_frio(pauli_op, estados_encontrados, ansatz, N, nivel,
                         semillas, metodo, maxiter, beta):
    """Arranque en frio multi-semilla (solo para el primer punto)."""

    def funcion_costo(parametros):
        circuito = construir_circuito(ansatz, N, parametros)
        psi = Statevector(circuito)
        energia = psi.expectation_value(pauli_op).real
        penal = sum(beta * abs(p.inner(psi)) ** 2 for p in estados_encontrados)
        return energia + penal

    candidatos = []
    for semilla in semillas:
        resultado_opt, tiempo_s = optimizar_con_semilla(
            funcion_costo, ansatz.num_parameters, semilla, metodo, maxiter)
        psi = Statevector(construir_circuito(ansatz, N, resultado_opt.x))
        energia_real = float(psi.expectation_value(pauli_op).real)
        candidatos.append({
            'theta': resultado_opt.x, 'energia': energia_real, 'psi': psi,
            'tiempo_s': tiempo_s, 'n_evaluaciones': int(resultado_opt.nfev),
            'convergio': bool(resultado_opt.success),
            'origen': f'frio:semilla={semilla}',
        })

    return min(candidatos, key=lambda c: c['energia'])


def resolver_nivel_warm(pauli_op, estados_encontrados, ansatz, N, theta_previo,
                         metodo, maxiter, beta, etiqueta_punto_previo):
    """Warm-start: un unico intento, partiendo de la solucion del punto anterior."""

    def funcion_costo(parametros):
        circuito = construir_circuito(ansatz, N, parametros)
        psi = Statevector(circuito)
        energia = psi.expectation_value(pauli_op).real
        penal = sum(beta * abs(p.inner(psi)) ** 2 for p in estados_encontrados)
        return energia + penal

    resultado_opt, tiempo_s = optimizar_desde(funcion_costo, theta_previo, metodo, maxiter)
    psi = Statevector(construir_circuito(ansatz, N, resultado_opt.x))
    energia_real = float(psi.expectation_value(pauli_op).real)
    return {
        'theta': resultado_opt.x, 'energia': energia_real, 'psi': psi,
        'tiempo_s': tiempo_s, 'n_evaluaciones': int(resultado_opt.nfev),
        'convergio': bool(resultado_opt.success),
        'origen': f'warm_start:{etiqueta_punto_previo}',
    }



# Función principal


def calcular_bandas_vqd(pkl_hamiltonianos=PKL_HAMILTONIANOS, n_estados=N_ESTADOS,
                         reps=REPS, metodo=METODO, maxiter=MAXITER, beta=BETA,
                         salida_json=SALIDA_JSON):

    with open(pkl_hamiltonianos, 'rb') as f:
        puntos = pickle.load(f)  # lista de {etiqueta, kcart, pauli_op, H}, en orden del camino

    N = puntos[0]['H'].shape[0]
    ansatz = crear_ansatz(N, reps=reps)
    n_params = ansatz.num_parameters

    print(f"{len(puntos)} puntos k  |  N={N} qubits  |  reps={reps}  |  "
          f"{n_params} parametros  |  metodo={metodo}  |  maxiter={maxiter}  |  beta={beta}\n")

    # ---------------- cargar progreso previo, si existe (reanudacion) ----------------
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
        pauli_op, H, etiqueta = punto['pauli_op'], punto['H'], punto['etiqueta']
        ev_directo = np.linalg.eigvalsh(0.5 * (H + H.conj().T))

        # theta del punto anterior, por nivel (None si este es el primer punto de todos)
        if idx_punto > 0:
            theta_previo_por_nivel = [np.array(r['theta']) for r in resultados_por_punto[-1]['niveles']]
            etiqueta_previa = resultados_por_punto[-1]['etiqueta']
        else:
            theta_previo_por_nivel = None
            etiqueta_previa = None

        # si estamos reanudando a medias de este punto, recuperar lo ya hecho
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

            if theta_previo_por_nivel is not None:
                mejor = resolver_nivel_warm(pauli_op, estados_encontrados, ansatz, N,
                                             theta_previo_por_nivel[nivel], metodo, maxiter,
                                             beta, etiqueta_previa)
            else:
                semillas = [SEMILLA_ESTADO0_FRIO] if nivel == 0 else list(range(1, N_SEMILLAS_FRIO + 1))
                mejor = resolver_nivel_frio(pauli_op, estados_encontrados, ansatz, N, nivel,
                                             semillas, metodo, maxiter, beta)

            estados_encontrados.append(mejor['psi'])

            energia_exacta = float(ev_directo[nivel])
            error = abs(mejor['energia'] - energia_exacta)
            solapamientos_previos = [
                float(abs(estados_encontrados[j].inner(mejor['psi'])) ** 2)
                for j in range(len(estados_encontrados) - 1)
            ]
            solap_max = max(solapamientos_previos) if solapamientos_previos else 0.0

            print(f"  nivel {nivel}: energia={mejor['energia']:.6f} "
                  f"(exacto={energia_exacta:.6f}) error={error:.2e} "
                  f"solap_max={solap_max:.2e} tiempo={mejor['tiempo_s']:.1f}s "
                  f"evals={mejor['n_evaluaciones']} convergio={mejor['convergio']} "
                  f"origen={mejor['origen']}")

            niveles_resultado.append({
                'nivel': nivel,
                'energia_vqd': mejor['energia'],
                'energia_exacta': energia_exacta,
                'error': error,
                'solapamiento_maximo': solap_max,
                'theta': mejor['theta'].tolist(),
                'tiempo_s': mejor['tiempo_s'],
                'n_evaluaciones': mejor['n_evaluaciones'],
                'convergio': mejor['convergio'],
                'origen_semilla': mejor['origen'],
            })

            # ---- GUARDADO INCREMENTAL: tras cada nivel, no solo al final del punto ----
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
    p.add_argument('--reps', type=int, default=REPS)
    p.add_argument('--metodo', default=METODO)
    p.add_argument('--maxiter', type=int, default=MAXITER)
    p.add_argument('--beta', type=float, default=BETA)
    p.add_argument('--salida', default=SALIDA_JSON)
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    calcular_bandas_vqd(pkl_hamiltonianos=args.pkl, n_estados=args.n_estados,
                         reps=args.reps, metodo=args.metodo, maxiter=args.maxiter,
                         beta=args.beta, salida_json=args.salida)
