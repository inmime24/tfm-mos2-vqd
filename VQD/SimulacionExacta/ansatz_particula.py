"""
Ansatz variacional restringido al sector de una particula del espacio de
Fock (circuito tipo excitation-preserving, con puertas que conservan el
numero de excitaciones). Incluye una comprobacion numerica de que el
numero de particulas se mantiene en 1 para distintos valores de sus
parametros.
"""

import numpy as np
from qiskit.circuit.library import excitation_preserving
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector


def crear_ansatz(n_qubits, reps=2, entanglement='full'):
    return excitation_preserving(n_qubits, reps=reps, entanglement=entanglement)


def operador_numero_total(n_qubits):
    """N_total = sum_i (I - Z_i)/2, como SparsePauliOp."""
    etiquetas = ['I' * n_qubits]
    coefs = [n_qubits / 2]
    for i in range(n_qubits):
        soporte = ['I'] * n_qubits
        soporte[n_qubits - 1 - i] = 'Z'  # qubit i -> posicion desde la derecha (convenio Qiskit)
        etiquetas.append(''.join(soporte))
        coefs.append(-0.5)
    return SparsePauliOp(etiquetas, coeffs=np.array(coefs, dtype=complex))


def verificar_conservacion_particula(circuito, n_qubits, n_pruebas=20,
                                      semilla=42, tolerancia=1e-8):
    """Prepara el estado de 1 particula |10...0>, aplica el ansatz con
    varios juegos de parametros aleatorios, y comprueba que <N_total> = 1
    en todos los casos."""
    N_total = operador_numero_total(n_qubits)
    rng = np.random.default_rng(semilla)
    n_params = circuito.num_parameters

    print(f"{'Prueba':>8}  {'<N_total>':>14}  Resultado")
    error_max = 0.0
    for prueba in range(n_pruebas):
        valores = rng.uniform(0, 2 * np.pi, n_params)

        estado_inicial = QuantumCircuit(n_qubits)
        estado_inicial.x(0)  # qubit 0 ocupado: estado de 1 particula

        circuito_completo = estado_inicial.compose(
            circuito.assign_parameters(valores)
        )

        psi = Statevector(circuito_completo)
        valor_N = psi.expectation_value(N_total).real

        error = abs(valor_N - 1.0)
        error_max = max(error_max, error)
        ok = 'OK' if error < tolerancia else 'REVISAR'
        print(f"{prueba:>8}  {valor_N:>14.10f}  {ok}")

    print(f"\nError maximo respecto a N=1: {error_max:.2e} "
          f"(tolerancia {tolerancia:.0e})")
    return error_max


if __name__ == '__main__':
    N_QUBITS = 11
    REPS = 2
    ENTANGLEMENT = 'full'

    ansatz = crear_ansatz(N_QUBITS, reps=REPS, entanglement=ENTANGLEMENT)
    print(f"Ansatz: {N_QUBITS} qubits, reps={REPS}, entanglement='{ENTANGLEMENT}', "
          f"{ansatz.num_parameters} parametros\n")

    verificar_conservacion_particula(ansatz, N_QUBITS, n_pruebas=20)
