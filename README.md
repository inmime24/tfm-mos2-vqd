# Simulación cuántica de la estructura electrónica del MoS₂ mediante VQD

Código desarrollado para el Trabajo de Fin de Máster *"Simulación cuántica de la
estructura electrónica del MoS₂ mediante VQD"*, presentado en la Facultade de
Informática de la Universidade da Coruña dentro del Máster en Quantum Information
Science and Technologies (curso 2025).

Este repositorio recoge el proceso completo seguido en el trabajo: el cálculo de
referencia mediante DFT, el ajuste del modelo tight-binding sobre dicha referencia,
y la resolución del Hamiltoniano resultante mediante el algoritmo cuántico
variacional VQD, tanto en simulación exacta como con muestreo estadístico.

## Construcción general del repositorio

El repositorio se ha organizado en tres carpetas principales, siguiendo el mismo
orden en el que se presentan los métodos en el capítulo 4 de la memoria: primero
el cálculo de referencia (`DFT/`), después el ajuste del modelo tight-binding a
dicho cálculo (`AJUSTE/`), y por último la resolución cuántica del Hamiltoniano
ya ajustado (`VQD/`). Cada carpeta es, en la práctica, autocontenida: toma como
entrada la salida de la etapa anterior y produce los resultados que se discuten
en el capítulo 5.

```
DFT/        Cálculo de referencia mediante SIESTA
AJUSTE/     Ajuste del modelo tight-binding a las bandas de DFT
VQD/        Resolución del Hamiltoniano mediante VQD
```

Dentro de cada carpeta se ha mantenido, en la medida de lo posible, la misma
subdivisión: un subdirectorio `datos/` con las entradas necesarias, un
subdirectorio `resultados/` con las salidas generadas, y los scripts de Python
correspondientes al mismo nivel. Esta convención se sigue de forma estricta en
`AJUSTE/` y en `VQD/SimulacionExacta/`; en `VQD/Muestreo/` los ficheros `.pkl`
que definen los puntos k de cada experimento se han dejado directamente en la
raíz de la carpeta en lugar de en `datos/`, ya que así es como se referencian
desde los scripts de envío al clúster (`slurm/`).

## `DFT/`

Contiene la configuración y los resultados de los cálculos DFT realizados con
SIESTA (sección 4.1 de la memoria), organizados según el tipo de cálculo:

- `OPT/`: optimización de la base de orbitales, del radio de corte, del
  `MeshCutoff` y de la malla de puntos k, así como la relajación estructural
  (`RESULTADOS_OPTIMIZACIÓN/`, `RESULTADOS_RELAJACIÓN/`) y el barrido del
  parámetro de red (`OPT_CTERED/`, una subcarpeta por valor de `a`).
- `SIN_SPIN/`: estructura de bandas sin polarización de espín, empleada como
  referencia para el ajuste (corresponde a la Figura 5.1 de la memoria).
- `SPIN/` y `SPIN_ORBITA/`: cálculos de comprobación con polarización de espín
  y con acoplamiento espín-órbita, respectivamente (sección 5.1).

Cada subcarpeta de cálculo conserva el `input.fdf` empleado, el `output.out`
correspondiente y los ficheros de resultados de SIESTA con valor propio
(`mos2.bands`, `mos2.EIG`, `mos2.XV`, `mos2.STRUCT_OUT`, etc.), junto con los
scripts y figuras de post-procesado (`plot_siesta_bands.py`, `scan_energy.py`,
`*.png`). Se han eliminado los ficheros de reinicio y de estado interno propios
de SIESTA (matrices `.DM`/`.HSX`, ficheros `.ion`/`.psml`, `INPUT_TMP.*`, logs
de `fdf`, etc.), ya que no aportan información reproducible por sí mismos y su
peso conjunto era considerable.

## `AJUSTE/`

Contiene el ajuste de los doce parámetros del modelo tight-binding a las bandas
de referencia de DFT, descrito en la sección 4.2:

- `sk_model.py`: define la geometría del MoS₂ monocapa y construye H(k) según el
  formalismo de Slater-Koster (apéndices A y B de la memoria).
- `fit.py`: plantea el ajuste como un problema de mínimos cuadrados no lineales
  (`scipy.optimize.least_squares`, método TRF) con varios arranques aleatorios,
  y guarda el mejor resultado en `resultados/fit_result.npz`.
- `plot_gap.py`: genera la comparación entre las bandas ajustadas y las de DFT
  (Figura 5.2).
- `datos/kpath.npz`: puntos k y bandas de referencia de DFT usados como entrada.
- `resultados/`: parámetros ajustados (`fitted_parameters.txt`, `fit_result.npz`)
  y la figura de comparación de gap (`gap_comparison.png`).

## `VQD/`

Contiene la resolución del Hamiltoniano tight-binding mediante VQD (capítulo 3
y sección 4.3), dividida en los dos regímenes descritos en la sección 5.3.

### `VQD/SimulacionExacta/`

VQD evaluado mediante multiplicación de matrices (`statevector`), sin ruido de
muestreo:

- `ansatz_particula.py`: construye el *ansatz* `ExcitationPreserving` que
  preserva el número de partículas (sección 4.3).
- `generar_hamiltonianos.py`: traduce H(k) a `SparsePauliOp` mediante la
  transformación de Jordan-Wigner, para los puntos k solicitados.
- `calcular_bandas_vqd.py`: ejecuta VQD con optimizador L-BFGS-B y *warm-start*
  entre puntos k consecutivos, con guardado incremental por si el trabajo se
  interrumpe.
- `representar_bandas.py`, `tiempo_por_punto.py`: generación de las Figuras 5.3
  y 5.6.
- `datos/`: `fit_result.npz` y `kpath.npz`, entrada de `generar_hamiltonianos.py`.
- `resultados/`: Hamiltonianos precalculados (`hamiltonianos_40.pkl`,
  `hamiltonianos_precalculados.pkl`) y figuras correspondientes.

### `VQD/Muestreo/`

VQD evaluado mediante muestreo estadístico (simulador Qulacs, optimizador SPSA),
correspondiente a los tres experimentos descritos en la sección 4.3: la cadena
de *warm-start*, los puntos especiales alrededor de Γ, K y M, y el barrido del
número de *shots*:

- `calcular_bandas_vqd_ruido.py`: script principal, parametrizado por número de
  *shots*, presupuesto de evaluaciones, optimizador y fichero `.pkl` de entrada.
- `warmstart_40.pkl`, `alrededor_Gamma.pkl`, `alrededor_K.pkl`, `alrededor_M.pkl`:
  conjuntos de Hamiltonianos ya precalculados (en formato `SparsePauliOp`) para
  cada uno de los experimentos.
- `datos/`: `fit_result.npz` y `kpath.npz`, comunes al resto de scripts del
  proyecto.
- `slurm/`: un script de envío por experimento (`enviar_ruido_1_warmstart.sh` a
  `enviar_ruido_9_shots16384_K.sh`), con la configuración exacta (shots,
  presupuesto, β) empleada en cada uno.
- `resultados/`: salida en formato JSON de cada experimento, con la energía VQD,
  la energía exacta de referencia, el error, el solapamiento máximo y el vector
  de parámetros θ obtenidos en cada (punto, nivel).
- `representar_bandas_shots.py`, `representar_error_vs_shots.py`: generación de
  las Figuras 5.4 y 5.5.

## Reproducción

Los tres bloques se ejecutan en orden: `AJUSTE/` genera `fit_result.npz`, que es
consumido por `VQD/SimulacionExacta/generar_hamiltonianos.py` para construir los
Hamiltonianos en segunda cuantización, que a su vez son la entrada de
`calcular_bandas_vqd.py` (simulación exacta) o `calcular_bandas_vqd_ruido.py`
(con muestreo).

El entorno de cómputo cuántico requiere `qiskit`, `qiskit-algorithms`, `qulacs`,
además de `numpy` y `scipy`. Los experimentos con muestreo se lanzaron en un
clúster HPC mediante SLURM, como reflejan los scripts de `VQD/Muestreo/slurm/`.

## Referencia

Montiel Estévez, I. (2026). *Simulación cuántica de la estructura electrónica
del MoS₂ mediante VQD*. Trabajo de Fin de Máster, Universidade da Coruña.
