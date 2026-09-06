#!/bin/bash
#SBATCH -p ilk
#SBATCH --mem-per-cpu=1G
#SBATCH --time=20:00:00
#SBATCH -o ruido_9_shots16384_K_%j.out
#SBATCH -e ruido_9_shots16384_K_%j.err

source /home/usc/cursos/curso181/miniconda3/etc/profile.d/conda.sh
conda activate tfm_env
cd $SLURM_SUBMIT_DIR
pip install qulacs qiskit-algorithms --quiet 2>&1 | tail -5

python3 -u calcular_bandas_vqd_ruido.py \
    --pkl datos/hamiltonianos_precalculados.pkl \
    --n-estados 11 \
    --shots 16384 \
    --metodo SPSA \
    --presupuesto 3000 \
    --beta 10 \
    --salida resultados/ruido_9_shots16384_K.json
