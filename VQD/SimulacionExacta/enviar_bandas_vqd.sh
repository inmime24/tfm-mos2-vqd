#!/bin/bash
#SBATCH -p ilk
#SBATCH --mem-per-cpu=1G
#SBATCH --time=10:00:00
#SBATCH -o bandas_vqd11_%j.out
#SBATCH -e bandas_vqd11_%j.err

source /home/usc/cursos/curso181/miniconda3/etc/profile.d/conda.sh
conda activate tfm_env
cd $SLURM_SUBMIT_DIR

# Si el trabajo se corta por el limite de tiempo, simplemente vuelve a
# lanzar este mismo script con "sbatch" -- el guardado incremental hace
# que continue exactamente donde se quedo, sin repetir nada.
python3 calcular_bandas_vqd.py \
    --pkl datos/hamiltonianos_40.pkl \
    --n-estados 11 \
    --reps 2 \
    --metodo L-BFGS-B \
    --maxiter 5000 \
    --beta 10 \
    --salida resultados/bandas_vqd11.json
