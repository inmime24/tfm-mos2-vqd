#!/bin/bash

# =========================================================
# OPTIMIZACIÓN RCUT PRINCIPALES (r1)
# MoS2 - SIESTA
# =========================================================

SIESTA="siesta"
TEMPLATE="mos2.fdf"

mkdir -p res_basis

RESULTS="res_basis_op.dat"

# tolerancia convergencia (meV)
TOL=5.0

# =========================================================
# VALORES A PROBAR
# =========================================================

Mo4s_vals=(2.3 2.4 2.5 2.6 2.7)

Mo5s_vals=(5.8 5.9 6.0 6.1 6.2)

Mo4p_vals=(2.6 2.7 2.8 2.9 3.0)

Mo5p_vals=(5.8 5.9 6.0 6.1 6.2)

Mo4d_vals=(4.0 4.2 4.4 4.6 4.8)

S3s_vals=(3.2 3.4 3.6 3.8 4.0)

S3p_vals=(4.0 4.2 4.4 4.6 4.8)

# =========================================================
# CABECERA
# =========================================================

echo "# Mo4s Mo5s Mo4p Mo5p Mo4d S3s S3p Energy_eV" \
> $RESULTS

counter=0

# =========================================================
# LOOPS
# =========================================================

for Mo4s in "${Mo4s_vals[@]}"
do
for Mo5s in "${Mo5s_vals[@]}"
do
for Mo4p in "${Mo4p_vals[@]}"
do
for Mo5p in "${Mo5p_vals[@]}"
do
for Mo4d in "${Mo4d_vals[@]}"
do
for S3s in "${S3s_vals[@]}"
do
for S3p in "${S3p_vals[@]}"
do

    runname="run_${counter}"

    rundir="runs/${runname}"

    mkdir -p $rundir

    echo
    echo "================================="
    echo "$runname"
    echo "================================="

    # -----------------------------------------------------
    # GENERAR INPUT
    # -----------------------------------------------------

    sed \
    -e "s/__Mo4s_r1__/${Mo4s}/g" \
    -e "s/__Mo5s_r1__/${Mo5s}/g" \
    -e "s/__Mo4p_r1__/${Mo4p}/g" \
    -e "s/__Mo5p_r1__/${Mo5p}/g" \
    -e "s/__Mo4d_r1__/${Mo4d}/g" \
    -e "s/__S3s_r1__/${S3s}/g" \
    -e "s/__S3p_r1__/${S3p}/g" \
    $TEMPLATE > ${rundir}/mos2.fdf

    # -----------------------------------------------------
    # COPIAR PSEUDOS
    # -----------------------------------------------------

    cp *.psml $rundir/

    # -----------------------------------------------------
    # EJECUTAR
    # -----------------------------------------------------

    cd $rundir

    $SIESTA < mos2.fdf > mos2.out

    cd ../../

    # -----------------------------------------------------
    # EXTRAER ENERGÍA
    # -----------------------------------------------------

    energy=$(grep "siesta: E_KS(eV)" \
    ${rundir}/mos2.out | tail -1 | awk '{print $NF}')

    # fallback
    if [ -z "$energy" ]
    then
        energy=999999
    fi

    # -----------------------------------------------------
    # GUARDAR
    # -----------------------------------------------------

    echo "$Mo4s $Mo5s $Mo4p $Mo5p $Mo4d $S3s $S3p $energy" \
    >> $RESULTS

    counter=$((counter+1))

done
done
done
done
done
done
done

# =========================================================
# ORDENAR
# =========================================================

sort -k8 -n $RESULTS > sorted_results.dat

# =========================================================
# MEJOR ENERGÍA
# =========================================================

best=$(awk 'NR==2 {print $8}' sorted_results.dat)

# =========================================================
# TABLA FINAL
# =========================================================

echo "# Mo4s Mo5s Mo4p Mo5p Mo4d S3s S3p Energy_eV DeltaE_meV Conv" \
> final_results.dat

awk -v best="$best" -v tol="$TOL" '

BEGIN{
OFS=" "
}

NR>1 {

delta=($8-best)*1000.0

if(delta<0)
delta=-delta

conv="NO"

if(delta<tol)
conv="YES"

print $1,$2,$3,$4,$5,$6,$7,$8,delta,conv

}

' sorted_results.dat >> final_results.dat

# =========================================================
# MOSTRAR
# =========================================================

echo
echo "================================="
echo "MEJOR ENERGÍA"
echo "================================="

echo "$best eV"

echo
echo "================================="
echo "CONFIGURACIONES CONVERGIDAS"
echo "================================="

grep YES final_results.dat

echo
echo "================================="
echo "TOP 20"
echo "================================="

head -20 final_results.dat
