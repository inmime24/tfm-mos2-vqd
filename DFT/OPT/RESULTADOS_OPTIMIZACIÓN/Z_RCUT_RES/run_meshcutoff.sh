#!/bin/bash

export LANG=C
export LC_ALL=C

INPUT="mos2.fdf"
SIESTA="siesta"

OUTFILE="meshcutoff_scan.dat"

echo "# MeshCutoff(Ry) Energy(eV)" > $OUTFILE

# Barrido de MeshCutoff desde 50 hasta 500 Ry
for MCUT in $(seq 50 25 500)
do

    DIR=$(printf "meshcutoff_%03d" $MCUT)

    echo "======================================"
    echo "Probando MeshCutoff = ${MCUT} Ry"
    echo "======================================"

    mkdir -p $DIR

    cp *.psml $DIR/

    # Sustituye el marcador MESH_CUTOFF en el .fdf
    sed "s/MESH_CUTOFF/${MCUT} Ry/g" $INPUT > $DIR/input.fdf

    cd $DIR

    $SIESTA < input.fdf > siesta.out

    ENERGY=$(grep "siesta:         Total =" siesta.out | tail -1 | awk '{print $4}')

    if [ -z "$ENERGY" ]; then
        ENERGY="ERROR"
    fi

    echo "$MCUT   $ENERGY" >> ../$OUTFILE

    cd ..

done

echo ""
echo "Scan de MeshCutoff terminado."
