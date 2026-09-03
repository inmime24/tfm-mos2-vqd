#!/bin/bash

export LANG=C
export LC_ALL=C

INPUT="mos2.fdf"
SIESTA="siesta"

OUTFILE="kgrid_scan.dat"

echo "# Kgrid Energy(eV)" > $OUTFILE

# Barrido isotrópico de Monkhorst-Pack
for K in $(seq 2 2 16)
do

    DIR=$(printf "kgrid_%02dx%02dx%02d" $K $K 1)

    echo "======================================"
    echo "Probando k-grid = ${K}x${K}x1"
    echo "======================================"

    mkdir -p $DIR

    cp *.psml $DIR/

    # Sustituye el bloque KGRID_BLOCK en el .fdf
    sed "s/KGRID_VALUE/${K}/g" $INPUT > $DIR/input.fdf

    cd $DIR

    $SIESTA < input.fdf > siesta.out

    ENERGY=$(grep "siesta:         Total =" siesta.out | tail -1 | awk '{print $4}')

    if [ -z "$ENERGY" ]; then
        ENERGY="ERROR"
    fi

    echo "${K}x${K}x1   $ENERGY" >> ../$OUTFILE

    cd ..

done

echo ""
echo "Scan de k-grid terminado."
