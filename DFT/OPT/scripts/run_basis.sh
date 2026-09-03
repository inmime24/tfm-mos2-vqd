#!/bin/bash

export LANG=C
export LC_ALL=C

INPUT="mos2.fdf"
SIESTA="siesta"

OUTFILE="rcut_scan.dat"

echo "# rcut Energy(eV)" > $OUTFILE

for RCUT in $(seq 2.0 0.05 6.8)
do

    DIR=$(printf "rcut_%3.1f" $RCUT)

    echo "======================================"
    echo "Probando rcut = $RCUT"
    echo "======================================"

    mkdir -p $DIR

    cp *.psml $DIR/

    sed "s/RCUT_MO/${RCUT}/g" $INPUT > $DIR/input.fdf

    cd $DIR

    $SIESTA < input.fdf > siesta.out

    ENERGY=$(grep "siesta:         Total =" siesta.out | tail -1 | awk '{print $4}')

    if [ -z "$ENERGY" ]; then
        ENERGY="ERROR"
    fi

    echo "$RCUT   $ENERGY" >> ../$OUTFILE

    cd ..

done

echo ""
echo "Scan terminado."
