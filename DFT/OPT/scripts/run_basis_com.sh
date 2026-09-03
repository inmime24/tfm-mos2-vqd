#!/bin/bash

export LANG=C
export LC_ALL=C

INPUT="mos2.fdf"
SIESTA="siesta"

# Valores originales
BASE=(
2.522
6.155
3.648
2.827
6.155
4.308
2.111
3.662
1.928
4.343
2.265
)

# rango de scan
VALUES=$(seq 2.0 0.2 6.0)

# Loop sobre los 11 parámetros
for IDX in $(seq 1 11)
do

    OUTFILE="scan_RCUT${IDX}.dat"

    echo "# rcut energy(eV)" > $OUTFILE

    echo "===================================="
    echo "Escaneando RCUT${IDX}"
    echo "===================================="

    for RCUT in $VALUES
    do

        DIR=$(printf "RCUT%d_%3.1f" $IDX $RCUT)

        mkdir -p $DIR

        cp *.psml $DIR/

        # construir sed dinámico
        CMD="sed "

        for J in $(seq 1 11)
        do

            if [ "$J" -eq "$IDX" ]; then
                VALUE=$RCUT
            else
                VALUE=${BASE[$((J-1))]}
            fi

            CMD="$CMD -e 's/RCUT${J}/${VALUE}/g'"
        done

        CMD="$CMD $INPUT > $DIR/input.fdf"

        eval $CMD

        cd $DIR

        $SIESTA < input.fdf > siesta.out

        ENERGY=$(grep "siesta:         Total =" siesta.out | tail -1 | awk '{print $4}')

        if [ -z "$ENERGY" ]; then
            ENERGY="ERROR"
        fi

        echo "$RCUT   $ENERGY" >> ../$OUTFILE

        cd ..

        echo "RCUT${IDX} = $RCUT   ->   $ENERGY"

    done

done

echo ""
echo "Todos los scans terminados."
