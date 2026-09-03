#!/bin/bash

OUTPUT="scf_results.csv"
INPUT="mos2.fdf"
TMP="tmp.fdf"
OUT="out.out"

echo "mix_type,method,weight,history,n_iter" > $OUTPUT

methods=("Linear" "Pulay" "Broyden")
weights=("0.05" "0.1" "0.2" "0.4" "0.6" "0.9")
history=("3" "5" "8")

mix_types=("Hamiltonian" "Density")

for mix in "${mix_types[@]}"; do
for m in "${methods[@]}"; do
for w in "${weights[@]}"; do
for h in "${history[@]}"; do

echo "Running $mix $m w=$w h=$h"

# copiar base input
cp $INPUT $TMP

# ----------------------------
# modificar mixer
# ----------------------------

if [ "$m" = "Linear" ]; then
    mixer_line="DM.Mixer Method Linear"
elif [ "$m" = "Pulay" ]; then
    mixer_line="DM.Mixer Method Pulay"
else
    mixer_line="DM.Mixer Method Broyden"
fi

# reemplazar o añadir líneas
sed -i "/DM.Mixer/d" $TMP
sed -i "/DM.MixingWeight/d" $TMP
sed -i "/DM.NumberPulay/d" $TMP

echo "$mixer_line" >> $TMP
echo "DM.MixingWeight $w" >> $TMP
echo "DM.NumberPulay $h" >> $TMP

# ----------------------------
# run SIESTA
# ----------------------------
siesta < $TMP > $OUT

# ----------------------------
# extract iterations
# ----------------------------
iter=$(grep -i "iter" $OUT | tail -1 | awk '{print $NF}')

if [ -z "$iter" ]; then
iter="NA"
fi

echo "$mix,$m,$w,$h,$iter" >> $OUTPUT

done
done
done
done

echo "DONE -> $OUTPUT"
