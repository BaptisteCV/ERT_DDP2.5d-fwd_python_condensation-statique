#!/usr/bin/env bash
set -euo pipefail

# Script de benchmark pour testEPM_EigenCS
# Génère des fichiers tmp_Simu.txt à partir du template et lance le binaire

TEMPLATE="../../input/testEPM_EigenCS/Simu.txt"
TMP="tmp_Simu.txt"
EXE="./ElectricalResistivity"
OUT="bench_eigencs_times.csv"
LOGDIR="bench_logs"
mkdir -p "$LOGDIR"

sizes=(50 100 150 200 250 300 350 400)
reps=3

echo "size,rep,real,user,sys,exit_code" > "$OUT"

for s in "${sizes[@]}"; do
  for ((r=1; r<=reps; r++)); do
    # reconstruire tmp_Simu.txt : garder 2 premières lignes, remplacer lignes 3-4 par la taille
    head -n2 "$TEMPLATE" > "$TMP"
    echo "$s" >> "$TMP"
    echo "$s" >> "$TMP"
    tail -n +5 "$TEMPLATE" >> "$TMP" || true

    TIMEFILE="time_${s}_${r}.txt"
    LOGFILE="$LOGDIR/run_${s}_${r}.log"

    pushd ../.. >/dev/null
    /usr/bin/time -f "%e %U %S" -o "electrical_resistivity/Release/$TIMEFILE" "electrical_resistivity/Release/ElectricalResistivity" "electrical_resistivity/Release/$TMP" >"electrical_resistivity/Release/$LOGFILE" 2>&1 || true
    exit_code=$?
    popd >/dev/null
    read real user sys < "electrical_resistivity/Release/$TIMEFILE" || { real=""; user=""; sys=""; }
    echo "${s},${r},${real},${user},${sys},${exit_code}" >> "$OUT"
  done
done

echo "Bench terminé — résultats dans $OUT, logs dans $LOGDIR"
