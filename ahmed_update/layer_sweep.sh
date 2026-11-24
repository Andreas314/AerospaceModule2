#!/bin/bash
# layer_sweep.sh - how many prism layers can we really add?  (meshing only, no solve)
#   ./layer_sweep.sh [-m S|M|L] [-c cores] [body:legs ...]
#   ./layer_sweep.sh -m M -c 6 10:2 15:2 20:2 22:4 22:8
# Each config runs ./Allrun.pre in its own scratch copy (sweep/<tag>/), so the real
# case is never touched. Results -> layer_sweep.csv + a table at the end.
cd "$(dirname "$0")" || exit 1
mesh=M; cores=
while [ "$#" -gt 0 ]; do case "$1" in
    -m) mesh=$2; shift 2 ;;  -c) cores=$2; shift 2 ;;  -h|--help) sed -n '2,7p' "$0"; exit 0 ;;  *) break ;;
esac; done
configs=("$@"); [ "${#configs[@]}" -eq 0 ] && configs=(10:2 15:2 20:2 22:2)
src=$PWD; out=$src/layer_sweep.csv
[ -f "$out" ] || echo "tag,mesh,nBody,nLegs,cells,attempts,coverage_pct,body_layers,body_thk_m,body_thk_pct,legs_layers,legs_thk_m,legs_thk_pct,stack_mm,last_over_cell,checkMesh,hard,soft,volratio_faces,nonorth_max,seconds,verdict" > "$out"

for cfg in "${configs[@]}"; do
    nb=${cfg%%:*}; nl=${cfg##*:}; tag="${mesh}_b${nb}_l${nl}${LAYER_SET:+_set$(echo $LAYER_SET | tr " " "-")}"; dir=$src/sweep/$tag
    rm -rf "$dir"; mkdir -p "$dir"
    tar -C "$src" --exclude=./sweep --exclude='./processor*' --exclude=./postProcessing \
        --exclude=./results --exclude='./log.*' --exclude=./0 -cf - . | tar -C "$dir" -xf -
    sed -i -E "s/^nLayersBody=[0-9]+/nLayersBody=$nb/; s/^nLayersLegs=[0-9]+/nLayersLegs=$nl/" "$dir/case.inputs"
    echo ">>> $tag  (body $nb, legs $nl layers, mesh $mesh)"
    t0=$(date +%s)
    ( cd "$dir" && ./Allrun.pre -m "$mesh" ${cores:+-c $cores} > sweep.out 2>&1 )
    sec=$(( $(date +%s) - t0 ))
    cd "$dir" || continue
    cells=$(awk '/^ *cells:/ {print $2; exit}' log.checkMesh 2>/dev/null)
    att=$(ls log.snappyHexMesh.try* 2>/dev/null | wc -l)
    cov=$(grep -E 'Added [0-9]+ out of [0-9]+ cells' log.snappyHexMesh 2>/dev/null | tail -1 | sed -E 's/.*\(([0-9.eE+-]+)%\).*/\1/')
    # per-patch row of the ACHIEVED layer table ("layers overall thickness"): avg layers, thickness [m], [%]
    # (snappy also prints a table of the REQUESTED layers earlier in the log - do not use that one)
    row() { awk -v p="$1" '/layers +overall +thickness/ {on=NR; next} on && NR<=on+8 && $1==p && NF>=5 {r=$3","$4","$5} END{print (r=="" ? ",," : r)}' log.snappyHexMesh 2>/dev/null; }
    b=$(row body); l=$(row legs)
    info=$(grep -E '^  body  n=' sweep.out | head -1)
    stack=$(sed -nE 's/.*total=([0-9.e+-]+).*/\1/p' <<< "$info" | awk '{printf "%.3f", $1*1000}')
    ratio=$(sed -nE 's/.*= ([0-9.]+) x surface cell.*/\1/p' <<< "$info")
    # checkMesh: "soft" = stiffness warnings that every y+~1 prism-layer mesh triggers
    # (small determinant, face tets, concave, interpolation weight, volume ratio);
    # "hard" = any other *** line (skewness, openness, negative volume, non-orthogonality...)
    soft_re='face tets|determinant|concave|interpolation weight|volume ratio|flatness|warpage'
    hard=$(grep -E '^ *\*\*\*' log.checkMesh 2>/dev/null | grep -Eiv "$soft_re" | wc -l)
    soft=$(grep -E '^ *\*\*\*' log.checkMesh 2>/dev/null | grep -Eic "$soft_re")
    vr=$(grep 'small volume ratio' log.checkMesh 2>/dev/null | sed -E 's/.*faces: *([0-9]+).*/\1/'); vr=${vr:-0}
    no=$(sed -nE 's/.*non-orthogonality Max: *([0-9.]+).*/\1/p' log.checkMesh 2>/dev/null | head -1)
    chk=$([ -s log.checkMesh ] && [ "$hard" -eq 0 ] && echo OK || echo FAIL)
    # verdict: coverage >= 90 %, body gets >= 90 % of the requested layers, legs >= 70 %, no HARD checkMesh failure
    bl=${b%%,*}; ll=${l%%,*}
    v=$(awk -v c="${cov:-0}" -v bl="${bl:-0}" -v nb="$nb" -v ll="${ll:-0}" -v nl="$nl" -v k="$chk" \
        'BEGIN{print (c>=90 && bl>=0.9*nb && ll>=0.7*nl && k=="OK") ? "PASS" : "FAIL"}')
    echo "$tag,$mesh,$nb,$nl,${cells:-},$att,${cov:-},$b,$l,$stack,$ratio,$chk,$hard,$soft,$vr,${no:-},$sec,$v" >> "$out"
    rm -rf processor* 0 constant/polyMesh    # keep only logs (small)
    cd "$src" || exit 1
done

echo; echo "layer_sweep.csv:"
column -s, -t "$out" 2>/dev/null || cat "$out"
echo; echo "Pick the largest PASS row; set nLayersBody / nLayersLegs in case.inputs. Logs: sweep/<tag>/"
