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
echo "tag,mesh,nBody,nLegs,cells,attempts,coverage_pct,body_layers,body_thk_a,body_thk_b,legs_layers,legs_thk_a,legs_thk_b,stack_mm,last_over_cell,checkMesh,seconds,verdict" > "$out"

for cfg in "${configs[@]}"; do
    nb=${cfg%%:*}; nl=${cfg##*:}; tag="${mesh}_b${nb}_l${nl}"; dir=$src/sweep/$tag
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
    # per-patch rows of the layer table (either snappy format): name faces layers colA colB
    row() { grep -A8 -E 'layers +(overall|avg) +thickness' log.snappyHexMesh 2>/dev/null | awk -v p="$1" '$1==p && NF>=5 {print $3","$4","$5; f=1; exit} END{if(!f) print ",,"}'; }
    b=$(row body); l=$(row legs)
    info=$(grep -E '^  body  n=' sweep.out | head -1)
    stack=$(sed -nE 's/.*total=([0-9.e+-]+).*/\1/p' <<< "$info" | awk '{printf "%.3f", $1*1000}')
    ratio=$(sed -nE 's/.*= ([0-9.]+) x surface cell.*/\1/p' <<< "$info")
    chk=$(grep -q 'Mesh OK' log.checkMesh 2>/dev/null && echo OK || echo FAIL)
    # verdict: coverage >= 90 %, body gets >= 90 % of the requested layers, legs >= 70 %, mesh OK
    bl=${b%%,*}; ll=${l%%,*}
    v=$(awk -v c="${cov:-0}" -v bl="${bl:-0}" -v nb="$nb" -v ll="${ll:-0}" -v nl="$nl" -v k="$chk" \
        'BEGIN{print (c>=90 && bl>=0.9*nb && ll>=0.7*nl && k=="OK") ? "PASS" : "FAIL"}')
    echo "$tag,$mesh,$nb,$nl,${cells:-},$att,${cov:-},$b,$l,$stack,$ratio,$chk,$sec,$v" >> "$out"
    rm -rf processor* 0 constant/polyMesh    # keep only logs (small)
    cd "$src" || exit 1
done

echo; echo "layer_sweep.csv:"
column -s, -t "$out" 2>/dev/null || cat "$out"
echo; echo "Pick the largest PASS row; set nLayersBody / nLayersLegs in case.inputs. Logs: sweep/<tag>/"
