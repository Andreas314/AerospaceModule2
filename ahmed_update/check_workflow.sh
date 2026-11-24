#!/bin/bash
# check_workflow.sh - pass/fail audit of the Ahmed-body workflow.
# Run it (1) before anything, (2) after ./Allrun.pre, (3) after ./Allrun. It checks
# what exists at that moment and prints  OK / FAIL / WARN  with a hint for every FAIL.
cd "$(dirname "$0")" || exit 1
fails=0; warns=0
ok()   { printf '  [ OK ] %s\n' "$1"; }
fail() { printf '  [FAIL] %s\n         -> %s\n' "$1" "$2"; fails=$((fails+1)); }
warn() { printf '  [WARN] %s\n         -> %s\n' "$1" "$2"; warns=$((warns+1)); }
between() { awk -v x="$1" -v a="$2" -v b="$3" 'BEGIN{exit !(x+0>=a && x+0<=b)}'; }

echo "== 1. environment"
[ "${WM_PROJECT_VERSION}" = "14" ] && ok "OpenFOAM 14 loaded" || warn "WM_PROJECT_VERSION='${WM_PROJECT_VERSION}'" "source ~/OpenFOAM/OpenFOAM-14/etc/bashrc"
for t in blockMesh surfaceFeatures decomposePar snappyHexMesh checkMesh potentialFoam foamRun foamDictionary reconstructPar mpirun; do
    command -v $t >/dev/null && ok "tool $t" || fail "tool $t not found" "OpenFOAM environment not loaded?"
done
python3 -c "import numpy, matplotlib" 2>/dev/null && ok "python3 numpy + matplotlib" || fail "numpy/matplotlib missing" "pip install numpy matplotlib"

echo "== 2. inputs and geometry"
if source ./case.inputs 2>/dev/null; then ok "case.inputs readable (U=$U, y+=$yplus, layers body/legs=$nLayersBody/$nLayersLegs, mesh=$meshSize)"; else fail "case.inputs" "missing or syntax error"; fi
read -r L A <<< "$(python3 stl_info.py constant/geometry/ahmedBody.stl body 2>/dev/null)"
between "$L" 1.03 1.06 && ok "body length L=$L m (Ahmed: 1.044)" || fail "body length '$L'" "STL missing / wrong units"
between "$A" 0.054 0.058 && ok "half-model frontal area Aref=$A m2 (Ahmed full body 0.112 / 2)" || warn "Aref='$A'" "expected ~0.056 for the half model"

if [ -f log.snappyHexMesh ]; then
echo "== 3. mesh (after Allrun.pre)"
    grep -q 'Finished meshing' log.snappyHexMesh && ok "snappyHexMesh finished" || fail "snappyHexMesh did not finish" "see log.snappyHexMesh*"
    tgt=$(grep -m1 '^coverageTarget=' Allrun.pre | sed -E 's/coverageTarget=([0-9.]+).*/\1/')
    cov=$(grep -E 'Added [0-9]+ out of [0-9]+ cells' log.snappyHexMesh | tail -1 | sed -E 's/.*\(([0-9.eE+-]+)%\).*/\1/')
    between "${cov:-0}" "${tgt:-90}" 100.1 && ok "layer coverage ${cov}% >= ${tgt}%" || fail "layer coverage ${cov:-?}% < ${tgt}%" "fewer layers, or see layer_sweep.sh"
    ntry=$(ls log.snappyHexMesh.try* 2>/dev/null | wc -l)
    [ "$ntry" -le 1 ] && ok "layers added at the first attempt" || warn "needed $ntry snappy attempts" "later attempts use featureAngle up to 180 / relaxed limits -> check cell quality"
    soft_re='face tets|determinant|concave|interpolation weight|volume ratio|flatness|warpage'
    hardn=$(grep -E '^ *\*\*\*' log.checkMesh 2>/dev/null | grep -Eiv "$soft_re" | wc -l)
    softn=$(grep -E '^ *\*\*\*' log.checkMesh 2>/dev/null | grep -Eic "$soft_re")
    if [ ! -s log.checkMesh ]; then fail "log.checkMesh missing" "checkMesh did not run"
    elif [ "$hardn" -gt 0 ]; then fail "checkMesh: $hardn HARD failure(s) (skewness / openness / negative volume / non-orthogonality)" "grep -n '\*\*\*' log.checkMesh"
    elif [ "$softn" -gt 0 ]; then warn "checkMesh: $softn soft stiffness warning(s) (small determinant / face tets / volume ratio ...)" "typical for y+~1 prism layers; fewer with more layers (smaller size jump)"
    else ok "checkMesh: Mesh OK"; fi
    cells=$(awk '/^ *cells:/ {print $2; exit}' log.checkMesh 2>/dev/null)
    [ "${cells:-0}" -gt 30000 ] && ok "mesh has $cells cells (more than the bare background mesh)" || fail "only ${cells:-0} cells" "the layered mesh is not the one checked: add -overwrite to snappyHexMesh"
    if ls -d processor0/constant/polyMesh >/dev/null 2>&1; then
        solverCells=$(grep -ho 'nCells:[0-9]*' processor*/constant/polyMesh/owner 2>/dev/null | cut -d: -f2 | awk '{s+=$1} END{print s+0}')
        [ "$solverCells" = "${cells}" ] && ok "solver mesh = checked mesh ($solverCells cells)" || fail "solver mesh has $solverCells cells, checkMesh saw ${cells:-?}" "snappyHexMesh wrote 1/ 2/ 3/ time dirs: add -overwrite to 'runp snappyHexMesh' in Allrun.pre"
        for p in inlet outlet bottom top midplane sides body legs; do grep -q "^ *$p\$" processor0/constant/polyMesh/boundary || fail "patch $p missing" "check blockMesh/snappy patch names vs 0.orig"; done
        grep -q 'symmetryPlane' processor0/constant/polyMesh/boundary && ok "all 8 patches present, symmetryPlane found"
    fi
fi

if [ -f log.foamRun ]; then
echo "== 4. solution (after Allrun)"
    grep -qE '^End' log.foamRun && ok "foamRun finished" || fail "foamRun did not reach 'End'" "tail log.foamRun"
    grep -qi 'nan\|floating point exception' log.foamRun && fail "NaN / FPE in log.foamRun" "lower relaxation factors / check mesh quality" || ok "no NaN / FPE"
    cd_=$(awk '$1=="Cd"{print $2; exit}' convergence_report.txt 2>/dev/null)
    between "${cd_:-0}" 0.25 0.33 && ok "Cd = $cd_ (experiment, 25 deg slant: ~0.285)" || fail "Cd = '${cd_:-none}' outside 0.25-0.33" "wrong Aref/CofR/symmetry, or not converged"
    grep -q 'NOT converged' convergence_report.txt 2>/dev/null && warn "force coefficients not converged by the 1% rule" "more iterations (endTime) or accept with std shown"
    yp=$(grep 'y+ on body' convergence_report.txt 2>/dev/null | sed -nE 's/.* avg +([0-9.eE+-]+).*/\1/p' | head -1)
    if [ -n "$yp" ]; then between "$yp" 0.3 2 && ok "mean y+ on body = $yp (target $yplus)" || warn "mean y+ on body = $yp (target $yplus)" "re-mesh with yplus_new = yplus * target/achieved"; else warn "no y+ in report" "yPlus function object / log pattern"; fi
fi

if [ -f summary.csv ]; then
echo "== 5. bookkeeping"
    awk -F, 'NR==1{n=NF} NF!=n{bad=1} END{exit bad}' summary.csv && ok "summary.csv rows all have the same number of columns" || fail "summary.csv misaligned" "use the fixed plot_convergence.py (commit 1)"
fi
echo; echo "RESULT: $fails FAIL, $warns WARN"; [ "$fails" -eq 0 ]
