#!/bin/bash
# run_study.sh - full S -> M -> L study in one command (meshing + solve + convergence + checks)
#   ./run_study.sh [cores]            MESHES="S M" ./run_study.sh 6      (subset)
# Per mesh the key files are archived in results/<mesh>/ ; summary.csv gets one row per run.
cd "$(dirname "$0")" || exit 1
cores=$1; mkdir -p results
for m in ${MESHES:-S M L}; do
    echo "=================== mesh $m  ($(date +%H:%M))"
    ./Allclean > /dev/null 2>&1                      # summary.csv is kept on purpose
    ./Allrun -m "$m" ${cores:+-c $cores} > "results/run_$m.log" 2>&1 || echo "  mesh $m: Allrun returned an error (see results/run_$m.log)"
    d=results/$m; mkdir -p "$d"
    cp -f convergence.png convergence_report.txt log.checkMesh log.snappyHexMesh log.layerSetup \
          system/snappyHexMeshDict case.inputs constant/caseSettings "$d"/ 2>/dev/null
    ./check_workflow.sh > "$d/check.txt" 2>&1; tail -1 "$d/check.txt"
    tail -n 12 convergence_report.txt 2>/dev/null
done
python3 mesh_study.py
