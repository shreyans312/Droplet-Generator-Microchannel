#!/usr/bin/env bash
cd /home/shreyans/Desktop/microfluidic/Droplet-Generator-Microchannel || exit 1
source /opt/openfoam13/etc/bashrc
set -eo pipefail
for name in coarse_uncapped coarse_capped fine_capped; do
    cfd_case="studies/roshan2025-cfd/static-drop/$name"
    if ! test -d "$cfd_case/processor0"; then
        blockMesh -case "$cfd_case" > "$cfd_case/log.blockMesh" 2>&1
        checkMesh -case "$cfd_case" > "$cfd_case/log.checkMesh" 2>&1
        setFields -case "$cfd_case" > "$cfd_case/log.setFields" 2>&1
        decomposePar -case "$cfd_case" > "$cfd_case/log.decomposePar" 2>&1
    fi
    mpirun --oversubscribe -np 8 foamRun -case "$cfd_case" -parallel -solver incompressibleVoF > "$cfd_case/log.foamRun" 2>&1
    reconstructPar -case "$cfd_case" -noFunctionObjects > "$cfd_case/log.reconstruct" 2>&1
done
/usr/bin/python3 studies/roshan2025-cfd/static_drop.py analyse
