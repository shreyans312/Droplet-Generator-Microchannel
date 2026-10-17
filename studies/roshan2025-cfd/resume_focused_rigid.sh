#!/usr/bin/env bash

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
case_dir="studies/roshan2025-cfd/runs-focused/unstretched_robust"
exec > "$case_dir/log.launch.service" 2>&1

source /opt/openfoam13/etc/bashrc
set -eo pipefail
printf 'Starting CFD from configured checkpoint in %s\n' "$case_dir"
mpirun --oversubscribe -np 8 foamRun \
    -case "$case_dir" -parallel -solver incompressibleVoF \
    > "$case_dir/log.foamRun.service" 2>&1

python3 studies/roshan2025-cfd/finalize.py \
    --case "$case_dir" --log log.foamRun.service --parallel \
    --deadline-hours 24 > "$case_dir/log.finalize.service" 2>&1
