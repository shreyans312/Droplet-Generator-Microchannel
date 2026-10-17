#!/usr/bin/env bash
cd /home/shreyans/Desktop/microfluidic/Droplet-Generator-Microchannel || exit 1
source /opt/openfoam13/etc/bashrc
set -eo pipefail
export LD_LIBRARY_PATH="$PWD/studies/roshan2025-cfd/lib:$LD_LIBRARY_PATH"
exec /usr/bin/python3 studies/roshan2025-cfd/automatic_followup.py
