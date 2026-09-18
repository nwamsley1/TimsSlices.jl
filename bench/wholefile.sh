#!/bin/bash
# Whole-file conversion with wall time and peak RSS. Usage: bench/wholefile.sh <run.d> <out_dir> [options]
set -u
cd "$(dirname "$0")/.."
/usr/bin/time -l julia -t ${JULIA_THREADS:-12} --project=. bin/tims_convert.jl "$@" 2>&1 | grep -E "^source|^frames|^CPU|^wrote|real|maximum resident|/ .* frames"
