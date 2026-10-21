#!/bin/sh
# Adapted launch recipe from the historical WT/Anc 50-ns PJM scripts.
# Supply site-specific PJM resource/account/module settings outside this file.
# Historical allocation: 4 nodes, 32 MPI processes, 6 OpenMP threads/process.
# This file was syntax-checked only; no job is submitted by this update.
set -eu
: "${GENESIS_BIN:?Set GENESIS_BIN to the approved production executable directory}"
: "${GREST_INPUT:?Set GREST_INPUT to a private reviewed WT or Anc production input}"
test -x "$GENESIS_BIN/spdyn"
test -f "$GREST_INPUT"
export OMP_NUM_THREADS=6
export PLE_MPI_STD_EMPTYFILE=off
exec mpiexec "$GENESIS_BIN/spdyn" "$GREST_INPUT"
