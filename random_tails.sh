#!/bin/bash

declare base

function random_tails {
  local prob

  # start of tail binning
  local tstart="1.5"

  local mode=$1
  local steps=$2

  seq 1 1 99 | while read arg; do
    prob=$(echo "scale=2; ${arg} / 100" | bc)
    ./random --normal $mode $steps "0${prob}" | \
      ./histogram --normal 1 0.1 $tstart > "${base}_${arg}.tail-${tstart}.env"
  done
}

if [[ -n "${1}" ]] && [[ -n "${2}" ]] && [[ -n "${3}" ]]; then
  base="${1}"

  random_tails $2 $3 2> "${base}.status"
fi
