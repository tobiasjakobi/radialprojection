#!/bin/bash

declare base

function cyclotomic_random_envelopes {
  local prob

  local mode=$1
  local steps=$2

  seq 1 1 99 | while read arg; do
    prob=$(echo "scale=2; ${arg} / 100" | bc)
    ./cyclotomic_random --normal $mode $steps "0${prob}" | \
      ./histogram --normal > "${base}_${arg}.env"
  done
}

if [[ -n "${1}" ]] && [[ -n "${2}" ]] && [[ -n "${3}" ]]; then
  base="${1}"

  cyclotomic_random_envelopes $2 $3 2> "${base}.status"
fi
