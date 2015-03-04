#!/bin/bash

base="out/decagonal_visrnd"
mode=2
steps=800

function cyclotomic_random_envelopes {
  local prob

  seq 1 1 99 | while read arg; do
    prob=$(echo "scale=2; ${arg} / 100" | bc)
    ./cyclotomic_random --normal $mode $steps "0${prob}" | \
      ./histogram --normal > "${base}_${arg}.env"
  done
}

cyclotomic_random_envelopes 2> "${base}.status"
