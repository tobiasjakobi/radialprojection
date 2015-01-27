#!/bin/bash

function randomhist {
  local workdir="${HOME}/PhD_Work/radialprojection/"
  local retval=0
  local tempfile

  local mode=1
  local steps=1000
  local prob="0.5"

  [[ -z "${1}" ]] || mode=$1
  [[ -z "${2}" ]] || steps=$2
  [[ -z "${3}" ]] || prob="${3}"

  if [ ! -d "${workdir}" ]; then
    echo "error: work directory not found"
    return 2
  fi

  pushd "${workdir}" > /dev/null

  if [[ -f ./random ]] && [[ -f ./histogram ]]; then
    tempfile=$(mktemp --tmpdir=/dev/shm)
    if [ $? -ne 0 ]; then
      echo "error: temp file could not be created"
      retval=4
    else
      ./random --single $mode $steps "${prob}" > "${tempfile}"
      
      if [ $? -ne 0 ]; then
        echo "error: call to 'random' failed"
        retval=5
      fi
    fi
  else
    echo "error: required binaries not found"
    retval=3
  fi

  popd > /dev/null

  echo "${tempfile}"

  return $retval
}

retstring=$(randomhist "$@" 2> /dev/null)

if [ $? -eq 0 ]; then
  echo "\"${retstring}\""
else
  [[ "${DEBUG}" -eq 1 ]] && echo "${retstring}"
fi
