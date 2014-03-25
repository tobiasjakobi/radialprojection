#!/bin/bash

pushd $(dirname $0) > /dev/null

export PATH=$(pwd):$PATH

source "mmachine_scripts"

if [ -z $1 ]; then
  echo "status: fetching patch data..."
  mmachine_fetchdata
  if [ $? -ne 0 ]; then
    echo "error: fetching patch data failed, aborting..."
    exit 1
  else
    echo "status: patch data successfully fetched"
  fi
fi

mmachine_slave_main $*
if [ $? -ne 0 ]; then
  echo "error: slave main routine failed"
  exit 1
else
  echo "status: slave main routine successfully called"
fi

popd > /dev/null
