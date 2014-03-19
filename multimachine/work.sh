#!/bin/bash

pushd $(dirname $0) > /dev/null

export PATH=$(pwd):$PATH

source "mmachine_scripts"
mmachine_slave_main $*

popd > /dev/null

