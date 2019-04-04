#!/bin/sh

echo "configuring setup for wls manchester"

source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup
setup root v6_12_04e -q e15:prof

export PATH=${ROOTSYS}/bin:${PWD}/tools:${PATH}
export LD_LIBRARY_PATH=${GEANT4_FQ_DIR}/lib64:${ROOTSYS}/lib:${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}

