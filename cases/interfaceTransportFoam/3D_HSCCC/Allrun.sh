#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
cp -r 0.org 0
cp -r ../coil_blockMesh/constant/polyMesh constant/.
setFields

decomposePar

runParallel interfaceTransportFoam


#------------------------------------------------------------------------------
