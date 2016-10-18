#!/bin/bash

export PATH=$PACKAGES_DIR/mpich-3.2/bin:$PATH
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export C_INCLUDE_PATH=$PACKAGES_DIR/mpich-3.2/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$PACKAGES_DIR/mpich-3.2/include:$CPLUS_INCLUDE_PATH
export FPATH=$PACKAGES_DIR/mpich-3.2/include:$FPATH
export MANPATH=$PACKAGES_DIR/mpich-3.2/share/man:$MANPATH
export LD_LIBRARY_PATH=$PACKAGES_DIR/mpich-3.2/lib:$PACKAGES_DIR/mpich-3.2/lib/openmpi:$LD_LIBRARY_PATH