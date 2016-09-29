#!/bin/bash

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90

mkdir -p repos
pushd repos

git clone --depth 100 -b maint https://bitbucket.org/petsc/petsc.git petsc
pushd petsc
./configure \
    PETSC_DIR=${REPOS_DIR}/petsc \
    --prefix=${DEPS_DIR}/petsc \
    --with-c2html=0 \
    --with-x=0 \
    --with-clanguage=C \
    --with-mpicompilers=1 \
    --with-debugging=1 \
    --with-shared-libraries=1 \
    --with-64-bit-points=1 \
    --with-large-file-io=1 \
    --download-f2cblaslapack=1
make PETSC_DIR=${REPOS_DIR}/petsc all
make PETSC_DIR=${REPOS_DIR}/petsc install
make PETSC_DIR=${DEPS_DIR}/petsc PETSC_ARCH="" test
popd

git clone --depth 100 -b v1.0.0 https://github.com/libMesh/libmesh.git libmesh
pushd libmesh
METHODS="dbg" ./configure --prefix=${DEPS_DIR}/libmesh
make -j${MAKE_THREADS}
make install
popd

# $HOME
popd
