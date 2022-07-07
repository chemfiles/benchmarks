#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

set -exu

conda create -n chemfiles-benchmarks -y
conda activate chemfiles-benchmarks

conda install MDAnalysis openbabel ase cython chemfiles mdtraj parmed rdkit -y
conda update --all -y

pip install --upgrade tabulate

ROOT=`pwd`

mkdir -p pytraj-build && cd pytraj-build

if [[ ! -d cpptraj ]]; then
    git clone https://github.com/Amber-MD/cpptraj
fi
cd cpptraj

if [[ `uname -s` == "Darwin" ]]; then
    yes | bash configure -shared -macAccelerate -noarpack -nosanderlib clang
elif [[ `uname -s` == "Linux" ]]; then
    yes | bash configure -shared -noarpack -nosanderlib -openmp gnu
else
    echo "unknown machine"
    exit 1
fi

make -j8 libcpptraj

export CPPTRAJHOME=`pwd`

cd $ROOT/pytraj-build
if [[ ! -d pytraj ]]; then
    git clone https://github.com/Amber-MD/pytraj
fi
cd pytraj

python setup.py install
