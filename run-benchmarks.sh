#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

set -exu

conda activate chemfiles-benchmarks

python bench.py
