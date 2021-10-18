#!/bin/bash -l

module load daint-gpu
module load cray-python
module load cudatoolkit
module load Boost

NVCC_PATH=$(which nvcc)
CUDA_PATH=$(echo $NVCC_PATH | sed -e "s/\/bin\/nvcc//g")
export CUDA_HOME=$CUDA_PATH
export export LD_LIBRARY_PATH=$CUDA_PATH/lib64:$LD_LIBRARY_PATH

export GT_CACHE_ROOT=/scratch/snx3000tds/subbiali/burgers/python/gt4py/gt_cache
export GT_CACHE_DIR_NAME=.gt_cache
