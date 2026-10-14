#!/bin/bash
set -eo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
source /scratchfs2/cms/liugeliang/miniconda3/bin/activate unbin
cd "$(dirname "${BASH_SOURCE[0]}")"
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python -u produce_reduced_stream.py --era 2018 --cat merged --config-dir cards --sf-dir sf --output-dir /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet --step-size 100000
