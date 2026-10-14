#!/bin/bash
set -eo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
source /scratchfs2/cms/liugeliang/miniconda3/bin/activate unbin
cd "$(dirname "${BASH_SOURCE[0]}")"
python -u produce_signal_reduced.py --era 2018 --cat resolved --input-dir /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/SignalWithDeepJet_20261012 --output-dir /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_signal_20261012
