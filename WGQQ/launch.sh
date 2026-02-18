#!/usr/bin/env bash

cd /grid_mnt/data__data.polcms/cms/song/CMSSW_14_1_0_pre4/src/EFT2Obs
source env.sh

python3 WGQQ/launch.py -p 200.0 --step run