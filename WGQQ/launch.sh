#!/usr/bin/env bash

cd /grid_mnt/data__data.polcms/cms/song/CMSSW_16_0_0/src/EFT2Obs
source env.sh

python3 WGQQ/launch.py -p 200.0,400.0,600.0,800.0,1000.0 --step run