#!/usr/bin/env bash

cd /grid_mnt/data__data.polcms/cms/song/CMSSW_14_1_0_pre4/src/EFT2Obs
source env.sh
pta_min=$1
pta_max=$2

if [[ "$pta_max" == "-1.0" ]]; then
    name="ptg${pta_min%.*}"
else
    name="ptg${pta_min%.*}to${pta_max%.*}"
fi
PROCESS=${name}
TEST=${3-0}
RIVET="WGQQ"

if [[ $TEST == "1" ]]; then
  python scripts/run_gridpack.py --gridpack WGQQ/${PROCESS}/gridpack.tar.gz -s 101 -e 500 \
  -p WGQQ \
  -o WGQQ/${PROCESS}
else
  rm -f WGQQ/${PROCESS}/*.yoda
  python scripts/launch_jobs.py --gridpack WGQQ/${PROCESS}/gridpack.tar.gz -j 100 -s 1 -e 20000 \
  -p ${RIVET} -o WGQQ/${PROCESS} --job-mode condor --dir jobs --task-name WGQQ-${PROCESS} \
  --sub-opts '+JobFlavour = "longlunch"\nT3Queue = long\ngetenv = False\nuse_x509userproxy = True\nWNTag=el9\n+SingularityCmd = ""\ninclude : /opt/exp_soft/cms/t3_tst/t3queue |'
fi