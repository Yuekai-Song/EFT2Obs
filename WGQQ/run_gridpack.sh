#!/usr/bin/env bash

set -e

cd "$(dirname "${BASH_SOURCE[0]}")/.."
source env.sh
WGQQ_OUTPUT_DIR=$(python3 WGQQ/storage.py)
pta_min=$1
pta_max=$2

if [[ "$pta_max" == "-1.0" ]]; then
    name="ptg${pta_min%.*}"
else
    name="ptg${pta_min%.*}to${pta_max%.*}"
fi
PROCESS=${name}
OUTPUT_DIR=${WGQQ_OUTPUT_DIR}/${PROCESS}
TEST=${3-0}
RIVET="WGQQ"

if [[ $TEST == "1" ]]; then
  python3 scripts/run_gridpack.py --gridpack "${OUTPUT_DIR}/gridpack.tar.gz" -s 101 -e 500 \
  -p WGQQ \
  -o "${OUTPUT_DIR}"
else
  rm -f "${OUTPUT_DIR}"/*.yoda
  mkdir -p jobs/WGQQ_${PROCESS}
  python3 scripts/launch_jobs.py --gridpack "${OUTPUT_DIR}/gridpack.tar.gz" -j 100 -s 1 -e 20000 \
  -p ${RIVET} -o "${OUTPUT_DIR}" --job-mode condor --dir jobs/WGQQ_${PROCESS} --task-name WGQQ-${PROCESS} \
  --sub-opts '+MaxRuntime = 28800\ngetenv = False\nuse_x509userproxy = True\nshould_transfer_files = NO'
fi
