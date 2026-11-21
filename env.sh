#!/usr/bin/env bash
eftobs_env_cwd=${PWD}
eftobs_env_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
# EFT2Obs is checked out directly under the CMSSW src directory.
cd "${eftobs_env_dir}/.." || return 1
eval `scramv1 runtime -sh`

cd "${eftobs_env_dir}" || return 1
export EFTOBS_LOCAL_LHAPDF=1
export LHAPDF_CONFIG_PATH="${PWD}/lhapdf/bin/lhapdf-config"
export PYTHONPATH="${PWD}/$(echo lhapdf/lib64/python*/site-packages):${PYTHONPATH}"
export LD_LIBRARY_PATH="${PWD}/lhapdf/lib:${LD_LIBRARY_PATH}"
export RIVET_ANALYSIS_PATH=${PWD}/RivetPlugins
export MG_DIR="MG5_aMC_v2_9_16"
export MG_TARBALL="MG5_aMC_v2.9.16.tar.gz"
export RIVET_VERSION="4.1.4"
export DEBUG_SCRIPTS=0

if [ -f "local/rivetenv.sh" ]; then
  source local/rivetenv.sh
fi

if [ "$DEBUG_SCRIPTS" -eq "1" ]; then
	set -x
fi

if [[ ! -z "$PYTHIA8DATA" ]]; then
        export PYTHIA8DATA=""
fi
cd "${eftobs_env_cwd}" || return 1
unset eftobs_env_cwd eftobs_env_dir
#[[ ":$PYTHONPATH:" != *":$PWD/${MG_DIR}:"* ]] && PYTHONPATH="$PWD/${MG_DIR}:${PYTHONPATH}"
