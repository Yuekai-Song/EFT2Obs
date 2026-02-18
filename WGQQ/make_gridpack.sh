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
PROCESS=WGQQ-${name}
GRID_DIR=WGQQ/${name}

if [[ -d "${MG_DIR}/${PROCESS}" ]]; then
    echo "Directory ${MG_DIR}/${PROCESS} already exists. Remove it!"
    rm -r ${MG_DIR}/${PROCESS}
fi
if [[ -d "cards/${PROCESS}" ]]; then
    echo "Directory cards/${PROCESS} already exists. Remove it!"
    rm -r cards/${PROCESS}
fi
cp -r cards/WGQQ cards/${PROCESS}

pushd cards/${PROCESS}
sed -i "119s/^ 200\.0  = pta / ${pta_min}  = pta /" run_card.dat
sed -i "121s/^ -1\.0  = ptamax / ${pta_max}  = ptamax /" run_card.dat
sed -i "8s/^output WGQQ/output ${PROCESS}/" proc_card.dat
popd

./scripts/setup_process.sh $PROCESS

python3 scripts/make_config.py -p ${PROCESS} -o cards/${PROCESS}/config.json \
 --pars SMEFT:2 SMEFTcpv:2  --def-val 0.01 --def-sm 0.0 --def-gen 0.0

python3 scripts/make_param_card.py -p ${PROCESS} -c cards/${PROCESS}/config.json -o cards/${PROCESS}/param_card.dat

python3 scripts/make_reweight_card.py cards/${PROCESS}/config.json cards/${PROCESS}/reweight_card.dat --prepend "change helicity false"

if [[ -d "${GRID_DIR}" ]]; then
    echo "Directory of gridpack ${GRID_DIR} already exists. Remove it!"
    rm -r ${GRID_DIR}
fi

./scripts/make_gridpack.sh ${PROCESS} 0 0 ${GRID_DIR}
rm -r ${MG_DIR}/${PROCESS}
rm -r cards/${PROCESS}