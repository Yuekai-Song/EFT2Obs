#!/usr/bin/env bash

if [[ $# -lt 2 ]]; then
    echo "Insufficient number of arguments, usage is ./run_NLO.sh events seed"
    exit 1
fi

EVTS=$1
SEED=$2

### SET ENVIRONMENT VARIABLES HERE
RUNLABEL="GridRun"

{
  echo "generate_events -oxR --nb_core=1 -n ${RUNLABEL}"
  echo "shower=ON"
  echo "reweight=ON"
  echo "done"
  echo "set nevents ${EVTS}"
  echo "set iseed ${SEED}"
  echo "done"
} > mgrunscript

./bin/aMCatNLO < mgrunscript
