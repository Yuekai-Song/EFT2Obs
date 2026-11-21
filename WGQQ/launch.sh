#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source env.sh

python3 WGQQ/launch.py -p 200.0,250.0,300.0,400.0,500.0,600.0,700.0,800.0,1000.0 --step run
