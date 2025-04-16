For `WGto2QG`, the `proc_card.dat` is 
```
import model loop_sm-ckm_no_b_mass

define ww = w+ w-
generate p p > ww a [QCD] @0
add process p p > ww a j [QCD] @1
```
while for `GJ_4Jets`, it's
```
import model sm-ckm_no_b_mass

generate p p > a j @0
add process p p > a j j @1
add process p p > a j j j @2
add process p p > a j j j j @3
```
it looks like there are some overlap between the processes of `p p > ww a` and `p p > a j j`, but actually this is due to the model.

for WGQQ we should use a `proc_card.dat` like: 
```
import model SMEFTsim_topU3l_MwScheme_UFO-massless_cW

define ww = w+ w-

generate p p > ww a NP<=1, ww > j j @0
add process p p > ww a j NP<=1, ww > j j @1

output WGQQ
```

## Then make gridpack:

```shell
./scripts/setup_process.sh WGQQ

python scripts/make_config.py -p WGQQ -o cards/WGQQ/config.json \
 --pars SMEFT:2 SMEFTcpv:2  --def-val 0.01 --def-sm 0.0 --def-gen 0.0

python scripts/make_param_card.py -p WGQQ -c cards/WGQQ/config.json -o cards/WGQQ/param_card.dat

python scripts/make_reweight_card.py cards/WGQQ/config.json cards/WGQQ/reweight_card.dat

./scripts/make_gridpack.sh WGQQ 0 64
```
We have to make some modifications: for `reweight_card.dat`, change the first line to 
```
change rwgt_dir rwgt
change helicity false
```
although we don't think it matters too much.

And for `run_card.dat`, two changes that we have to make are:
- `True = use_syst`: change to `False = use_syst`
- `systematics = systematics_program`: change to `none = systematics_program`

and change some cut criteria:
- `200.0 = pta`
- `6800.0     = ebeam1  ! beam 1 total energy in GeV`
- `6800.0     = ebeam2  ! beam 2 total energy in GeV`






## To run gridpack:
```shell
python scripts/run_gridpack.py --gridpack WGQQ/gridpack.tar.gz -s 1 -e 500 \
  -p WGQQ \
  -o WGQQ
```
or use `condor` to submit jobs
```shell
python scripts/launch_jobs.py --gridpack WGQQ/gridpack.tar.gz -j 100 -s 1 -e 10000 \
-p WGQQ -o WGQQ --job-mode condor --dir jobs --task-name wgqq \
--sub-opts '+JobFlavour = "longlunch"\nT3Queue = long\ngetenv = False\nuse_x509userproxy = True\nWNTag=el9\n+SingularityCmd = ""\ninclude : /opt/exp_soft/cms/t3_tst/t3queue |'
```

## To draw plots:
```shell
yodamerge -o WGQQ/Rivet.yoda WGQQ/Rivet_*

python scripts/get_scaling.py -c cards/WGQQ/config.json \
  -i WGQQ/Rivet.yoda --hist "/WGQQ/p_pt" \
  --save json,txt,tex --translate-tex resources/translate_tex.json --dir WGQQ

python scripts/makePlot.py --hist WGQQ/WGQQ_p_pt.json -c cards/WGQQ/config.json --x-title "p_{T}(\gamma) [GeV]" --title-left "W\gamma #rightarrow qq\gamma" --title-right "SMEFTsim3" --ratio 0.0,5.0 --draw cw=0.5:4 --show-unc --y-min 1E-9 --translate resources/translate_root.json
```

