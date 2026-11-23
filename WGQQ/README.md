# WGQQ output directory

Run the scripts from the EFT2Obs directory. Gridpacks, YODA files, merged
histograms, and scaling outputs use `WGQQ/output`. When running on lxplus
from an AFS checkout, set up `WGQQ/output` as a symlink to EOS:

```bash
ln -s /eos/user/y/yuekai/WGQQ/EFT WGQQ/output
```

The EOS directory must already exist. For other setups, an existing local
`WGQQ/output` directory also works. Git ignores `WGQQ/output`.
`storage.py` provides this shared check: Python scripts import it, and Bash
scripts run it with `python3 WGQQ/storage.py` to get the output path.

For example:

```bash
./WGQQ/make_gridpack.sh 200.0 300.0
./WGQQ/run_gridpack.sh 200.0 300.0 1  # test: 500 events
./WGQQ/run_gridpack.sh 200.0 300.0    # submit production jobs
python3 WGQQ/makePlot.py ptg200to300 ptg
```

This bin's gridpack and YODA files are stored under
`/eos/user/y/yuekai/WGQQ/EFT/ptg200to300/`. Existing files from the old layout
must be moved or copied into the matching bin directories before using them.
Rebuilding a gridpack preserves the bin's YODA files; submitting production jobs
retains the existing behavior of removing that bin's `*.yoda` files first.

MadGraph process builds, generated cards, and submission files still use their
existing locations. Batch workers need the EOS directory mounted and accessible
through the same path; an AFS symlink alone does not grant EOS access.

