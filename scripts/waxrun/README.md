# waxrun

`waxrun` configures the two WACCM-X cases (spinup and coupled) of a GRX run on derecho. It is
normally driven by `engage.py --component grx`; the user guide is the kaiju documentation
(MAGE - With WACCM-X (GRX), build page and quick start).

## The CAM tree

Clone one CAM tree; waxrun creates and builds both cases from it.

```bash
git clone --branch waccmx_mage_cplng https://github.com/fvitt/CAM.git <srcroot>
cd <srcroot> && git checkout da3bc8d52f652b64b12d9978ba97b50fbcb750c7
git apply <KAIJUHOME>/scripts/waxrun/config/cmeps_pin.patch     # pins components/cmeps in .gitmodules
./bin/git-fleximod update && ./bin/git-fleximod status          # python3 >= 3.9 first on PATH
export WAXHOME=<srcroot>
```

The cmeps patch is required: without it the coupled job hangs at ESMF init. It has to be in
`.gitmodules` because `git-fleximod update` reverts a hand checkout.

### WAXHOME

`WAXHOME` is the default `simulation.srcroot`; `--srcroot`, the `-wo` JSON or the engage JSON win.

## The tree check

Before anything is written waxrun checks the tree and refuses on a failure, printing the command
that fixes it. A HEAD other than the pinned commit is a warning.

```bash
python tree_solver.py <srcroot>      # the report (--json: the whole spec)
python tree_solver.py --recipe       # the clone recipe
```

`--skip-tree-check` skips the check. The IC strip needs `ncks` (`module load nco`, or
`WAXRUN_NCKS=/path/to/ncks`).

## Options JSON and replay

A `-wo` JSON may hold only the options you set; the rest take their defaults. waxrun saves the full
set as `waxrun_parameters.json` in the run directory; its paths are relative to the file, so it
replays from anywhere:

```bash
engage.py --component grx -eo engage_parameters.json -mo makeitso_parameters.json \
          -wo waxrun_parameters.json
```

A replay keeps the recorded IC pick. A moved run directory is regenerated from the new place with
`--clobber --recreate-cases` added, which recreates and rebuilds both cases (CIME records them by
absolute path).
