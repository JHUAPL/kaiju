# waxrun SourceMods

waxrun ships no source modifications. Put your own sets here; git ignores everything in this
folder except this README (`git clean -x` deletes them, so keep a copy).

    config/sourcemods/<name>/src.cam/<file>.F90
    config/sourcemods/<name>/src.<component>/...     (any CIME SourceMods folder)
    config/sourcemods/<name>/tree.md5                (optional)

## Select a set

`spinup.sourcemods` and `coupled.sourcemods`: `none` (default), a set name, or the absolute
path of a set kept elsewhere. At `--mode INTERMEDIATE` or above the prompts are:

    Spinup SourceMods (none|<name> or a path) [none]:
    Coupled SourceMods (none|<name> or a path) [none]:

waxrun `-wo` JSON:

    {"spinup": {"sourcemods": "<name>"}, "coupled": {"sourcemods": "<name>"}}

engage `-eo` JSON:

    {"waxrun": {"spinup": {"sourcemods": "<name>"}, "coupled": {"sourcemods": "<name>"}}}

## tree.md5

Optional: the tree files your mods were derived from (`md5sum src/.../<file>.F90 > <set>/tree.md5`
in `$WAXHOME`). waxrun refuses when one of them has changed in the tree.
