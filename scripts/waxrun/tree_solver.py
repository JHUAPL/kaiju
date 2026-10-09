"""Check the user's CAM tree against the pins in config/tree_pins.json.

The coupled tree needs cmeps pinned in .gitmodules to fvitt/CMEPS waccmx-mage, or the MPMD job
hangs at ESMF init.
"""

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG_DIR = HERE / "config"
TREE_PINS_JSON = CONFIG_DIR / "tree_pins.json"
CMEPS_PIN_PATCH = CONFIG_DIR / "cmeps_pin.patch"

# needed by create_newcase and git-fleximod
MIN_PYTHON = (3, 9)

KINDS = ("coupled", "spinup")
# names the user's CAM tree, the default simulation.srcroot
WAXHOME_ENV = "WAXHOME"
CMEPS = "components/cmeps"
ESMAPP = f"{CMEPS}/cesm/driver/esmApp.F90"
MAGE_MODULE = "src/ionosphere/waccmx/mage_module.F90"

# esmApp.F90 must split MPI_COMM_WORLD on mage_module.F90's myAppId, or the job hangs at ESMF init
_SPLIT_CALL = re.compile(r"MPI_COMM_split\s*\(\s*MPI_COMM_WORLD\s*,\s*(\d+)\s*,", re.I)
_MYAPPID = re.compile(r"\bmyAppId\s*=\s*(\d+)")


def finding(level, check, detail, fix=None):
    """A case_solver.finding with an optional fix line."""
    f = {"level": level, "check": check, "detail": detail}
    if fix:
        f["fix"] = fix
    return f


def has_fail(findings):
    return any(f["level"] == "FAIL" for f in findings)


def load_pins(path=None):
    with open(path or TREE_PINS_JSON) as fh:
        return json.load(fh)


def _git(path, *args):
    """(rc, stdout) of git in path; safe.directory because the tree may belong to another user."""
    path = str(path)
    p = subprocess.run(["git", "-c", f"safe.directory={path}", "-C", path, *args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout


def git_head(srcroot):
    rc, out = _git(srcroot, "rev-parse", "HEAD")
    return out.strip() if rc == 0 else None


def git_describe(srcroot):
    rc, out = _git(srcroot, "describe", "--tags", "--always")
    return out.strip() if rc == 0 else None


def gitmodules(srcroot):
    path = Path(srcroot) / ".gitmodules"
    mods, name = {}, None
    if not path.is_file():
        return mods
    for line in path.read_text().split("\n"):
        line = line.strip()
        if line.startswith("[submodule "):
            name = line.split('"')[1]
            mods[name] = {}
        elif name and "=" in line and not line.startswith(("#", ";")):
            key, _, value = line.partition("=")
            mods[name][key.strip()] = value.strip()
    return mods


def contains_uncommented(path, marker):
    """True when marker appears on a line of path that is not a Fortran comment."""
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return False
    return any(marker in l and not l.lstrip().startswith("!") for l in text.splitlines())


def _first_group(pattern, path):
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        if line.lstrip().startswith("!"):
            continue
        m = pattern.search(line)
        if m:
            return int(m.group(1))
    return None


def split_color(path):
    """The colour of esmApp.F90's MPI_COMM_split(MPI_COMM_WORLD, ...), or None."""
    return _first_group(_SPLIT_CALL, path)


def app_id(path):
    return _first_group(_MYAPPID, path)


def colour_mismatch(esmapp, mage_module):
    """(split colour, myAppId) when both are read and disagree, else None."""
    colour, appid = split_color(esmapp), app_id(mage_module)
    if colour is None or appid is None or colour == appid:
        return None
    return colour, appid


def _clone_fix(pin, srcroot):
    return (f"git -C {srcroot} fetch {pin['remote']} {pin['branch']} && "
            f"git -C {srcroot} checkout {pin['commit']}")


def _pin_fix(srcroot):
    return (f"cd {srcroot} && git apply {CMEPS_PIN_PATCH} && {sys.executable} ./bin/git-fleximod update")


def _cmeps_fix(pin, srcroot):
    """The .gitmodules pin, then a stopgap checkout that the next fleximod update reverts."""
    ovr = pin["cmeps_override"]
    return (_pin_fix(srcroot) + "\n          or, until the next fleximod update: "
            f"cd {srcroot}/{CMEPS} && git fetch {ovr['checkout_remote']} {ovr['checkout_branch']}"
            f" && git checkout {ovr['fxtag']}")


def _fleximod_fix(srcroot):
    return f"cd {srcroot} && {sys.executable} ./bin/git-fleximod update"


def remediation(findings):
    out = []
    for f in findings:
        fix = f.get("fix")
        if fix and fix not in out:
            out.append(fix)
    return out


def report_text(spec_or_findings):
    findings = (spec_or_findings["findings"] if isinstance(spec_or_findings, dict)
                else spec_or_findings)
    lines = []
    for f in findings:
        lines.append(f"{f['level']:<4} {f['check']}: {f['detail']}")
        if f.get("fix"):
            lines.append(f"     fix: {f['fix']}")
    return "\n".join(lines)


def waxhome(environ=None):
    """$WAXHOME as a normalised absolute path, or None when unset or blank."""
    raw = str((os.environ if environ is None else environ).get(WAXHOME_ENV) or "").strip()
    if not raw:
        return None
    return os.path.normpath(os.path.abspath(os.path.expanduser(raw)))


def recipe_text(pins, srcroot="<srcroot>", spinup_srcroot="<srcroot_spinup>", kind=None):
    """The shell commands that produce a tree check_tree accepts."""
    blocks = []
    if kind in (None, "coupled"):
        c = pins["coupled"]
        blocks.append("\n".join([
            "# coupled tree (required)",
            f"git clone --branch {c['branch']} {c['remote']} {srcroot}"
            f"      # or {c['remote_https']}",
            f"cd {srcroot}",
            f"git checkout {c['commit']}",
            f"git apply {CMEPS_PIN_PATCH}"
            f"   # pin cmeps to fvitt/CMEPS@{c['cmeps_override']['fxtag'][:9]} before fleximod",
            "./bin/git-fleximod update          # with python >= %d.%d first on PATH" % MIN_PYTHON,
            "./bin/git-fleximod status          # no 's' (out-of-sync) row; 'cmeps at hash %s'"
            % c["cmeps_override"]["fxtag"][:8],
            f"export {WAXHOME_ENV}={srcroot}"
            "   # waxrun and engage GRX then use this tree by default (put it in your shell setup)",
        ]))
    if kind in (None, "spinup"):
        s = pins["spinup"]
        blocks.append("\n".join([
            "# spinup tree (only when spinup.srcroot is set; no cmeps pin)",
            f"git clone --branch {s['branch']} {s['remote']} {spinup_srcroot}",
            f"cd {spinup_srcroot}",
            f"git checkout {s['commit']}",
            "./bin/git-fleximod update && ./bin/git-fleximod status",
        ]))
    return "\n\n".join(blocks) + "\n"


def check_tree(srcroot, kind, pins=None, sourcemods="none", single_tree=False):
    """Findings for srcroot as the tree of the kind phase; waxrun refuses on any FAIL.

    With single_tree=True the spinup is checked against the coupled pin."""
    if kind not in KINDS:
        raise ValueError(f"kind={kind!r}: valid {KINDS}")
    pins = pins or load_pins()
    pin_kind = "coupled" if (single_tree and kind == "spinup") else kind
    pin = pins[pin_kind]
    srcroot = Path(srcroot)
    out = []

    create_newcase = srcroot / "cime" / "scripts" / "create_newcase"
    fleximod = srcroot / "bin" / "git-fleximod"
    missing = [str(p) for p in (create_newcase, fleximod) if not p.is_file()]
    if missing:
        # cime/ is a fleximod submodule
        if fleximod.is_file():
            detail = f"externals not fetched: {', '.join(missing)} missing (run git-fleximod update)"
            fix = _pin_fix(srcroot) if pin.get("cmeps_override") else _fleximod_fix(srcroot)
        elif (srcroot / ".git").exists():
            detail = f"not a CAM tree: missing {', '.join(missing)}"
            fix = _clone_fix(pin, srcroot)
        else:
            detail = f"{srcroot} is not a git work tree: missing {', '.join(missing)}"
            fix = recipe_text(pins, srcroot=str(srcroot), kind=pin_kind)
        out.append(finding("FAIL", "tree layout", detail, fix))
        return out
    if not os.access(str(create_newcase), os.X_OK):
        out.append(finding("FAIL", "tree layout", f"{create_newcase} is not executable",
                           f"chmod +x {create_newcase}"))
    else:
        out.append(finding("OK", "tree layout", f"{create_newcase} and {fleximod} present"))

    head = git_head(srcroot)
    describe = git_describe(srcroot)
    if head == pin["commit"]:
        out.append(finding("OK", "commit", f"{head[:9]} == {pin_kind} pin ({describe})"))
    else:
        got = f"HEAD {head[:9]} ({describe})" if head else "HEAD cannot be read"
        out.append(finding("WARN", "commit", f"{got}, not the {pin_kind} pin {pin['commit'][:9]} "
                           f"({pin['branch']})"))

    # read markers from the file that will be compiled (a SourceMod shadows its tree file by name)
    shadow = _sourcemod_shadows(kind, sourcemods)
    for rel, marker in pin.get("markers", {}).items():
        path = shadow.get(Path(rel).name) or srcroot / rel
        if not Path(path).is_file():
            out.append(finding("FAIL", "marker", f"{rel} is missing from the {pin_kind} tree",
                               _clone_fix(pin, srcroot)))
            continue
        fix = _pin_fix(srcroot) if rel.startswith(CMEPS) else _clone_fix(pin, srcroot)
        if not contains_uncommented(path, marker):
            out.append(finding("FAIL", "marker", f"{path} does not contain {marker!r} on an "
                               "uncommented line", fix))
            continue
        mismatch = colour_mismatch(path, srcroot / MAGE_MODULE) if rel == ESMAPP else None
        if mismatch:
            out.append(finding("FAIL", "marker", f"{path} splits MPI_COMM_WORLD on colour "
                               f"{mismatch[0]}, {MAGE_MODULE} has myAppId = {mismatch[1]} (the job "
                               "would hang)", fix))
        else:
            out.append(finding("OK", "marker", f"{rel} contains {marker!r}"))

    # outside .gitmodules the next fleximod update reverts the cmeps pin
    ovr = pin.get("cmeps_override")
    if ovr:
        stanza = gitmodules(srcroot).get(ovr["submodule"], {})
        got = (stanza.get("url"), stanza.get("fxtag"))
        if got == (ovr["url"], ovr["fxtag"]):
            out.append(finding("OK", "cmeps pin", f".gitmodules pins {ovr['submodule']} to "
                               f"{ovr['url']} @ {ovr['fxtag'][:9]}"))
        else:
            out.append(finding("FAIL", "cmeps pin",
                               f".gitmodules does not pin {ovr['submodule']} to {ovr['url']} @ "
                               f"{ovr['fxtag'][:9]} (has url={got[0]} fxtag={got[1]})",
                               _cmeps_fix(pin, srcroot)))

    if sourcemods not in ("none", None):
        sm = _sourcemods(kind, sourcemods)
        listed = os.path.join(sm["dir"], "tree.md5")
        if not sm["bases"]:
            out.append(finding("INFO", "sourcemod base", f"no tree.md5 in {sm['dir']}: your SourceMods are "
                               "applied without a tree check"))
        for rel, digest in sm["bases"].items():
            tree_file = srcroot / rel
            if not tree_file.is_file():
                out.append(finding("FAIL", "sourcemod base", f"{rel} is missing from the tree (listed in "
                                   f"{listed})", f"set {kind}.sourcemods=none"))
            elif _md5(tree_file) != digest:
                out.append(finding("FAIL", "sourcemod base", f"{rel} changed since {listed} was written: "
                                   "the tree moved under your SourceMods",
                                   f"re-derive your SourceMods from {tree_file} and rewrite tree.md5, or set "
                                   f"{kind}.sourcemods=none"))
            else:
                out.append(finding("OK", "sourcemod base", f"{rel} equals tree.md5"))

    if single_tree and kind == "spinup":
        out.append(finding("INFO", "single tree", "spinup on the MAGE tree "
                           f"({pins['coupled']['branch']}); a few atm_in and component defaults differ "
                           "from the reference spinup"))
    return out


def _sourcemods(kind, policy):
    # local import: case_solver owns the SourceMods sets
    import case_solver
    return case_solver.sourcemods(kind, policy)


def _md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def _sourcemod_bases(kind, policy):
    """{tree path: md5} from the set's tree.md5."""
    return _sourcemods(kind, policy)["bases"]


def _sourcemod_shadows(kind, policy):
    """{file name: the SourceMod compiled instead of the tree file of that name}."""
    if policy in ("none", None):
        return {}
    return {Path(rel).name: f for rel, f in _sourcemods(kind, policy)["files"].items()}


def tree_spec(srcroot, kind, pins=None, sourcemods="none", single_tree=False):
    """What waxrun stores in spec['tree']."""
    pins = pins or load_pins()
    findings = check_tree(srcroot, kind, pins=pins, sourcemods=sourcemods, single_tree=single_tree)
    srcroot = Path(srcroot)
    return {
        "srcroot": str(srcroot),
        "kind": kind,
        "pin_kind": "coupled" if (single_tree and kind == "spinup") else kind,
        "head": git_head(srcroot),
        "describe": git_describe(srcroot),
        "sourcemods": sourcemods,
        "single_tree": single_tree,
        "findings": findings,
        "remediation": remediation(findings),
        "ok": not has_fail(findings),
    }


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description="verify a CAM tree against config/tree_pins.json")
    p.add_argument("srcroot", nargs="?",
                   help=f"the CAM clone to check (default: ${WAXHOME_ENV}; with neither, the recipe is printed)")
    p.add_argument("--kind", choices=KINDS, default="coupled")
    p.add_argument("--sourcemods", default="none", help="none, a set name or an absolute path")
    p.add_argument("--single-tree", action="store_true")
    p.add_argument("--pins", default=None)
    p.add_argument("--recipe", action="store_true", help="print the clone recipe and exit")
    p.add_argument("--json", action="store_true", help="print the tree_spec dict")
    args = p.parse_args(argv)
    pins = load_pins(args.pins)
    args.srcroot = args.srcroot or (None if args.recipe else waxhome())
    if args.recipe or not args.srcroot:
        print(recipe_text(pins), end="")
        return 0
    try:
        spec = tree_spec(args.srcroot, args.kind, pins=pins, sourcemods=args.sourcemods,
                         single_tree=args.single_tree)
    except ValueError as e:
        raise SystemExit(f"tree_solver: {e}")
    if args.json:
        print(json.dumps(spec, indent=2))
    else:
        print(report_text(spec))
    return 0 if spec["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
