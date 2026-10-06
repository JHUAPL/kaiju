"""WACCM-X initial-condition (cam.i) selection and validation for waxrun.

The IC is a native restart re-dated through RUN_STARTDATE/START_TOD; the only change ever made to
the file is stripping physics-grid variables CAM cannot read back.
"""

import json
import os
import re
import shutil
import subprocess
import sys
from calendar import isleap
from collections import namedtuple
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4

HERE = Path(__file__).resolve().parent
GRIDS_JSON = HERE / "config" / "grids.json"

# 3-D variables CAM reads back from the physics grid (ncol); any other on ncol must be stripped
PHYS_OK = frozenset([
    "HO2", "N2D", "N2p", "NOp", "Np", "O1D", "O2_1D", "O2_1S", "O2p", "OH",
    "Op", "Op2D", "Op2P", "OpTM1", "TElec", "TIon", "UI", "VI", "WI", "e",
])

STRIPPED_SUFFIX = "_noBrCl"

ICCandidate = namedtuple("ICCandidate", "path month day tod year")
ICInfo = namedtuple("ICInfo", "path dims attrs var_dims date datesec varset")

_STAMP_RE = re.compile(r"(\d{4})-(\d{2})-(\d{2})-(\d{5})")
_CAMI_RE = re.compile(r"\.cam\.i\.(\d{4})-(\d{2})-(\d{2})-(\d{5})")


def _warn(s):
    print(s, file=sys.stderr)


def _dt(v):
    if isinstance(v, datetime):
        return v
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%m-%d-%Y %H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            continue
    raise ValueError(f"unrecognized datetime: {v!r}")


def parse_stamp(name):
    """(year, month, day, tod) from a YYYY-MM-DD-SSSSS stamp in name, or None."""
    m = _CAMI_RE.search(name) or _STAMP_RE.search(name)
    if not m:
        return None
    y, mo, d, tod = (int(g) for g in m.groups())
    return y, mo, d, tod


def candidate(path, stamp=None):
    """ICCandidate for one cam.i file, stamped from its name or else its date/datesec."""
    path = str(path)
    stamp = stamp or parse_stamp(os.path.basename(path))
    if stamp is None:
        info = inspect_ic(path)
        y, mo, d = info.date // 10000, (info.date // 100) % 100, info.date % 100
        stamp = (y, mo, d, info.datesec)
    y, mo, d, tod = stamp
    return ICCandidate(path, mo, d, tod, y)


def _cami_files(d):
    return sorted(p for p in Path(d).iterdir() if p.is_file() and _CAMI_RE.search(p.name))


def enumerate_candidates(source):
    """All cam.i restarts in a file, a dir, a CESM rest/ dir or an archive root, in time order."""
    src = Path(source)
    if src.is_file():
        return [candidate(src)]
    if not src.is_dir():
        raise FileNotFoundError(f"ic.start_file_dir not found: {source}")
    found = _cami_files(src)
    if not found:
        subdirs = sorted(p for p in src.iterdir() if p.is_dir() and _STAMP_RE.fullmatch(p.name))
        for sd in subdirs:
            found.extend(_cami_files(sd))
    if not found and (src / "rest").is_dir():
        return enumerate_candidates(src / "rest")
    if not found:
        raise FileNotFoundError(f"no cam.i files under ic.start_file_dir {source} "
                                f"(expected a cam.i file, a CESM rest/ dir or a dir of cam.i files)")
    cands = [candidate(p) for p in found]
    return sorted(cands, key=lambda c: (c.year, c.month, c.day, c.tod))


def _redated(c, coupling_start):
    """(datetime, clamped): c moved to the latest year that puts it before coupling_start.

    02-29 becomes 02-28 in a non-leap year (clamped=True)."""
    def mk(y):
        mo, d = c.month, c.day
        clamped = (mo, d) == (2, 29) and not isleap(y)
        if clamped:
            mo, d = 2, 28
        return datetime(y, mo, d) + timedelta(seconds=c.tod), clamped
    r, clamped = mk(coupling_start.year)
    if r >= coupling_start:
        r, clamped = mk(coupling_start.year - 1)
    return r, clamped


def _hours(coupling_start, redated):
    return (coupling_start - redated).total_seconds() / 3600.0


def stamp_of(c):
    return f"{c.year:04d}-{c.month:02d}-{c.day:02d}-{c.tod:05d}"


def stamp_iso(c):
    return (f"{c.year:04d}-{c.month:02d}-{c.day:02d}T{c.tod // 3600:02d}:{c.tod % 3600 // 60:02d}:"
            f"{c.tod % 60:02d}")


def series_stamps(candidates):
    """The distinct restart times of the candidates, in time order."""
    return sorted({stamp_iso(c) for c in candidates})


def candidate_table(candidates, coupling_start, target_start, min_h, max_h):
    """The refusal text listing each candidate's re-dated time and spinup length."""
    L = [f"No restart in ic.start_file_dir gives a spinup of {min_h:g}..{max_h:g} h before the coupled "
         f"start {coupling_start:%Y-%m-%d %H:%M:%S} (target spinup start "
         f"{target_start:%Y-%m-%d %H:%M:%S}):",
         f"  {'restart':<19} {'re-dated':<19} {'spinup [h]':>11}  bounds"]
    for c in candidates:
        r, _ = _redated(c, coupling_start)
        h = _hours(coupling_start, r)
        ok = "ok" if min_h <= h <= max_h else ("too short" if h < min_h else "too long")
        L.append(f"  {stamp_iso(c):<19} {r:%Y-%m-%dT%H:%M:%S}  {h:>10.1f}  {ok}")
    L.append(f"Raise spinup.max_hours (now {max_h:g}), use an ic.start_file_dir with a restart from the right "
             f"season, or set the start file date and the spinup start date.")
    return "\n".join(L)


START_FILE_DATE_RE = re.compile(r"(?:(\d{4})-)?(\d{2})-(\d{2})(?:[T ](\d{2}):(\d{2}):(\d{2}))?")
START_FILE_DATE_FORMS = "YYYY-MM-DDThh:mm:ss, YYYY-MM-DD or MM-DD"


def matches_date(c, start_file_date):
    """True when c is a restart of start_file_date; a date without a time matches the whole day."""
    m = START_FILE_DATE_RE.fullmatch(str(start_file_date).strip())
    if not m:
        raise ValueError(f"ic.start_file_date must be {START_FILE_DATE_FORMS}, got {start_file_date!r}")
    y, mo, d, hh, mi, ss = m.groups()
    if (c.month, c.day) != (int(mo), int(d)) or (y is not None and c.year != int(y)):
        return False
    return hh is None or c.tod == int(hh) * 3600 + int(mi) * 60 + int(ss)


def series_span(candidates, limit=20):
    """The restart times of the candidates on one line, abbreviated past limit."""
    return series_span_line(series_stamps(candidates), limit)


def series_span_line(days, limit=20):
    if len(days) <= limit:
        return ", ".join(days)
    return f"{len(days)} restarts: {', '.join(days[:3])}, ... {', '.join(days[-3:])}"


def select_restart(candidates, coupling_start, target_start, min_h=48, max_h=480,
                   start_file_date=None):
    """The restart whose re-dated time is nearest target_start within min_h..max_h of spinup.

    Ties go to the earlier one. start_file_date picks that restart with no bounds."""
    coupling_start, target_start = _dt(coupling_start), _dt(target_start)
    pool = list(candidates)
    if start_file_date is not None:
        pool = [c for c in pool if matches_date(c, start_file_date)]
        if not pool:
            raise SystemExit(f"Start file date {start_file_date} is not in the IC series "
                             f"(restarts: {series_span(candidates)}).")
        times = series_stamps(pool)
        if len(times) > 1:
            raise SystemExit(f"Start file date {start_file_date} matches {len(times)} restarts "
                             f"({', '.join(times)}); give the date and time.")
        min_h, max_h = float("-inf"), float("inf")
    picks = []
    for c in pool:
        r, clamped = _redated(c, coupling_start)
        h = _hours(coupling_start, r)
        if min_h <= h <= max_h:
            picks.append((abs(r - target_start), r, c, clamped))
    if not picks:
        raise SystemExit(candidate_table(pool, coupling_start, target_start, min_h, max_h))
    picks.sort(key=lambda p: (p[0], p[1]))
    _, r, c, clamped = picks[0]
    if clamped:
        _warn(f"[waxrun] WARNING: restart {stamp_of(c)} is re-dated to {r:%Y-%m-%d}: {r.year} is not a "
              f"leap year")
    return c


def spinup_start_error(spinup_start, coupling_start, min_h, max_h, dtime=None):
    """None, or why spinup_start is outside min_h..max_h or off the spinup timestep."""
    h = _hours(_dt(coupling_start), _dt(spinup_start))
    stamp = f"{_dt(spinup_start):%Y-%m-%dT%H:%M:%S}"
    if h < min_h:
        return (f"Spinup start {stamp} gives a {h:g} h spinup, shorter than spinup.min_hours "
                f"({min_h:g} h).")
    if h > max_h:
        return (f"Spinup start {stamp} gives a {h:g} h spinup, longer than spinup.max_hours "
                f"({max_h:g} h).")
    if dtime and (h * 3600) % float(dtime):
        return f"Spinup start {stamp} is not on the {float(dtime):g} s spinup timestep."
    return None


def redated(candidate, coupling_start):
    return _redated(candidate, _dt(coupling_start))[0]


def default_spinup_start(candidate, coupling_start, min_h, max_h, dtime):
    """(start, moved): the candidate's re-dated time, moved into min_h..max_h when outside.

    A moved start is a whole number of dtime steps before the coupled start."""
    cs_ = _dt(coupling_start)
    r = redated(candidate, cs_)
    h = _hours(cs_, r)
    if min_h <= h <= max_h:
        return r, False
    dt = float(dtime)
    if h > max_h:
        steps = int(max_h * 3600 // dt)
    else:
        steps = -int(-(min_h * 3600) // dt)
    return cs_ - timedelta(seconds=steps * dt), True


def redate(candidate, coupling_start, spinup_start=None):
    """(run_startdate, start_tod, stop_seconds) of a spinup from spinup_start to coupling_start."""
    coupling_start = _dt(coupling_start)
    r = _dt(spinup_start) if spinup_start is not None else _redated(candidate, coupling_start)[0]
    secs = int(round((coupling_start - r).total_seconds()))
    return r.strftime("%Y-%m-%d"), r.hour * 3600 + r.minute * 60 + r.second, secs


def stop_n(stop_seconds):
    """(STOP_OPTION, STOP_N) in whole hours when possible, else seconds."""
    if stop_seconds % 3600 == 0:
        return "nhours", stop_seconds // 3600
    return "nseconds", stop_seconds


def _py(v):
    return v.item() if hasattr(v, "item") else v


def inspect_ic(path):
    """Header-only read of a cam.i file."""
    path = str(path)
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        dims = {k: len(v) for k, v in ds.dimensions.items()}
        attrs = {a: _py(ds.getncattr(a)) for a in ds.ncattrs()}
        var_dims = {k: tuple(v.dimensions) for k, v in ds.variables.items()}
        date = int(ds.variables["date"][0]) if "date" in ds.variables else None
        datesec = int(ds.variables["datesec"][0]) if "datesec" in ds.variables else 0
    return ICInfo(path, dims, attrs, var_dims, date, datesec, frozenset(var_dims))


def load_grids(path=GRIDS_JSON):
    return json.loads(Path(path).read_text())


def expected_grid(res, grids=None):
    """The grids.json row for res (a full CIME resolution or its atm grid)."""
    grids = grids or load_grids()
    if res in grids:
        return grids[res]
    for key, row in grids.items():
        if key.split("_")[0] == res:
            return row
    raise KeyError(f"no grid definition for res {res!r} in {GRIDS_JSON} (known: {sorted(grids)})")


def assert_grid(info, res, grids=None):
    """Raise ValueError listing every dim/attr of info that differs from grids.json for res."""
    exp = expected_grid(res, grids)
    bad = []
    for dim, n in exp.items():
        if dim == "attrs":
            continue
        got = info.dims.get(dim)
        if got != n:
            bad.append(f"dim {dim}: found {got}, expected {n}")
    for a, v in exp.get("attrs", {}).items():
        got = info.attrs.get(a)
        if got != v:
            bad.append(f"attr {a}: found {got!r}, expected {v!r}")
    if bad:
        raise ValueError(f"{info.path} is not a native {res} initial file:\n  "
                         + "\n  ".join(bad))


# what CAM reads from every ncdata; constituents are defaulted when absent
REQUIRED_IC_VARS = ("U", "V", "T", "hyai", "hybi", "hyam", "hybm")
REQUIRED_IC_ANY = (("PS", "PSDRY"),)


def assert_required(info):
    """Raise ValueError when info lacks a field CAM reads from every ncdata."""
    have = info.varset
    missing = [v for v in REQUIRED_IC_VARS if v not in have]
    missing += ["|".join(alts) for alts in REQUIRED_IC_ANY if not any(a in have for a in alts)]
    if missing:
        raise ValueError(f"{info.path} lacks {', '.join(missing)}: CAM reads PS (or PSDRY), U, V, T and "
                         f"the hybrid coefficients hyai/hybi/hyam/hybm from ncdata and aborts "
                         f"without them")


def physgrid_offenders(info, phys_ok=PHYS_OK, constituents=None):
    """3-D ncol variables to strip: those not in phys_ok, plus any transported constituent.

    CAM aborts on a constituent stored on ncol ("Could not find BR field") but defaults a missing one."""
    phys_ok = set(phys_ok if phys_ok is not None else PHYS_OK)
    on_ncol = {v for v, dims in info.var_dims.items() if len(dims) == 3 and "ncol" in dims}
    out = on_ncol - phys_ok
    if constituents is not None:
        out |= on_ncol & set(constituents)
    return out


_CHEM_OPT = re.compile(r"(?:^|\s)-chem\s+(\S+)")
_MECH_CACHE = {}


def _compset_lname(srcroot, compset):
    """The long name of a compset alias from the tree, or None."""
    if "_" in str(compset) and "%" in str(compset):
        return str(compset)
    path = os.path.join(str(srcroot), "cime_config", "config_compsets.xml")
    try:
        text = Path(path).read_text()
    except OSError:
        return None
    m = re.search(r"<alias>\s*" + re.escape(str(compset)) + r"\s*</alias>\s*<lname>\s*(\S+?)\s*</lname>", text)
    return m.group(1) if m else None


def chem_mechanism(srcroot, compset, opts_append=""):
    """The chemistry package (last -chem in CAM_CONFIG_OPTS) a compset configures, or None.

    With no -chem, -waccmx defaults to waccm_ma."""
    lname = _compset_lname(srcroot, compset)
    if lname is None:
        return None
    path = os.path.join(str(srcroot), "cime_config", "config_component.xml")
    try:
        text = Path(path).read_text()
    except OSError:
        return None
    m = re.search(r'<entry id="CAM_CONFIG_OPTS">(.*?)</entry>', text, re.S)
    if not m:
        return None
    opts = []
    for rx, val in re.findall(r'<value\s+compset="([^"]*)"\s*>([^<]*)</value>', m.group(1)):
        if re.search(rx, lname):
            opts.append(val)
    opts.append(str(opts_append or ""))
    chem = None
    for o in opts:
        for c in _CHEM_OPT.findall(o):
            chem = c.lower()
    if chem is None and any(re.search(r"(?:^|\s)-waccmx\b", o) for o in opts):
        chem = "waccm_ma"
    return chem


def transported_constituents(srcroot, mech):
    """The Solution species of pp_<mech>/chem_mech.in minus its Not-Transported list."""
    key = (os.path.realpath(str(srcroot)), mech)
    if key in _MECH_CACHE:
        return _MECH_CACHE[key]
    path = os.path.join(str(srcroot), "src", "chemistry", f"pp_{mech}", "chem_mech.in")
    text = re.sub(r"\*.*", "", Path(path).read_text())

    def section(name):
        m = re.search(rf"^\s*{name}\s*$(.*?)^\s*End\s+{name}\s*$", text, re.M | re.S | re.I)
        if not m:
            return []
        items = [x.strip() for x in m.group(1).replace("\n", ",").split(",")]
        return [re.split(r"\s|->|=", x)[0] for x in items if x]

    solution = set(section("Solution"))
    if not solution:
        raise ValueError(f"{path}: no Solution species section")
    out = frozenset(solution - set(section("Not-Transported")))
    _MECH_CACHE[key] = out
    return out


def constituents_for(srcroot, compset, opts_append=""):
    """(transported constituents or None, source or reason it could not be derived)."""
    try:
        mech = chem_mechanism(srcroot, compset, opts_append)
    except re.error as e:
        return None, f"config_component.xml CAM_CONFIG_OPTS regex: {e}"
    if not mech:
        return None, f"the chemistry package of compset {compset!r} could not be read from {srcroot}/cime_config"
    try:
        return transported_constituents(srcroot, mech), f"pp_{mech}/chem_mech.in"
    except (OSError, ValueError) as e:
        return None, str(e)


class NcksNotFound(FileNotFoundError):
    """No usable ncks (NCO) executable."""


NCKS_MISSING = "ncks (module nco) is not on PATH: module load nco"


def resolve_ncks(ncks=None):
    """The ncks executable: the argument, else $WAXRUN_NCKS, else ncks on PATH."""
    for src, exe in (("ncks=", ncks), ("WAXRUN_NCKS=", os.environ.get("WAXRUN_NCKS"))):
        if exe:
            if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
                raise NcksNotFound(f"{src}{exe} is not an executable ncks")
            return exe
    exe = shutil.which("ncks")
    if not exe:
        raise NcksNotFound(NCKS_MISSING)
    return exe


def require_ncks():
    """Exit with one ERROR line unless ncks resolves; run before a GRX run writes anything."""
    try:
        return resolve_ncks()
    except NcksNotFound as e:
        raise SystemExit(f"ERROR: {e}") from None


def _ncks_cmd(args, ncks=None):
    return [resolve_ncks(ncks)] + list(args)


def stripped_name(path):
    """<stem>_noBrCl.nc for a cam.i path, without doubling the suffix."""
    stem = os.path.basename(str(path))
    stem = stem[:-3] if stem.endswith(".nc") else stem
    if stem.endswith(STRIPPED_SUFFIX):
        stem = stem[:-len(STRIPPED_SUFFIX)]
    return stem + STRIPPED_SUFFIX + ".nc"


def strip_offenders(path, offenders, shared_dir, ncks=None):
    """Strip offenders from path into shared_dir with ncks and return the new path.

    An existing output with the right variable set is reused."""
    path = str(path)
    offenders = sorted(offenders)
    if not offenders:
        raise ValueError("strip_offenders called with no offenders")
    out = os.path.join(str(shared_dir), stripped_name(path))
    if os.path.abspath(out) == os.path.abspath(path):
        raise ValueError(f"refusing to strip {path} onto itself; use another ic.shared_dir")
    want = inspect_ic(path).varset - set(offenders)
    if os.path.isfile(out) and inspect_ic(out).varset == want:
        _warn(f"[waxrun] reusing stripped IC {out}")
        return out
    os.makedirs(str(shared_dir), exist_ok=True)
    cmd = _ncks_cmd(["-O", "-x", "-v", ",".join(offenders), path, out], ncks)
    print("+ " + " ".join(cmd))
    r = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if r.stderr:
        sys.stderr.write(r.stderr)
    if r.returncode:
        raise RuntimeError(f"ncks failed (rc={r.returncode}) stripping {offenders} from {path}: "
                           f"{(r.stderr or '').strip()[-2000:] or '(no stderr)'}")
    got = inspect_ic(out).varset
    if got != want:
        raise ValueError(f"stripped IC {out} has the wrong variable set: "
                         f"unexpected {sorted(got - want)}, missing {sorted(want - got)}")
    return out


def validate_ic(path, res, phys_ok=PHYS_OK, shared_dir=None, reference=None, ncks=None,
                grids=None, offenders=None):
    """Check the grid, strip any offenders into shared_dir, and return the IC pick dict.

    A clean file is used as-is. offenders defaults to physgrid_offenders(info, phys_ok)."""
    path = str(path)
    info = inspect_ic(path)
    assert_grid(info, res, grids)
    if offenders is None:
        offenders = physgrid_offenders(info, phys_ok)
    else:
        offenders = set(offenders)
        absent = offenders - info.varset
        if absent:
            raise ValueError(f"{path} lacks the offenders {sorted(absent)} the spec planned to strip")
    stripped = None
    if offenders:
        if shared_dir is None:
            raise ValueError(f"{path} carries physgrid offenders {sorted(offenders)}; "
                             f"ic.shared_dir is required to strip them")
        stripped = strip_offenders(path, offenders, shared_dir, ncks)
        print(f"[waxrun] IC {path}: stripped {sorted(offenders)} -> {stripped}")
    else:
        print(f"[waxrun] IC {path}: no physgrid offenders, used as-is ({len(info.varset)} vars)")
    if reference:
        try:
            ref = inspect_ic(reference).varset
            final = info.varset - offenders
            print(f"[waxrun] IC vs reference {reference}: +{sorted(final - ref)} "
                  f"-{sorted(ref - final)} (informational)")
        except OSError as e:
            _warn(f"[waxrun] reference IC not readable ({e}); skipping the varset diff")
    return {
        "path": stripped or path,
        "raw": path,
        "stripped": stripped,
        "offenders": sorted(offenders),
        "grid": {k: info.dims.get(k) for k in ("ncol", "ncol_d", "lev", "ilev")},
        "nvars": len(info.varset) - len(offenders),
        "date": info.date,
        "datesec": info.datesec,
    }
