"""MERRA2 nudging data for a GRX run: find, generate, validate and assemble the slot files.

CAM reads every 3-hourly slot from one Nudge_Path and aborts at the first missing one.
"""

import argparse
import getpass
import glob
import json
import math
import os
import re
import shlex
import sys
from datetime import date, datetime, timedelta
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEN_CONFIG = HERE / "config" / "nudging_gen.json"
GRIDS = HERE / "config" / "grids.json"

FILE_TEMPLATE = "MERRA2_ne30_L145.cam2.i.%Y-%m-%d-%05d.nc"     # strftime + second-of-day
TIMES_PER_DAY = 8
NUDGE_T_REMOVE = {"remove": ["Nudge_T"]}
GEN_SCRIPT = "nudge_gen.pbs"
VWIN_LINDEX, VWIN_HINDEX = 97, 146       # the reference Nudge_Vwin band


def _parse(v):
    if isinstance(v, datetime):
        return v
    if isinstance(v, date):
        return datetime(v.year, v.month, v.day)
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d", "%m-%d-%Y %H:%M:%S"):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            continue
    raise ValueError(f"unrecognized date/datetime: {v!r}")


def _as_date(v):
    return _parse(v).date()


def _expand_year(path, year):
    return str(path).replace("{YYYY}", f"{year:04d}").rstrip("/")


def _cam_template(template):
    """CAM's Nudge_File_Template spelling of FILE_TEMPLATE (%y = 4-digit year, %s = tod)."""
    return template.replace("%Y", "%y").replace("%05d", "%s")


def load_config(path=None):
    return json.loads(Path(path or GEN_CONFIG).read_text())


def grid_dims(res=None, path=None):
    """(ncol, lev) of a resolution from config/grids.json; the single entry when res is None."""
    grids = json.loads(Path(path or GRIDS).read_text())
    if res is None:
        if len(grids) != 1:
            raise ValueError(f"config/grids.json holds {len(grids)} grids; name one")
        res = next(iter(grids))
    g = grids[res]
    return int(g["ncol"]), int(g["lev"])


def window(spinup_start, stop_date):
    """(Nudge_Beg, Nudge_End): the case start day through the day after the stop."""
    return _as_date(spinup_start), _as_date(stop_date) + timedelta(days=1)


def slots(beg, end, times_per_day=TIMES_PER_DAY):
    """Every nudging instant from beg through end (inclusive dates)."""
    if 86400 % times_per_day:
        raise ValueError(f"times_per_day {times_per_day} does not divide 86400")
    step = 86400 // times_per_day
    beg, end = _as_date(beg), _as_date(end)
    out = []
    d = beg
    while d <= end:
        out += [datetime(d.year, d.month, d.day) + timedelta(seconds=i * step) for i in range(times_per_day)]
        d += timedelta(days=1)
    return out


def filename(slot, template=FILE_TEMPLATE):
    tod = slot.hour * 3600 + slot.minute * 60 + slot.second
    return slot.strftime(template.replace("%05d", f"{tod:05d}"))


_OK_CACHE = {}


def file_ok(path, ncol=None, nlev=None):
    """True when path is a complete slot file (nudge_gen.output_ok); memoized on mtime and size."""
    try:
        st = os.stat(path)
    except OSError:
        return False
    if ncol is None or nlev is None:
        d_ncol, d_lev = grid_dims()
        ncol = d_ncol if ncol is None else ncol
        nlev = d_lev if nlev is None else nlev
    key = (os.path.realpath(path), st.st_mtime_ns, st.st_size, int(ncol), int(nlev))
    if key not in _OK_CACHE:
        import nudge_gen
        _OK_CACHE[key] = bool(nudge_gen.output_ok(path, int(ncol), int(nlev)))
    return _OK_CACHE[key]


def coverage(search_dirs, slots, template=FILE_TEMPLATE, check=True, dims=None):
    """{slot: first complete file in search_dirs, or None}; '{YYYY}' in a dir is the slot's year."""
    if isinstance(search_dirs, (str, Path)):
        search_dirs = [search_dirs]
    ncol, nlev = dims if dims else (None, None)
    isdir = {}
    out = {}
    for s in slots:
        name = filename(s, template)
        out[s] = None
        for d in search_dirs:
            d = _expand_year(d, s.year)
            if d not in isdir:
                isdir[d] = os.path.isdir(d)
            if not isdir[d]:
                continue
            p = os.path.join(d, name)
            if not os.path.exists(p):
                continue
            if check and not file_ok(p, ncol, nlev):
                continue
            out[s] = p
            break
    return out


def missing(coverage_map):
    return sorted(s for s, p in coverage_map.items() if p is None)


def missing_spans(coverage_map):
    """Contiguous whole-day spans [(first, last), ...] holding any unmapped slot."""
    spans = []
    for d in sorted({s.date() for s in missing(coverage_map)}):
        if spans and d == spans[-1][1] + timedelta(days=1):
            spans[-1][1] = d
        else:
            spans.append([d, d])
    return [tuple(s) for s in spans]


def plan_generated(coverage_map, outdir, template=FILE_TEMPLATE):
    """Replace every None with the path the generator will write, so links can be made before it runs."""
    return {s: p if p is not None else os.path.join(_expand_year(outdir, s.year), filename(s, template))
            for s, p in coverage_map.items()}


def assemble(dest_dir, coverage_map, run_root):
    """Symlink every slot into dest_dir and return it as Nudge_Path (trailing /).

    Targets inside run_root are linked relatively, others absolutely."""
    gaps = missing(coverage_map)
    if gaps:
        raise ValueError(f"cannot assemble {dest_dir}: {len(gaps)} slot(s) unmapped "
                         f"({gaps[0]:%Y-%m-%d %H:%M} .. {gaps[-1]:%Y-%m-%d %H:%M}); "
                         "generate them or plan_generated() first")
    dest = str(dest_dir).rstrip("/")
    os.makedirs(dest, exist_ok=True)
    for s in sorted(coverage_map):
        target = os.path.abspath(coverage_map[s])
        link = os.path.join(dest, os.path.basename(target))
        spelled = rel_or_abs(target, dest, run_root)
        target = spelled[2:] if spelled.startswith("./") else spelled
        if os.path.islink(link):
            if os.readlink(link) == target:
                continue
            os.remove(link)
        elif os.path.exists(link):
            raise FileExistsError(f"{link} exists and is not a symlink")
        os.symlink(target, link)
    return dest + "/"


def namelist(dtime_s, beg, end, enabled, path=None, strip_nudge_t=None,
             times_per_day=TIMES_PER_DAY, template=FILE_TEMPLATE):
    """The user_nl_cam patch: Nudge_* plus Model_Times_Per_day = 86400 // dtime.

    fincl2 loses 'Nudge_T' when nudging is off, since the field does not exist then."""
    dtime_s = int(dtime_s)
    if dtime_s <= 0 or 86400 % dtime_s:
        raise ValueError(f"dtime {dtime_s} s does not divide 86400; Model_Times_Per_day must be an integer")
    beg, end = _as_date(beg), _as_date(end)
    nl = {"Nudge_Model": ".true." if enabled else ".false."}
    if path is not None:
        nl["Nudge_Path"] = f"'{str(path).rstrip('/')}/'"
    elif enabled:
        raise ValueError("nudging enabled but no Nudge_Path (the assembled directory) given")
    nl.update({
        "Nudge_File_Template": f"'{_cam_template(template)}'",
        "Nudge_Times_Per_Day": times_per_day,
        "Model_Times_Per_day": 86400 // dtime_s,
        "Nudge_Beg_Year": beg.year, "Nudge_Beg_Month": beg.month, "Nudge_Beg_Day": beg.day,
        "Nudge_End_Year": end.year, "Nudge_End_Month": end.month, "Nudge_End_Day": end.day,
    })
    if strip_nudge_t is None:
        strip_nudge_t = not enabled
    if strip_nudge_t:
        nl["fincl2"] = dict(NUDGE_T_REMOVE)
    return nl


def split_fincl(nl):
    """(nl_patches, fincl_edits) from a namelist() dict."""
    return ({k: v for k, v in nl.items() if not isinstance(v, dict)},
            {k: v for k, v in nl.items() if isinstance(v, dict)})


def _res(cfg):
    return cfg["resources"]


def gen_workers(workers=None, cfg=None):
    """The nudge_gen.py worker count, capped at the node's cores."""
    cfg = cfg or load_config()
    n = int(workers if workers not in (None, "", "auto") else _res(cfg)["workers"])
    if n < 1:
        raise ValueError(f"nudging.gen_workers must be >= 1, got {n}")
    cap = int(cfg["pbs"].get("max_ncpus") or n)
    if n > cap:
        raise ValueError(f"nudging.gen_workers={n} exceeds the node limit of {cap} cpus")
    return n


def gen_mem_gb(workers, cfg=None):
    """Job memory in whole GB (measured 1.6 GB RSS per worker)."""
    cfg = cfg or load_config()
    r = _res(cfg)
    return int(math.ceil(r["base_mem_gb"] + int(workers) * r["per_worker_mem_gb"]))


def _mem_str(gb):
    return f"{int(gb)}GB"


def gen_walltime_s(num_days, cfg=None, workers=None, times_per_day=TIMES_PER_DAY, weights=False):
    """Walltime in seconds for num_days (measured 11 s per slot on one worker, 13.5 s on four)."""
    cfg = cfg or load_config()
    r = _res(cfg)
    w = gen_workers(workers, cfg)
    batches = int(math.ceil(int(num_days) * int(times_per_day) / w))
    secs = r["fixed_s"] + batches * r["per_slot_s"] * r["margin"] + (r["weights_gen_s"] if weights else 0)
    return int(math.ceil(secs / r["round_up_s"]) * r["round_up_s"])


def _hms(secs):
    return f"{secs // 3600:02d}:{secs % 3600 // 60:02d}:{secs % 60:02d}"


def max_days(cap_s, cfg=None, workers=None, times_per_day=TIMES_PER_DAY, weights=False):
    """Most days one job may cover under a walltime cap."""
    cfg = cfg or load_config()
    r = _res(cfg)
    w = gen_workers(workers, cfg)
    head = r["fixed_s"] + (r["weights_gen_s"] if weights else 0)
    n = int((cap_s - head) / (r["per_slot_s"] * r["margin"]) * w / times_per_day) + 1
    while n > 0 and gen_walltime_s(n, cfg, w, times_per_day, weights) > cap_s:
        n -= 1
    if n < 1:
        raise ValueError(f"walltime cap {cap_s} s cannot cover even one day "
                         f"({_hms(gen_walltime_s(1, cfg, w, times_per_day, weights))} needed)")
    return n


def gen_chunks(spans, max_days, split_years=True, allow_year_crossing=False):
    """[(ref_date, num_days), ...]: each span cut at Dec 31 and into pieces of at most max_days.

    The Dec 31 cut is required: nudge_gen keeps the NCL's year*365.25 time arithmetic, which is
    wrong past Jan 1 of the ref date's year."""
    out = []
    for beg, end in spans:
        beg, end = _as_date(beg), _as_date(end)
        if end < beg:
            raise ValueError(f"span {beg} .. {end} is reversed")
        if not split_years and beg.year != end.year and not allow_year_crossing:
            raise ValueError(f"span {beg} .. {end} crosses a year boundary and split_years=False: "
                             "nudge_gen reproduces the NCL's year*365.25 time bug, so one chunk "
                             "must stay inside one calendar year (pass allow_year_crossing=True "
                             "only if you know the time/ndcur artifact is acceptable)")
        d = beg
        while d <= end:
            stop = min(end, d + timedelta(days=max_days - 1))
            if split_years:
                stop = min(stop, date(d.year, 12, 31))
            out.append((d, (stop - d).days + 1))
            d = stop + timedelta(days=1)
    return out


def rda_inputs(spans, rda_input, cfg=None):
    """{date: path or None} of the daily MERRA2 f09 inputs nudge_gen.py reads."""
    cfg = cfg or load_config()
    out = {}
    for beg, end in spans:
        d = _as_date(beg)
        while d <= _as_date(end):
            p = os.path.join(rda_input, cfg["rda_template"].format(YYYY=f"{d.year:04d}", YYYYMMDD=d.strftime("%Y%m%d")))
            out[d] = p if os.path.exists(p) else None
            d += timedelta(days=1)
    return out


def weights_cached(tmpdir):
    """True when tmpdir already holds an ESMF weight file."""
    return bool(glob.glob(os.path.join(str(tmpdir), "WgtGrid_*_conserve_none.nc")))


def gen_header(jobname, account, queue, walltime_s, workers, mem, job_priority=None):
    """The #PBS block, all values literal."""
    return [f"#PBS -N {jobname}",
            "#PBS -j oe",
            f"#PBS -q {queue}",
            f"#PBS -A {account}",
            f"#PBS -l walltime={_hms(walltime_s)}"] + (
           [f"#PBS -l job_priority={job_priority}"] if job_priority else []) + [
            f"#PBS -l select=1:ncpus={int(workers)}:mem={mem}"]


def gen_tmpdir(cfg=None):
    """The job's TMPDIR as a literal path.

    {USER} is expanded here because PBS exports LOGNAME, not USER, and $USER aborts under set -u."""
    cfg = cfg or load_config()
    raw = str(cfg["pbs"].get("tmpdir") or "")
    return os.path.expandvars(raw.replace("{USER}", getpass.getuser())).replace("$USER", getpass.getuser())


def rel_or_abs(path, base, run_root):
    """path relative to base when it lies inside run_root and resolves to the same file, else absolute."""
    t = os.path.abspath(str(path))
    b, r = os.path.abspath(str(base)), os.path.abspath(str(run_root))
    if os.path.commonpath([r, t]) != r:
        return t
    rel = os.path.relpath(t, b)
    if os.path.realpath(os.path.join(b, rel)) != os.path.realpath(t):
        return t
    return rel if rel in (".", "..") or rel.startswith(("../", "./")) else "./" + rel


def gen_marker(script):
    """<script>.OK, touched after validation; the DAG gates on it because PBS history ages out."""
    return str(script) + ".OK"


def gen_wrapper(values, cfg=None):
    """The text of one PBS wrapper around nudge_gen.py; outdir and tmpdir are relative to gen_dir."""
    cfg = cfg or load_config()
    q = cfg["pbs"]
    inp = cfg["inputs"]
    v = values
    wdir = "${KAIJUHOME}/scripts/waxrun"
    gen = os.path.join(wdir, cfg["script"])
    args = [("--ref-date", v["ref_date"]), ("--num-days", v["num_days"]), ("--outdir", v["outdir"]),
            ("--inputdir", v.get("inputdir") or inp["inputdir"]), ("--grid-info", inp["grid_info"]),
            ("--topo", inp["topo"]), ("--dst-scrip", inp["dst_scrip"]), ("--weights", "auto"),
            ("--tmpdir", v["tmpdir"]), ("--times-per-day", v["times_per_day"]),
            ("--case", inp["case"]), ("--workers", v["workers"])]
    if inp.get("esmf_rwg"):
        args.append(("--esmf-bin", inp["esmf_rwg"]))
    L = ["#!/bin/bash"]
    L += gen_header(v["jobname"], v["account"], v["queue"], v["walltime_s"], v["workers"], v["mem"],
                    v.get("job_priority"))
    L += [f"# Nudging generator, chunk {v['chunk']} of {v['chunks']}: {v['num_days']} days from "
          f"{v['ref_date']} ({v['nslots']} slots) -> {v['outdir']}",
          "set -euo pipefail",
          ""]
    mods = list(q.get("modules") or [])
    if mods:
        L += ["# environment"]
        L += [f"module load {m}" for m in mods]
    else:
        L += ["# environment (no modules: the interpreter and --esmf-bin are pinned)"]
    L += [f'export TMPDIR="${{TMPDIR:-{v["tmpdir_default"]}}}"',
          'mkdir -p "$TMPDIR"']
    L += [f"export KAIJUHOME={v['kaijuhome']}"]
    L += [f"cd {v['gen_dir']}",
          f"rm -f {v['marker']}",
          "",
          "# generate (complete files are skipped)"]
    L += [f"{v['python']} {gen} \\"]
    L += [f"    {flag} {val} \\" for flag, val in args]
    L += ["    --strict-esmf \\", "    --brief", ""]
    L += ["# validate the output (a failure holds the spinup)",
          f"{v['python']} {os.path.join(wdir, 'nudge_solver.py')} validate "
          f"--outdir {v['outdir']}/ --ref-date {v['ref_date']} --num-days {v['num_days']} "
          f"--times-per-day {v['times_per_day']} --template {shlex.quote(str(v['template']))} "
          f"--lindex {v['lindex']} --hindex {v['hindex']} || exit 1",
          "",
          "# success marker",
          f"touch {v['marker']}",
          ""]
    return "\n".join(L)


def gen_job(spans, outdir, account, jobname, gen_dir, rda_input=None, cfg=None, times_per_day=None,
            workers=None, template=FILE_TEMPLATE, lindex=VWIN_LINDEX, hindex=VWIN_HINDEX, *, run_root,
            job_priority=None):
    """Write one PBS wrapper around nudge_gen.py per chunk of the missing spans; returns their paths.

    Several chunks each get their own TMP/NN so they can run concurrently. outdir may carry '{YYYY}'."""
    cfg = cfg or load_config()
    q = cfg["pbs"]
    times_per_day = int(times_per_day or cfg["inputs"]["times_per_day"])
    nworkers = gen_workers(workers, cfg)
    mem_gb = gen_mem_gb(nworkers, cfg)
    cap_gb = float(q.get("max_mem_gb") or 0)
    if cap_gb and mem_gb > cap_gb:
        raise ValueError(f"nudging.gen_workers={nworkers} needs {mem_gb}GB, above the node limit of "
                         f"{int(cap_gb)}GB")
    mem = _mem_str(mem_gb)
    cap_s = int(q["max_walltime_s"])
    chunks = gen_chunks(spans, max_days(cap_s, cfg, nworkers, times_per_day, weights=True))
    if not chunks:
        raise ValueError("gen_job: no spans to generate")
    if rda_input:
        gaps = [d for d, p in rda_inputs(spans, rda_input, cfg).items() if p is None]
        if gaps:
            raise FileNotFoundError(f"MERRA2 input missing under {rda_input} for {len(gaps)} day(s): "
                                    f"{gaps[0]} .. {gaps[-1]}")
    gen_dir = Path(gen_dir)
    gen_dir.mkdir(parents=True, exist_ok=True)
    multi = len(chunks) > 1
    paths = []
    for k, (ref, ndays) in enumerate(chunks, start=1):
        stem = Path(GEN_SCRIPT).stem
        script = gen_dir / (f"{stem}_{k:02d}.pbs" if multi else GEN_SCRIPT)
        tmpdir = gen_dir / "TMP" / f"{k:02d}" if multi else gen_dir / "TMP"
        tmpdir.mkdir(parents=True, exist_ok=True)
        values = {
            "jobname": f"{jobname}_{k:02d}" if multi else jobname,
            "account": account, "queue": q["queue"], "mem": mem, "workers": nworkers,
            "job_priority": job_priority if q.get("job_priority", True) else None,
            "walltime_s": gen_walltime_s(ndays, cfg, nworkers, times_per_day,
                                         weights=not weights_cached(tmpdir)),
            "ref_date": ref.strftime("%Y%m%d"), "num_days": ndays,
            "nslots": ndays * times_per_day, "chunk": f"{k:02d}", "chunks": f"{len(chunks):02d}",
            "tmpdir": rel_or_abs(tmpdir, gen_dir, run_root),
            "outdir": rel_or_abs(_expand_year(outdir, ref.year), gen_dir, run_root),
            "gen_dir": str(gen_dir),
            "times_per_day": times_per_day, "python": sys.executable,
            "inputdir": rda_input, "marker": os.path.basename(gen_marker(script)),
            "tmpdir_default": gen_tmpdir(cfg),
            "template": template, "lindex": int(lindex), "hindex": int(hindex),
            "kaijuhome": str(HERE.parent.parent),
        }
        with open(script, "w", newline="\n") as f:
            f.write(gen_wrapper(values, cfg))
        os.chmod(script, 0o755)
        paths.append(str(script))
    return paths


def submit_command(script, base, cfg=None):
    """The shell line, run from base, that cds to the wrapper's directory and submits it."""
    cfg = cfg or load_config()
    script = Path(script)
    sub = cfg["pbs"]["submit"]
    return f"cd {shlex.quote(rel_or_abs(script.parent, base, base))} && {sub} ./{shlex.quote(script.name)}"


def validate_generated(paths, dims=None, lindex=VWIN_LINDEX, hindex=VWIN_HINDEX,
                       ranges=None, fields=("U", "V", "T", "Q", "PS")):
    """Findings (empty = ok) for generated files; never raises.

    Every level must be finite (CAM multiplies all levels by the window weight, and NaN*0 = NaN);
    ranges are checked only inside the nudging window lindex..hindex."""
    import numpy as np
    from netCDF4 import Dataset

    dims = dims or dict(zip(("ncol", "lev"), grid_dims()))
    ranges = ranges or {"U": (-300.0, 300.0), "V": (-300.0, 300.0), "T": (100.0, 400.0),
                        "PS": (30000.0, 110000.0)}
    ncol, lev = dims.get("ncol"), dims.get("lev")
    findings = []
    for p in paths:
        if not os.path.exists(p):
            findings.append(f"{p}: missing")
            continue
        if ncol and lev:
            floor = 4 * ncol * (4 * lev + 2)
            size = os.path.getsize(p)
            if size < floor:
                findings.append(f"{p}: truncated ({size} bytes, at least {floor} expected)")
                continue
        try:
            with Dataset(p) as ds:
                for name, n in dims.items():
                    if name not in ds.dimensions or len(ds.dimensions[name]) != n:
                        found = len(ds.dimensions[name]) if name in ds.dimensions else None
                        findings.append(f"{p}: dim {name}={found}, expected {n}")
                nlev = len(ds.dimensions["lev"]) if "lev" in ds.dimensions else 0
                for v in fields:
                    if v not in ds.variables:
                        findings.append(f"{p}: variable {v} missing")
                        continue
                    var = ds.variables[v]
                    var.set_auto_mask(False)        # the raw values CAM reads
                    full = np.asarray(var[0], dtype="f8")
                    finite = np.isfinite(full)
                    if not np.all(finite):
                        if "lev" in var.dimensions:
                            ax = tuple(i for i in range(full.ndim) if i != var.dimensions[1:].index("lev"))
                            bad = np.nonzero(~np.all(finite, axis=ax))[0] + 1
                            where = f"level(s) {bad[0]}..{bad[-1]} ({len(bad)} of {nlev})"
                        else:
                            where = "the field"
                        findings.append(f"{p}: {v} has non-finite values in {where} (CAM multiplies every "
                                        f"level, NaN*0 = NaN)")
                        continue
                    data = full[lindex - 1:min(hindex, nlev)] if "lev" in var.dimensions else full
                    if v in ranges:
                        lo, hi = ranges[v]
                        if data.min() < lo or data.max() > hi:
                            findings.append(f"{p}: {v} range [{data.min():.4g}, {data.max():.4g}] outside [{lo}, {hi}]")
        except Exception as e:
            findings.append(f"{p}: unreadable ({type(e).__name__}: {e})")
    return findings


def vertical_mismatch(nudge_path, ic_path, names=("hyai", "hybi"), atol=1e-12):
    """Findings when a nudge file's hybrid coefficients differ from the IC's.

    Nudging applies file level k to model level k without interpolation."""
    import numpy as np
    from netCDF4 import Dataset

    out = []
    try:
        with Dataset(nudge_path) as dn, Dataset(ic_path) as di:
            for v in names:
                if v not in dn.variables or v not in di.variables:
                    out.append(f"{v} missing from {nudge_path if v not in dn.variables else ic_path}")
                    continue
                a = np.asarray(dn.variables[v][:], dtype="f8").ravel()
                b = np.asarray(di.variables[v][:], dtype="f8").ravel()
                if a.shape != b.shape:
                    out.append(f"{v}: {a.size} values in the nudge file, {b.size} in the IC")
                elif not np.allclose(a, b, rtol=0.0, atol=atol):
                    out.append(f"{v}: max|diff| = {np.max(np.abs(a - b)):.3g}")
    except OSError as e:
        out.append(f"unreadable ({e})")
    return out


def expected_paths(ref_date, num_days, outdir, times_per_day=TIMES_PER_DAY, template=FILE_TEMPLATE):
    """The files one generator chunk must produce."""
    ref = _as_date(datetime.strptime(str(ref_date), "%Y%m%d")) if re.fullmatch(r"\d{8}", str(ref_date)) \
        else _as_date(ref_date)
    last = ref + timedelta(days=int(num_days) - 1)
    return [os.path.join(str(outdir).rstrip("/"), filename(s, template)) for s in slots(ref, last, times_per_day)]


def main(argv=None):
    p = argparse.ArgumentParser(description="waxrun nudge_solver CLI")
    sub = p.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="check the files one generator chunk wrote (exit 1 on any finding)")
    v.add_argument("--outdir", required=True)
    v.add_argument("--ref-date", required=True, help="YYYYMMDD (the chunk's REF_DATE)")
    v.add_argument("--num-days", required=True, type=int)
    v.add_argument("--times-per-day", type=int, default=TIMES_PER_DAY)
    v.add_argument("--template", default=FILE_TEMPLATE)
    v.add_argument("--lindex", type=int, default=VWIN_LINDEX, help="Nudge_Vwin_Lindex of the case")
    v.add_argument("--hindex", type=int, default=VWIN_HINDEX, help="Nudge_Vwin_Hindex of the case")
    args = p.parse_args(argv)
    paths = expected_paths(args.ref_date, args.num_days, args.outdir, args.times_per_day, args.template)
    findings = validate_generated(paths, lindex=args.lindex, hindex=args.hindex)
    for f in findings:
        print(f"VALIDATE: {f}")
    print(f"nudge_solver validate: {len(paths)} file(s) under {args.outdir}: "
          f"{'OK' if not findings else f'{len(findings)} finding(s)'}")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
