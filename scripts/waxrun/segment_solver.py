"""Segment math for GRX (GAMERA + RAIJU + WACCM-X) runs, plus the shared cadence and walltime rules.

Segment 0 is the magnetosphere spinup, 0 < ii <= ic the coupling delay, ii > ic coupled.
CESM's STOP_N and REST_N must equal stephours, or the rpointer stamps the jobs relink to won't exist.
"""

import math
from datetime import datetime, timedelta

# standalone-waxrun defaults
DEFAULTS = {
    "runlabel": "msphere",
    "res": "quad",              # quad | oct | hex
    "stephours": 4,             # must equal CESM STOP_N
    "delaycouplinghrs": 4,      # coupling delay; coupled start = t0 + this
    "nhelpers": 2,
    "Ca": 10.0,
    "dtout": 60,                # seconds
    "ewmtype": "BW",            # RAIJU electron-loss model: BW | SS
    "walltime": "12:00:00",
    "queue": "main",
    "nwnodes": 16, "nwcpus": 128,
    "ngcpus": 2, "ngcpernode": 64,
    "bcwifn": "bcwind.h5",
}

# GAMERA grid file and node count per resolution
_RES_GRID = {
    "quad": ("lfmQ.h5", 8), "quadoct": ("lfmQ.h5", 8),
    "oct": ("lfmO.h5", 36), "hex": ("lfmH.h5", 6),
}


def _cfg(overrides=None):
    c = dict(DEFAULTS)
    if overrides:
        c.update({k: v for k, v in overrides.items() if v is not None})
    c["gridfn"], c["ngnodes"] = _RES_GRID.get(c["res"], _RES_GRID["quad"])
    c["nwaccmx"] = c["nwnodes"] * c["nwcpus"]
    c["nvcpus"] = c["ngnodes"] * c["ngcpus"] + c["nhelpers"] + 1    # voltron ranks
    return c


# the only user-facing run date/time spelling, shared with engage, makeitso and tiegcmrun
RUN_DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S"


def parse_run_datetime(value, field="date"):
    """Parse a user-given YYYY-MM-DDThh:mm:ss; anything else raises ValueError naming field."""
    if isinstance(value, datetime):
        return value
    try:
        return datetime.strptime(str(value).strip(), RUN_DATETIME_FORMAT)
    except (TypeError, ValueError):
        raise ValueError(f"{field} must be given as yyyy-mm-ddThh:mm:ss "
                         f"(e.g. 2013-03-17T00:00:00); got {value!r}") from None


def _parse_dt(v):
    if isinstance(v, datetime):
        return v
    for fmt in ("%m-%d-%Y %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            continue
    raise ValueError(f"unrecognized datetime: {v!r}")


# dtCouple (s) is the timestep of both coupled models, so spin-up, segment boundaries and output
# must all fall on multiples of it.
DTCOUPLE_DEFAULT = 5.0
# TIE-GCM needs a step dividing 60 s; CAM needs dtime dividing 3600 and 86400
VALID_DTCOUPLE = tuple(d for d in range(1, 61) if 60 % d == 0)
SPINUP_DTIME_DEFAULT = 300
_DTCOUPLE_WARNED = set()


def _seconds(value):
    """float(value), or None for a missing, non-numeric or non-finite value."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _fmt(value):
    v = float(value)
    return f"{v:.0f}" if v.is_integer() else f"{v:g}"


def _off(value, step):
    """True when value is not a whole multiple of step."""
    return abs(value - step * round(value / step)) > 1e-6


def _second_of_day(t):
    return t.hour * 3600 + t.minute * 60 + t.second + t.microsecond / 1e6


def dtcouple_error(dtCouple, name="voltron.coupling.dtCouple"):
    """None, or why dtCouple is not a whole number of seconds dividing 60."""
    v = _seconds(dtCouple)
    if v is not None and v.is_integer() and int(v) in VALID_DTCOUPLE:
        return None
    return (f"{name}={dtCouple} s is not a whole number of seconds that divides 60 (valid: "
            f"{', '.join(map(str, VALID_DTCOUPLE))}): both models step on dtCouple, which must "
            f"divide 60 s for TIE-GCM and 3600 s for WACCM-X")


def cadence_error(dtCouple, tSpin=None, gr_warm_up_time=None, segment_duration=None,
                  coupled_start=None, dtOut=None, window=None):
    """None, or the first coupled-cadence rule the given values break (None arguments are skipped).

    All times (s) must be multiples of dtCouple; with dtOut (GRT), segment and window of dtOut too."""
    err = dtcouple_error(dtCouple)
    if err:
        return err
    dtc = float(dtCouple)
    d = _fmt(dtc)
    hang = ("so voltron and the GCM would not exchange at the same instants and the coupled job "
            "would hang until walltime")
    v = _seconds(tSpin)
    if tSpin is not None and (v is None or _off(v, dtc)):
        return (f"voltron.spinup.tSpin={tSpin} s is not a multiple of dtCouple={d} s: voltron's "
                f"coupling clock starts at -tSpin, {hang}")
    v = _seconds(gr_warm_up_time)
    if gr_warm_up_time is not None and (v is None or _off(v, dtc)):
        return (f"coupling.gr_warm_up_time={gr_warm_up_time} s is not a multiple of dtCouple={d} "
                f"s: the coupled start (the GR start + gr_warm_up_time) falls between two of "
                f"voltron's exchanges, {hang}")
    v = _seconds(segment_duration)
    if segment_duration is not None and (v is None or _off(v, dtc)):
        return (f"simulation.segment_duration={segment_duration} s is not a multiple of "
                f"dtCouple={d} s: every coupled segment would end between two exchanges")
    if coupled_start is not None:
        try:
            cs = _parse_dt(coupled_start)
        except (TypeError, ValueError):
            cs = None                      # reported by the date's own validator
        if cs is not None and _off(_second_of_day(cs), dtc):
            return (f"the coupled start {cs:%Y-%m-%dT%H:%M:%S} is {_fmt(_second_of_day(cs))} s into "
                    f"its day, not a multiple of dtCouple={d} s: both models start on the dtCouple "
                    f"clock")
    v = _seconds(window)
    if window is not None and v is not None and _off(v, dtc):
        return (f"the coupled window (coupled start -> stop_date, {_fmt(v)} s) is not a multiple "
                f"of dtCouple={d} s: the last coupled segment would end between two exchanges")
    if dtOut is not None:
        out = _seconds(dtOut)
        if out is None or out <= 0 or _off(out, dtc):
            return (f"voltron.output.dtOut={dtOut} s is not a multiple of dtCouple={d} s: the "
                    f"coupled TIE-GCM writes a secondary history every dtOut, which must be a "
                    f"multiple of its step")
        seg = _seconds(segment_duration)
        if seg is not None and _off(seg, out):
            return (f"simulation.segment_duration={segment_duration} s is not a multiple of "
                    f"voltron.output.dtOut={_fmt(out)} s: each coupled TIE-GCM segment must hold "
                    f"a whole number of secondary histories")
        if v is not None and _off(v, out):
            return (f"the coupled window (coupled start -> stop_date, {_fmt(v)} s) is not a "
                    f"multiple of voltron.output.dtOut={_fmt(out)} s: the last coupled TIE-GCM "
                    f"segment would not hold a whole number of secondary histories")
    return None


def coupling_start_error(coupling_start, dtCouple=None, spinup_dtime=SPINUP_DTIME_DEFAULT,
                         name="simulation.coupling_start"):
    """None, or why the coupling start's second of day is not a multiple of spinup_dtime and dtCouple.

    The spin-up case must stop exactly there to write the coupled case's initial file."""
    try:
        cs = _parse_dt(coupling_start)
    except (TypeError, ValueError):
        return None
    sod = _second_of_day(cs)
    dt = _seconds(spinup_dtime) or SPINUP_DTIME_DEFAULT
    stamp = f"{cs:%Y-%m-%dT%H:%M:%S}"
    if _off(sod, dt):
        prev = int(sod // dt * dt)
        return (f"{name}={stamp} is {_fmt(sod)} s into its day, not a multiple of the WACCM-X "
                f"spin-up dtime ({_fmt(dt)} s): the spin-up case must stop exactly at the coupling "
                f"start (STOP_OPTION=nseconds) to write the coupled case's initial file "
                f"<case>.cam.i.{cs:%Y-%m-%d}-{int(sod):05d}.nc, and a STOP_N off its {_fmt(dt)} s "
                f"clock cannot end there; move the coupling start to a multiple of {_fmt(dt)} s "
                f"of the day (e.g. {prev // 3600:02d}:{prev % 3600 // 60:02d}:{prev % 60:02d})")
    if dtCouple is not None and dtcouple_error(dtCouple) is None and _off(sod, float(dtCouple)):
        return (f"{name}={stamp} is {_fmt(sod)} s into its day, not a multiple of dtCouple="
                f"{_fmt(dtCouple)} s: the coupled WACCM-X case starts (START_TOD) on its dtime = "
                f"dtCouple clock")
    return None


def dtcouple_warning(dtCouple):
    """A warning for a dtCouple other than the tested 5 s, else None."""
    v = _seconds(dtCouple)
    if v is None or v == DTCOUPLE_DEFAULT:
        return None
    s = _fmt(v)
    return (f"WARNING: dtCouple={s}: gamera's inner boundary keeps its built-in 5 s coupling "
            f"timescale whatever dtCouple is, and no coupled run at dtCouple={s} has ever executed")


def warn_dtcouple_once(dtCouple, emit=print):
    """Emit dtcouple_warning once per process per value; returns the text emitted, else None."""
    text = dtcouple_warning(dtCouple)
    if text is None or _seconds(dtCouple) in _DTCOUPLE_WARNED:
        return None
    _DTCOUPLE_WARNED.add(_seconds(dtCouple))
    emit(text)
    return text


def compute_segments(t0, t1, stephours, delaycouplinghrs):
    """(ic, iu, maxhours): segments are ii in range(0, iu), ic of them before coupling.

    The coupled window must be a whole number of stephours segments, else ValueError."""
    t0, t1 = _parse_dt(t0), _parse_dt(t1)
    total_s = (t1 - t0).total_seconds()
    maxhours = int(total_s // 3600)
    ic = int(math.ceil(delaycouplinghrs / stephours))
    window_s = total_s - delaycouplinghrs * 3600
    span = (f"the coupled window {t0 + timedelta(hours=delaycouplinghrs):%Y-%m-%dT%H:%M:%S} -> "
            f"{t1:%Y-%m-%dT%H:%M:%S} ({window_s / 3600:g} h)")
    if window_s <= 0:
        raise ValueError(f"{span} is empty: the stop must be after t0 + the {delaycouplinghrs} h "
                         f"coupling delay")
    if window_s % (stephours * 3600):
        valid = valid_stephours(window_s / 3600, delaycouplinghrs, math.inf) if not window_s % 3600 else []
        raise ValueError(f"{span} is not a whole number of {stephours}-h segments: the last "
                         f"{window_s % (stephours * 3600) / 3600:g} h would never run; valid segment "
                         f"lengths (hours, before the walltime rule): {valid}")
    iu = ((maxhours - delaycouplinghrs) // stephours) + ic + 1
    return ic, iu, maxhours


# Walltime rule: max segment hours = floor(min(walltime, queue max) / (margin x rate)), rate in
# wall-clock hours per simulated hour. An unrecorded rate falls back to RATE_UNKNOWN and only warns.
RATE_MARGIN = {"measured": 1.5, "estimate": 2.0}
RATE_UNKNOWN = {"rate": 3.0, "basis": "estimate"}
# measured GRX rate (Q grid, ne30)
WALL_H_PER_SIM_H = 1.25
WALLTIME_MARGIN = RATE_MARGIN["measured"]


def _rate_entry(value):
    """A rates entry as {rate, basis}; a bare number is a measured rate."""
    if isinstance(value, dict):
        return {"rate": float(value["rate"]), "basis": str(value.get("basis") or "measured")}
    return {"rate": float(value), "basis": "measured"}


def _same_res(key, res):
    try:
        return float(key) == float(res)
    except (TypeError, ValueError):
        return str(key).strip().lower() == str(res).strip().lower()


def lookup_rate(table, component, grid=None, res=None):
    """{rate, basis, margin, k, key, note} for component on grid/res from a machines.yaml rates table.

    Falls back to grid/'*', then the grid's largest rate when res is unknown, then RATE_UNKNOWN."""
    entries = (((table or {}).get(component) or {}).get(str(grid)) if grid not in (None, "") else None) or {}
    hit, key, note = None, None, None
    if res not in (None, ""):
        for k_, v in entries.items():
            if k_ != "*" and _same_res(k_, res):
                hit, key = v, f"{grid}/{k_}"
                break
    if hit is None and "*" in entries:
        hit, key = entries["*"], f"{grid}/*"
    if hit is None and res in (None, "") and entries:
        k_, hit = max(entries.items(), key=lambda kv: _rate_entry(kv[1])["rate"])
        key = f"{grid}/{k_} (the largest of grid {grid}'s rates: the resolution is not known yet)"
    if hit is None:
        hit, key = RATE_UNKNOWN, "unknown"
        note = (f"NOTE: no wall-clock rate is recorded for a {component} segment on grid {grid}"
                f"{'' if res in (None, '') else f' at resolution {res}'}: the walltime rule uses "
                f"{RATE_UNKNOWN['rate']:g} h per simulated hour ({RATE_UNKNOWN['basis']}), "
                f"conservatively.")
    e = _rate_entry(hit)
    e["margin"] = float(RATE_MARGIN.get(e["basis"], RATE_MARGIN["estimate"]))
    e["k"] = e["margin"] * e["rate"]
    e["key"], e["note"] = key, note
    return e


def rate_text(info):
    return f"rate {info['rate']:g} h/h {info['basis']} x margin {info['margin']:g}"


def hms(hours):
    secs = int(round(float(hours) * 3600))
    return f"{secs // 3600:02d}:{secs % 3600 // 60:02d}:{secs % 60:02d}"


def walltime_hours(walltime):
    """Hours of an HH:MM[:SS] walltime."""
    parts = [int(p) for p in str(walltime).strip().split(":")]
    if not 2 <= len(parts) <= 3:
        raise ValueError(f"not a walltime: {walltime!r}")
    parts += [0] * (3 - len(parts))
    return parts[0] + parts[1] / 60 + parts[2] / 3600


def fit_hours(walltime, walltime_max=None):
    """(walltime hours, queue max hours or None, the smaller of the two)."""
    wall_h = walltime_hours(walltime)
    try:
        max_h = walltime_hours(walltime_max) if walltime_max else None
    except (TypeError, ValueError):
        max_h = None
    return wall_h, max_h, (wall_h if max_h is None else min(wall_h, max_h))


def max_segment_hours(wall_h, k):
    """The longest whole-hour segment wall_h hours hold at k = margin x rate."""
    return int(math.floor(float(wall_h) / float(k) + 1e-9))


def advisory(info, max_h):
    """True when the rule only warns: the rate is unknown, or even the queue max holds no 1 h segment."""
    return info.get("key") == "unknown" or (max_h is not None and max_segment_hours(max_h, info["k"]) < 1)


def walltime_error(info, walltime, walltime_max=None, segment_seconds=None, *, queue=None, hpc="",
                   config="config/machines.yaml", where="", label="", what="segment_duration",
                   remedy="shorten segment_duration", warn=None):
    """None, or why walltime cannot run a segment at the rate info (from lookup_rate).

    An advisory rate never refuses; warn(text) is called instead."""
    try:
        wall_h, max_h, fit_h = fit_hours(walltime, walltime_max)
    except (TypeError, ValueError):
        return None
    if max_h is not None and wall_h > max_h + 1e-9:
        return (f"pbs.walltime {walltime} is more than queue {queue!r} on {hpc} allows ({walltime_max}, "
                f"{config} queues): lower pbs.walltime or use another queue")
    top = max_segment_hours(fit_h, info["k"])
    seg_h = None if segment_seconds in (None, "") else float(segment_seconds) / 3600
    if advisory(info, max_h):
        if warn is not None and (top < 1 or (seg_h is not None and seg_h > top + 1e-9)):
            hold = top if top >= 1 else round(float(fit_h) / float(info["k"]), 2)
            warn(f"at {where} ({rate_text(info)}) walltime {walltime} is expected to hold a {label}"
                 f"segment of at most {hold:g} h: not refused (no measured rate to refuse on); watch the "
                 f"first job's wall time.")
        return None
    if top < 1:
        if seg_h is not None and 0 < seg_h < 1:
            # a segment shorter than 1 h needs room for itself only
            if info["k"] * seg_h <= fit_h + 1e-9:
                return None
            return (f"{what} {_fmt(seg_h * 3600)} s does not fit walltime {walltime}: it needs "
                    f"{hms(info['k'] * seg_h)} ({where}, {rate_text(info)}, {config} rates"
                    + (f"; queue {queue!r} allows at most {walltime_max}" if max_h is not None else "")
                    + f"): raise pbs.walltime to at least {hms(info['k'] * seg_h)} or {remedy}.")
        need = hms(info["k"])
        return (f"pbs.walltime {walltime} cannot hold even a 1 h {label}segment at {where} "
                f"({rate_text(info)}, {config} rates): a 1 h segment needs {need}"
                + (f"; queue {queue!r} allows up to {walltime_max}" if max_h is not None else "")
                + f": raise pbs.walltime to at least {need}"
                + (" or use another queue" if max_h is not None and info["k"] > max_h else ""))
    if seg_h is None or seg_h <= 0 or seg_h <= top + 1e-9:
        return None
    return (f"{what} {_fmt(seg_h * 3600)} s does not fit walltime {walltime}: it needs "
            f"{hms(info['k'] * seg_h)} ({where}, {rate_text(info)}, {config} rates"
            + (f"; queue {queue!r} allows at most {walltime_max}" if max_h is not None else "")
            + f"). The longest segment that fits is {top * 3600} s ({top} h): {remedy}.")


def ncpus_error(ncpus_max, nodes, cores_per_node, queue, what, hpc="", config="config/machines.yaml"):
    """None, or why nodes x cores exceeds the queue's ncpus_max."""
    if ncpus_max is None or not cores_per_node or int(nodes) * int(cores_per_node) <= int(ncpus_max):
        return None
    n = int(nodes) * int(cores_per_node)
    return (f"{what}: {int(nodes)} nodes x {int(cores_per_node)} cores = {n} cores is more than queue "
            f"{queue!r} on {hpc} allows ({ncpus_max}, {config} queues): use another queue")


def segment_walltime_hours(stephours, rate=WALL_H_PER_SIM_H, margin=WALLTIME_MARGIN):
    """Expected walltime hours of a segment: margin x rate x stephours."""
    return float(margin) * float(rate) * float(stephours)


def segment_fits_walltime(stephours, wall_h, rate=WALL_H_PER_SIM_H, margin=WALLTIME_MARGIN):
    return segment_walltime_hours(stephours, rate, margin) <= float(wall_h) + 1e-9


def valid_stephours(window_h, delay_h, wall_h, rate=WALL_H_PER_SIM_H, margin=WALLTIME_MARGIN):
    """GRX segment hours that divide the coupled window, cover the delay and fit the walltime."""
    window_h = int(window_h)
    return [s for s in range(1, window_h + 1)
            if window_h % s == 0 and delay_h <= s and segment_fits_walltime(s, wall_h, rate, margin)]


def grt_valid_hours(window_h, wall_h, k, lattice_s=()):
    """GRT/standalone segment hours that fit the walltime and window and are multiples of lattice_s."""
    window_top = int(math.floor(float(window_h) + 1e-9))
    top = window_top if math.isinf(float(wall_h)) else min(max_segment_hours(wall_h, k), window_top)
    return [s for s in range(1, top + 1)
            if all(float(lat) <= 0 or (s * 3600) % float(lat) == 0 for lat in lattice_s if lat)]


def nearest_divisible_stop(coupling_start, window_h, delay_h, max_h):
    """(before, after): the nearest stop dates whose window some valid segment length divides."""
    start = _parse_dt(coupling_start)
    lo, hi = max(1, int(math.ceil(delay_h))), int(max_h)
    if hi < lo:
        return None, None
    cands = set()
    for s in range(lo, hi + 1):
        cands.add(int(window_h) // s * s)
        cands.add(-(-int(window_h) // s) * s)
    before = max((c for c in cands if 0 < c < window_h), default=None)
    after = min((c for c in cands if c > window_h), default=None)
    fmt = lambda h: None if h is None else start + timedelta(hours=h)
    return fmt(before), fmt(after)


def coupled_start(t0, delaycouplinghrs):
    """(RUN_STARTDATE, START_TOD) of the coupled case: t0 + delaycouplinghrs."""
    tc = _parse_dt(t0) + timedelta(hours=delaycouplinghrs)
    return tc.strftime("%Y-%m-%d"), tc.hour * 3600 + tc.minute * 60 + tc.second


def _check_ic(ic):
    if ic != 1:
        raise ValueError("gr_warm_up_time must be <= segment_duration")


def nres(ii, ic, delaycouplinghrs, stephours, dtRes=1800.0):
    """GAMERA restart index segment ii resumes from (-1 for spinup), in units of dtRes seconds.

    A segment start that is not a multiple of dtRes has no restart file: ValueError."""
    _check_ic(ic)
    if ii == 0:
        return -1
    if ii == 1:
        return 0
    start = (delaycouplinghrs + (ii - 2) * stephours) * 3600
    dtRes = float(dtRes)
    if dtRes <= 0 or start % dtRes:
        raise ValueError(f"GRX segment {ii - ic:02d} starts {start:.0f} s into the run, not a multiple of "
                         f"voltron.restart.dtRes ({dtRes:g} s), so it has no restart file: use a dtRes "
                         f"that divides the segment and warm-up lengths.")
    return int(start / dtRes)


def tfin(ii, ic, delaycouplinghrs, stephours, raiju=True, dtCouple=DTCOUPLE_DEFAULT):
    """voltron tFin (s from t0) of segment ii.

    Coupled: segment end + k*dtCouple + min(2, dtCouple/2), k = 2 on the first coupled segment and
    1 after, because CAM makes N+2 exchanges in a startup segment and N+1 in a continue segment."""
    _check_ic(ic)
    if ii <= ic:
        return float(min(ii, 1) * 3600 * delaycouplinghrs + 4)
    dtc = float(dtCouple)
    t1 = (delaycouplinghrs + (ii - 1 - ic) * stephours + stephours) * 3600
    k = 2 if ii == ic + 1 else 1
    # the pad keeps tFin inside the next coupling interval, clear of voltron's +0.25 s near-multiple bump (voltapp.F90)
    return float(t1 + k * dtc + min(2.0, dtc / 2.0))


def rpointer_suffix(ii, ic, coupling_start, stephours):
    """YYYY-MM-DD-SSSSS stamp of the CESM restart coupled segment ii resumes from."""
    t = _parse_dt(coupling_start) + timedelta(hours=(ii - 1 - ic) * stephours)
    return f"{t:%Y-%m-%d}-{t.hour * 3600 + t.minute * 60 + t.second:05d}"


def segment_id(job_name, ii, ic):
    """<job>-SPINUP, <job>-WARMUP-NN, then the coupled <job>-NN."""
    if ii == 0:
        return f"{job_name}-SPINUP"
    if ii <= ic:
        return f"{job_name}-WARMUP-{ii:02d}"
    return f"{job_name}-{ii - ic:02d}"


def generate(options, run_dir=None):
    """Retired."""
    raise NotImplementedError(
        "segment_solver.generate() is retired: the .xml decks come from makeitso.gen_grx_voltron "
        "and the .pbs jobs from engage.create_grx_pbs_scripts (engage.py --component grx)")
