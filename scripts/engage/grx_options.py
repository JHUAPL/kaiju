"""Defaults and checks of engage's coupled-run options (GRX and GRT), shared by engage's prompts and wizard.

Does not import engage (tui_form_engage runs inside engage's __main__).
"""
import datetime
import importlib
import json
import math
import os
import sys

import yaml

_HERE = os.path.dirname(os.path.abspath(__file__))
WAXRUN_DIR = os.path.normpath(os.path.join(_HERE, os.pardir, "waxrun"))
ENGAGE_MACHINES = os.path.join(_HERE, "config", "machines.yaml")

# Seconds: 4-h coupled segments and a 10-day spinup target.
GRX_SEGMENT_DURATION = "14400"
GRX_GCM_SPIN_UP_TIME = "864000"
RUN_DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S"


def segment_solver():
    """waxrun's segment_solver, from the same checkout."""
    if WAXRUN_DIR not in sys.path:
        sys.path.append(WAXRUN_DIR)
    return importlib.import_module("segment_solver")


def waxrun_descriptions():
    with open(os.path.join(WAXRUN_DIR, "config", "options_description.json"), encoding="utf-8") as f:
        return json.load(f)


def grx_machines():
    """waxrun's hpc_system valids that also have a stacks.grx entry in engage's machines.yaml."""
    wax = waxrun_descriptions()["simulation"]["hpc_system"].get("valids") or []
    with open(ENGAGE_MACHINES, encoding="utf-8") as f:
        machines = yaml.safe_load(f)["machines"]
    return [m for m in wax if "grx" in ((machines.get(m) or {}).get("stacks") or {})]


def grx_defaults(hpc_system=None):
    """{option path: default} of the engage options a GRX run defaults differently from GRT."""
    desc = waxrun_descriptions()
    machines = grx_machines()
    hpc = hpc_system if hpc_system in machines else (machines[0] if machines else "derecho")
    job = desc["job"].get(hpc) or desc["job"]["derecho"]
    return {
        "simulation.hpc_system": hpc,
        "simulation.start_date": desc["simulation"]["coupling_start"]["default"],
        "simulation.stop_date": desc["simulation"]["stop_time"]["default"],
        "simulation.use_segments": "Y",
        "simulation.segment_duration": GRX_SEGMENT_DURATION,
        "coupling.gcm_spin_up_time": GRX_GCM_SPIN_UP_TIME,
        "pbs.walltime": job["walltime"]["default"],
    }


def machine_error(hpc_system):
    machines = grx_machines()
    if hpc_system in machines:
        return None
    return (f"GRX runs only on {'|'.join(machines) or 'no machine'} (a WACCM-X stack in engage "
            f"config/machines.yaml); got hpc_system={hpc_system!r}")


def _dt(value):
    if isinstance(value, datetime.datetime):
        return value
    return datetime.datetime.strptime(str(value).strip(), RUN_DATETIME_FORMAT)


def _seconds(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def hours_error(name, value):
    s = _seconds(value)
    if s is None or s < 3600 or s % 3600:
        return f"GRX runs need {name} in whole hours (>= 3600 s); got {name}={value} s"
    return None


def _s(value):
    v = _seconds(value)
    if v is None:
        return str(value)
    return f"{v:.0f}" if float(v).is_integer() else f"{v:g}"


def seconds_list(hours):
    return "[" + ", ".join(f"{int(h) * 3600} ({int(h)} h)" for h in hours) + "]"


def window(start_date, stop_date):
    start, stop = _dt(start_date), _dt(stop_date)
    return start, stop, (stop - start).total_seconds()


def valid_segments(start_date, stop_date, gr_warm_up_time, walltime, hpc_system="derecho", queue=None):
    """Valid GRX segment lengths (hours) at engage's machines.yaml GRX rate."""
    _, _, window_s = window(start_date, stop_date)
    delay_h = (_seconds(gr_warm_up_time) or 0) / 3600
    info = rate_info(hpc_system, "grx", "Q", quiet=True)
    _, _, fit_h = fit_walltime(walltime, hpc_system, queue)
    return segment_solver().valid_stephours(int(window_s / 3600), delay_h, fit_h, rate=info["rate"],
                                            margin=info["margin"])


def segment_error(start_date, stop_date, segment_duration, gr_warm_up_time, walltime,
                  hpc_system="derecho", queue=None):
    """None, or why segment_duration is not whole hours dividing the GRX coupled window."""
    err = hours_error("segment_duration", segment_duration)
    if err:
        return err
    try:
        start, stop, window_s = window(start_date, stop_date)
    except (TypeError, ValueError):
        return None                                  # the dates' own validators report it
    step_h = int(_seconds(segment_duration)) // 3600
    if window_s <= 0 or window_s % (step_h * 3600):
        if window_s % 3600:
            return (f"the coupled window {start:%Y-%m-%dT%H:%M:%S} -> {stop:%Y-%m-%dT%H:%M:%S} is not a "
                    f"whole number of hours. Enter a stop_date with start_date's minutes and seconds "
                    f"({start:%M:%S}).")
        valid = valid_segments(start_date, stop_date, gr_warm_up_time, walltime, hpc_system, queue)
        return (f"segment_duration {_s(segment_duration)} s does not divide the coupled window "
                f"{start:%Y-%m-%dT%H:%M:%S} -> {stop:%Y-%m-%dT%H:%M:%S} ({window_s:.0f} s). "
                f"Valid values (seconds): {seconds_list(valid)}")
    return None


def warmup_error(gr_warm_up_time, segment_duration, start_date=None, stop_date=None, walltime=None,
                 hpc_system="derecho", queue=None):
    """None, or why gr_warm_up_time is not whole hours no longer than one segment."""
    err = hours_error("gr_warm_up_time", gr_warm_up_time)
    if err:
        return err
    seg = _seconds(segment_duration)
    if not seg or seg % 3600:
        return None                                  # segment_duration's own check reports it
    delay_h, step_h = int(_seconds(gr_warm_up_time)) // 3600, int(seg) // 3600
    if math.ceil(delay_h / step_h) != 1:
        valid = f" Enter a gr_warm_up_time of at most {_s(seg)} s."
        if start_date and stop_date and walltime:
            try:
                fits = valid_segments(start_date, stop_date, gr_warm_up_time, walltime, hpc_system, queue)
            except (TypeError, ValueError):
                fits = []
            if fits:
                valid += f" Segment values (seconds) that fit this warm-up: {seconds_list(fits)}"
        return (f"gr_warm_up_time {_s(gr_warm_up_time)} s is longer than segment_duration {_s(seg)} s."
                + valid)
    return None


def coupling_start_error(start_date, dtCouple=None, spinup_dtime=None):
    """None, or why start_date's second of day is not a multiple of the spin-up dtime and dtCouple."""
    ss = segment_solver()
    return ss.coupling_start_error(start_date, dtCouple,
                                   spinup_dtime if spinup_dtime else ss.SPINUP_DTIME_DEFAULT,
                                   name="simulation.start_date (the GRX coupling start)")


def cadence_error(dtCouple, tSpin=None, gr_warm_up_time=None, segment_duration=None, start_date=None):
    """None, or the first coupled-cadence rule a GRX run breaks (no dtOut rules: WACCM-X has its own)."""
    return segment_solver().cadence_error(dtCouple, tSpin=tSpin, gr_warm_up_time=gr_warm_up_time,
                                          segment_duration=segment_duration, coupled_start=start_date)


_DOC = {}


def _engage_doc():
    """engage config/machines.yaml, re-read only when the file changes."""
    key = os.stat(ENGAGE_MACHINES).st_mtime_ns
    if key not in _DOC:
        with open(ENGAGE_MACHINES, encoding="utf-8") as f:
            _DOC.clear()
            _DOC[key] = yaml.safe_load(f) or {}
    return _DOC[key]


def engage_machine(hpc_system):
    return (_engage_doc().get("machines") or {}).get(str(hpc_system)) or {}


def local_prefix_lines(prefix):
    """Bash lines putting a prefix of user-built libraries on the library, include and bin paths."""
    L = [f'export PREFIX="{prefix}"']
    for var, sub in (("LIBRARY_PATH", "lib"), ("LD_LIBRARY_PATH", "lib"), ("CPATH", "include")):
        L.append(f'export {var}="$PREFIX/{sub}${{{var}:+:${var}}}"')
    return L + ['export PATH="${PATH}:$PREFIX/bin"']


def load_engage_descriptions(path=None):
    """option_engage_descriptions.json with pbs.<machine>.local_modules defaults from machines.yaml."""
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "option_engage_descriptions.json")
    with open(path, encoding="utf-8") as f:
        od = json.load(f)
    for machine, block in (od.get("pbs") or {}).items():
        field = (block or {}).get("local_modules")
        if isinstance(field, dict):
            prefix = engage_machine(machine).get("local_prefix")
            if prefix:
                field["default"] = local_prefix_lines(prefix)
    return od


def queue_limits(hpc_system, queue):
    """{ncpus_max, walltime_max} of a queue in engage's machines.yaml; {} when not recorded."""
    queues = engage_machine(hpc_system).get("queues")
    if isinstance(queues, dict):
        return dict(queues.get(queue) or {})
    return {}


# Resolutions the rates tables are keyed by: coupled TIE-GCM per GAMERA grid, and WACCM-X.
GRT_COUPLED_RES = {"D": 2.5, "Q": 1.25, "O": 0.625, "H": 0.625}
GRX_RES = "ne30"
_NOTED = set()


def component_res(component, grid):
    if component == "grt":
        return GRT_COUPLED_RES.get(str(grid))
    if component == "grx":
        return GRX_RES
    return None


# Off by default: a print from inside the full-screen wizard would corrupt its display.
ANNOUNCE = False


def announce(on=True):
    global ANNOUNCE
    ANNOUNCE = bool(on)


def rate_info(hpc_system, component, grid=None, res=None, quiet=False):
    """Rate of a coupled `component` segment on `hpc_system` from engage's machines.yaml rates."""
    if res is None:
        res = component_res(component, grid)
    info = segment_solver().lookup_rate(engage_machine(hpc_system).get("rates") or {}, component,
                                        grid, res)
    if info.get("note") and not quiet and ANNOUNCE and (hpc_system, component, grid, res) not in _NOTED:
        _NOTED.add((hpc_system, component, grid, res))
        print(info["note"] + f" ({hpc_system}: engage config/machines.yaml rates)")
    return info


def fit_walltime(walltime, hpc_system, queue):
    wall_max = queue_limits(hpc_system, queue).get("walltime_max") if queue else None
    return segment_solver().fit_hours(walltime, wall_max)


def advisory(info, max_h):
    return segment_solver().advisory(info, max_h)


_ADVISED = set()


def _advise_once(key, text):
    if ANNOUNCE and key not in _ADVISED:
        _ADVISED.add(key)
        print(f"WARNING: {text}")


def walltime_error(segment_duration, walltime, hpc_system="derecho", queue=None, component="grx",
                   grid=None, res=None):
    """None, or why a job of `walltime` on `queue` cannot run this run's coupled segments."""
    info = rate_info(hpc_system, component, grid, res)
    label = {"grx": "GRX", "grt": "GRT"}.get(component, component)
    where = f"{label} grid {grid}/{info['key'].split('/', 1)[-1]} on {hpc_system}" if grid else f"{label} on {hpc_system}"
    seg = _seconds(segment_duration)
    wall_max = queue_limits(hpc_system, queue).get("walltime_max") if queue else None
    return segment_solver().walltime_error(
        info, walltime, wall_max, seg if seg and seg > 0 else None, queue=queue, hpc=hpc_system,
        config="engage config/machines.yaml", where=where, label=f"coupled {label} ",
        warn=lambda text: _advise_once((hpc_system, component, grid, str(walltime)), text))


def grt_window_error(start_date, stop_date):
    """None, or why a GRT coupled window is not a whole number of hours."""
    try:
        start, stop, window_s = window(start_date, stop_date)
    except (TypeError, ValueError):
        return None
    if window_s > 0 and window_s % 3600:
        return (f"the coupled window {start:%Y-%m-%dT%H:%M:%S} -> {stop:%Y-%m-%dT%H:%M:%S} "
                f"({window_s / 3600:g} h) is not a whole number of hours: coupled segments are whole "
                f"hours, so stop_date must share start_date's minute and second ({start:%M:%S})")
    return None


# Default voltron dtRes, dtCouple and dtOut (s): GRT segment lengths must be multiples of each.
GRT_LATTICE_S = (1800, 5, 60)


def segment_choices(component, hpc_system, queue, walltime, grid, start_date, stop_date,
                    gr_warm_up_time=None, lattice_s=GRT_LATTICE_S):
    """(valid whole-hour segment lengths, rate_info, max_seg_h) of a coupled run."""
    ss = segment_solver()
    try:
        _, max_h, fit_h = fit_walltime(walltime, hpc_system, queue)
        _, _, window_s = window(start_date, stop_date)
    except (TypeError, ValueError):
        return [], rate_info(hpc_system, component, grid, quiet=True), 0
    info = rate_info(hpc_system, component, grid)
    top = ss.max_segment_hours(fit_h, info["k"])
    # an advisory (unknown) rate caps nothing; the default stays conservative
    adv = advisory(info, max_h)
    cap_h, k = (math.inf, 1.0) if adv else (fit_h, info["k"])
    if component == "grx":
        delay_h = (_seconds(gr_warm_up_time) or 0) / 3600
        valid = ss.valid_stephours(int(window_s / 3600), delay_h, cap_h,
                                   rate=1.0 if adv else info["rate"],
                                   margin=1.0 if adv else info["margin"])
    else:
        valid = ss.grt_valid_hours(window_s / 3600, cap_h, k, lattice_s)
    if adv and valid:
        fit = [v for v in valid if v <= max(top, 1)]
        info = dict(info, suggested=max(fit) if fit else min(valid))
    elif valid:
        info = dict(info, suggested=max(valid))
    return valid, info, top


def suggested_seconds(valid, info, top=1):
    """The segment_duration default in seconds."""
    best = info.get("suggested", max(valid) if valid else max(top, 1))
    return str(int(best) * 3600)


def segment_suggestion(component, hpc_system, queue, walltime, grid, valid, info, top):
    best = int(suggested_seconds(valid, info, top)) // 3600
    return (f"Suggested segment for walltime {walltime} on {queue or hpc_system}: {best * 3600} s ({best} h); "
            f"valid values (seconds): {seconds_list(valid)}")


def segment_hours_error(component, hpc_system, queue, walltime, grid, start_date, stop_date,
                        segment_duration, gr_warm_up_time=None):
    """None, or why `segment_duration` (seconds) is not one of segment_choices' valid lengths."""
    label = {"grx": "GRX", "grt": "GRT"}.get(component, component)
    valid, info, top = segment_choices(component, hpc_system, queue, walltime, grid, start_date,
                                       stop_date, gr_warm_up_time)
    if hours_error("segment_duration", segment_duration):
        return (f"segment_duration {_s(segment_duration)} s is not whole hours (a multiple of 3600 s). "
                f"Valid values (seconds): {seconds_list(valid)}")
    step_h = int(_seconds(segment_duration)) // 3600
    if component != "grx":
        # GRT runs a segment longer than the window as one window-length segment
        try:
            window_h = int(window(start_date, stop_date)[2] // 3600)
        except (TypeError, ValueError):
            window_h = 0
        if 1 <= window_h < step_h:
            step_h = window_h
    if step_h in valid:
        if info.get("suggested") is not None and step_h > max(top, 1) and component != "grx":
            walltime_error(segment_duration, walltime, hpc_system, queue, component, grid)   # advisory WARNING
        return None
    if not valid and component == "grx":
        try:
            start, stop, window_s = window(start_date, stop_date)
            delay_h = (_seconds(gr_warm_up_time) or 0) / 3600
            before, after = segment_solver().nearest_divisible_stop(start, window_s / 3600, delay_h, top)
        except (TypeError, ValueError):
            before = after = None
        near = " or ".join(f"{t:%Y-%m-%dT%H:%M:%S}" for t in (before, after) if t is not None)
        return (f"no segment of at most {top * 3600} s fits walltime {walltime}, divides the coupled "
                f"window ({window_s:.0f} s) and covers the {delay_h * 3600:.0f} s warm-up. "
                + (f"Enter stop_date {near}." if near else "Change stop_date or pbs.walltime."))
    if component == "grx":
        why = (f"it must fit the walltime (at most {top * 3600} s), divide the coupled window and be at "
               f"least the warm-up")
    else:
        why = f"it must be whole hours that fit the walltime (at most {top * 3600} s)"
    return (f"segment_duration {_s(segment_duration)} s is not a valid coupled {label} segment for walltime "
            f"{walltime} on {queue or hpc_system}: {why}. Valid values (seconds): {seconds_list(valid)}")


def ncpus_error(hpc_system, queue, nodes, what):
    """None, or why nodes x cores exceed the queue's ncpus_max."""
    cap = (queue_limits(hpc_system, queue) if queue else {}).get("ncpus_max")
    cores = int(engage_machine(hpc_system).get("cores_per_node") or 0)
    return segment_solver().ncpus_error(cap, nodes, cores, queue, what, hpc_system,
                                        "engage config/machines.yaml")
