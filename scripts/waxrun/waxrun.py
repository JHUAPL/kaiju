#!/usr/bin/env python3
"""waxrun.py: WACCM-X (GRX) run configurator, the WACCM-X analog of tiegcmrun.

Builds the <job>_spinup and <job> (coupled) CAM cases with create_newcase and the jobs around them.
"""

import argparse
import copy
import glob
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path

import case_solver
import forcing_solver
import ic_solver
import nudge_solver
import segment_solver
import templates
import tree_solver

HERE = Path(__file__).resolve().parent
CONFIG_DIR = HERE / "config"
OPTIONS_DESCRIPTION = CONFIG_DIR / "options_description.json"

API_VERSION = "1.0.0"
LEVELS = ["BENCH", "BASIC", "INTERMEDIATE", "EXPERT"]
PHASES = ("spinup", "coupled", "both")
PARAMETERS_JSON = "waxrun_parameters.json"

RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
RESET = "\033[0m"
# a cam.* value that already is a quoted Fortran string (or a list of them)
_QUOTED_RHS = re.compile(r"""^\s*(['"]).*\1\s*$""", re.S)
_BOOL_TRUE = {"yes", "y", "true", "t", "1", ".true."}
_BOOL_FALSE = {"no", "n", "false", "f", "0", ".false."}

# set by waxrun; cam.* overrides may not touch them
DERIVED_NL_KEYS = {"ncdata", "inithist", "nudge_model", "nudge_path", "nudge_file_template",
                   "nudge_times_per_day", "model_times_per_day", "nudge_beg_year", "nudge_beg_month",
                   "nudge_beg_day", "nudge_end_year", "nudge_end_month", "nudge_end_day",
                   "nudge_vwin_lindex", "nudge_vwin_hindex"}
PRODUCT_COMMENT = "PRODUCT OF SPINUP - exists only after spinup runs"
# option keys that are not options_description.json leaves
INTERNAL_KEYS = {"api_version", "mode", "ic.pick", "ic.stripped"}
_WALLTIME_RE = re.compile(r"^\d+:\d\d:\d\d$")
TASKS_PER_NODE = 128
# dtime (s) at which spinup.wall_s_per_simday was measured; the estimate scales with the ratio
REFERENCE_DTIME = 300

LAST = {}     # the last bridge run: {"options", "spinup", "coupled", "path"}
_TREE_CACHE = {}
_OPTIONS_FILE = None
CLONE_STEPS = "scripts/waxrun/README.md"

check_tree = tree_solver.check_tree


def _warn(s):
    print(s, file=sys.stderr)


def _bool_of(value):
    """True/False for a bool or a yes/no spelling; None when `value` is neither."""
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        low = value.strip().lower()
        if low in _BOOL_TRUE:
            return True
        if low in _BOOL_FALSE:
            return False
    return None


def _as_bool(value, leaf):
    """A boolean leaf as a bool (bool("N") would be True); anything else is refused."""
    b = _bool_of(value)
    if b is None:
        raise SystemExit(f"[waxrun] {leaf}={value!r} is not a yes/no value (true/false, Y/N, yes/no)")
    return b


def _dt(v):
    return segment_solver._parse_dt(v)


def _iso(v, field="date"):
    """A run date/time as YYYY-MM-DDThh:mm:ss; anything else raises ValueError."""
    return segment_solver.parse_run_datetime(v, field).strftime(segment_solver.RUN_DATETIME_FORMAT)


def _run_dt(v, field):
    """Parse a run date from the options; SystemExit on a bad value."""
    try:
        return segment_solver.parse_run_datetime(v, field)
    except ValueError as e:
        raise SystemExit(f"[waxrun] {e}") from None


def _hms(secs):
    secs = int(secs)
    return f"{secs // 3600:02d}:{secs % 3600 // 60:02d}:{secs % 60:02d}"


def _secs(hms):
    parts = [int(p) for p in str(hms).split(":")]
    while len(parts) < 3:
        parts.append(0)
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


def _hours(value, name):
    s = int(round(float(value)))
    if s % 3600:
        raise SystemExit(f"{name}={value!r} is not a whole number of hours (GRX runs need whole-hour "
                         "gr_warm_up_time and segment_duration)")
    return s // 3600


def _whole_hours(value, leaf):
    if float(value) != int(float(value)):
        raise SystemExit(f"{leaf}={value!r} is not a whole number of hours (GRX runs need whole-hour "
                         "gr_warm_up_time and segment_duration)")
    return int(float(value))


def _tod(t):
    return t.hour * 3600 + t.minute * 60 + t.second


def _deep_update(dst, src, keep_none=False):
    """Merge src into dst; a None in src is skipped unless keep_none."""
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v, keep_none)
        elif v is not None or (keep_none and not isinstance(dst.get(k), dict)):
            dst[k] = v


def _leaf_paths(d, prefix=""):
    out = []
    for k, v in d.items():
        p = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict) and v and not k.startswith("_"):
            out += _leaf_paths(v, p)
        else:
            out.append(p)
    return out


def layout(job_name, run_dir):
    """The run-directory layout for (job_name, run_dir); no prompting or I/O."""
    run_dir = os.path.abspath(str(run_dir))
    return {
        "job_name": job_name,
        "run_dir": run_dir,
        "rundir": run_dir,
        "caseroot_spinup": os.path.join(run_dir, "wax_case", f"{job_name}_spinup"),
        "caseroot_coupled": os.path.join(run_dir, "wax_case", job_name),
        "cime_output_root": os.path.join(run_dir, "cesm"),
        "shared_dir": os.path.join(run_dir, "shared"),
        "nudging_gen_dir": os.path.join(run_dir, "nudging_gen"),
        "assembled_dir": os.path.join(run_dir, f"met_data_{job_name}"),
        "dryrun_dir": os.path.join(run_dir, "dryrun"),
    }


def _expand(value, lay):
    if isinstance(value, str):
        return value.replace("<run_directory>", lay["run_dir"]).replace("<job>", lay["job_name"])
    if isinstance(value, list):
        return [_expand(v, lay) for v in value]
    return value


ENGAGE_MACHINES = HERE.parent / "engage" / "config" / "machines.yaml"


def _engage_machines():
    """engage config/machines.yaml `machines` ({} when unreadable)."""
    try:
        import yaml
        return yaml.safe_load(ENGAGE_MACHINES.read_text())["machines"] or {}
    except (OSError, ImportError, KeyError, TypeError):
        return {}


def grx_rate(hpc, grid="Q", res="ne30"):
    """The coupled GRX segment rate from engage config/machines.yaml rates.grx."""
    rates = (_engage_machines().get(hpc) or {}).get("rates") or {}
    return segment_solver.lookup_rate(rates, "grx", grid, res)


def queue_max_walltime(hpc):
    """{queue: walltime_max} of a machine's queues."""
    queues = (_engage_machines().get(str(hpc)) or {}).get("queues") or {}
    if not isinstance(queues, dict):
        return {}
    return {q: lim["walltime_max"] for q, lim in queues.items() if (lim or {}).get("walltime_max")}


def project_required(hpc):
    """True when the machine charges every job to an account."""
    return bool((_engage_machines().get(str(hpc)) or {}).get("account_required"))


def _jobgen():
    """engage's jobgen module, loaded once by file path."""
    mod = sys.modules.get("waxrun_jobgen")
    if mod is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location("waxrun_jobgen", HERE.parent / "engage" / "jobgen.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["waxrun_jobgen"] = mod
        spec.loader.exec_module(mod)
    return mod


def project_choices(hpc):
    """The user's project codes to offer at the simulation.project prompt ([] if none apply)."""
    machine = _engage_machines().get(str(hpc)) or {}
    if not machine.get("account_groups") or not project_required(hpc) or machine.get("account_default"):
        return []
    return _jobgen().user_projects()


def project_prompt_default(hpc):
    choices = project_choices(hpc)
    return choices[0] if len(choices) == 1 else None


def project_error(value, hpc):
    """None, or why simulation.project cannot be blank on the machine."""
    if value not in (None, "") and str(value).strip() not in ("", "None", "null"):
        return None
    if not project_required(hpc):
        return None
    return (f"simulation.project is required on {hpc}: the cases, the nudging job and every GRX "
            f"job are charged to it, and there is no shared default (the account is per user)")


def _check_project(opt):
    sim = opt["simulation"]
    err = project_error(sim.get("project"), sim.get("hpc_system", "derecho"))
    if err:
        raise SystemExit(f"[waxrun] {err}. Add simulation.project to the options JSON (engage: "
                         f"pbs.account_name).")


CONFIG_BATCH_PATHS = ("ccs_config/machines/{m}/config_batch.xml", "ccs_config/machines/config_batch.xml")


def cime_queues(srcroot, machine):
    """(queues, path) for the machine in the tree's CIME config_batch.xml, or (None, None)."""
    if not srcroot:
        return None, None
    for rel in CONFIG_BATCH_PATHS:
        path = os.path.join(str(srcroot), rel.format(m=machine))
        if not os.path.isfile(path):
            continue
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            continue
        queues = [q.text.strip() for bs in root.iter("batch_system") if bs.get("MACH") == machine
                  for q in bs.iter("queue") if q.text and q.text.strip()]
        if queues:
            return queues, path
    return None, None


def queue_error(queue, srcroot, machine, hpc, nodes=None, what="the spinup case"):
    """None, or why `queue` cannot run the job (unknown to CIME, or too many cores for the queue)."""
    queues, path = cime_queues(srcroot, machine)
    if queues is not None and queue not in queues:
        return (f"job.{hpc}.queue={queue!r} (engage pbs.queue) is not a queue CIME defines for {machine} "
                f"({path}: {', '.join(queues)}): pick one of {', '.join(queues)}")
    limits = ((_engage_machines().get(str(hpc)) or {}).get("queues") or {})
    cap = limits.get(queue, {}).get("ncpus_max") if isinstance(limits, dict) else None
    if nodes is None:
        return None
    err = segment_solver.ncpus_error(cap, nodes, TASKS_PER_NODE, queue, what, hpc,
                                     "engage config/machines.yaml")
    return err and f"{err} (job.{hpc}.queue=main; develop is for the 1-node nudging generator)"


def load_descriptions(path=OPTIONS_DESCRIPTION):
    """options_description.json with the account and walltime defaults from engage's machines.yaml."""
    desc = json.loads(Path(path).read_text())
    try:
        import yaml
        machines = yaml.safe_load(ENGAGE_MACHINES.read_text())["machines"]
    except (OSError, ImportError, KeyError, TypeError):
        return desc
    sim = desc.get("simulation", {})
    hpc = (sim.get("hpc_system") or {}).get("default")
    account = _jobgen().account_value((machines.get(hpc) or {}).get("account_default"))
    if account and isinstance(sim.get("project"), dict):
        sim["project"]["default"] = account
    for name, job in (desc.get("job") or {}).items():
        walltime = (machines.get(name) or {}).get("walltime_default") if isinstance(job, dict) else None
        if walltime and isinstance(job.get("walltime"), dict):
            job["walltime"]["default"] = walltime
    return desc


def _defaults(desc):
    out = {}
    for k, v in desc.items():
        if k.startswith("_") or k == "version":
            continue
        if isinstance(v, dict) and "LEVEL" in v:
            out[k] = copy.deepcopy(v.get("default"))
        elif isinstance(v, dict):
            out[k] = _defaults(v)
    return out


def _at_or_below(level, mode):
    try:
        return LEVELS.index(level) <= LEVELS.index(mode)
    except ValueError:
        return True


# leaves that must be YYYY-MM-DDThh:mm:ss
_DATE_LEAVES = ("coupling_start", "stop_time", "spinup_start")


_NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
PATH_TYPES = ("file", "dir", "path")
_PATH_NONE = ("", "none", "null")


def _anchor(value, run_dir):
    """Expand ~ and join a relative path onto run_dir; placeholders and $VAR are left as they are."""
    s = os.path.expanduser(str(value).strip())
    if run_dir is None or os.path.isabs(s) or s.startswith("$") or "<" in s:
        return s
    return os.path.normpath(os.path.join(str(run_dir), s))


def _path_value(path, kind, spec, value, run_dir):
    """An anchored file/dir/path value; a file or dir is checked to exist when it can be."""
    if not isinstance(value, str):
        raise ValueError(f"{path}={value!r} is not a path")
    if value.strip().lower() in _PATH_NONE:
        return None
    tool = spec.get("anchor") == "tool"
    full = _anchor(value, HERE if tool else run_dir)
    checkable = os.path.isabs(full) and "<" not in full and "$" not in full
    if kind == "file" and checkable and not os.path.isfile(full):
        raise ValueError(f"{path}={value!r} is not a file" if os.path.exists(full)
                         else f"{path}={value!r} does not exist")
    if kind == "dir" and checkable and not os.path.isdir(full):
        raise ValueError(f"{path}={value!r} is not a directory")
    return value.strip() if tool else full


def _number(path, value, integral):
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{path}={value!r} is not a number")
    if isinstance(value, str):
        try:
            value = float(value.strip())
        except ValueError:
            raise ValueError(f"{path}={value!r} is not a number") from None
    if not math.isfinite(value):
        raise ValueError(f"{path}={value!r} is not a number")
    if integral:
        if not float(value).is_integer():
            raise ValueError(f"{path}={value!r} is not a whole number")
        return int(value)
    return value


def _range_error(path, spec, value):
    lo, hi = spec.get("min"), spec.get("max")
    if lo is not None and (value <= lo if spec.get("exclusive_min") else value < lo):
        return f"{path}={value!r} must be {'>' if spec.get('exclusive_min') else '>='} {lo}"
    if hi is not None and value > hi:
        return f"{path}={value!r} must be <= {hi}"
    return None


def sourcemods_choices():
    """none and the set names in config/sourcemods/ (offered at the prompts, not enforced)."""
    return ["none"] + case_solver.sourcemods_sets()


def _sourcemods_value(path, value):
    """none, a set name or an absolute path, checked by case_solver.sourcemods; ValueError with one line."""
    if not isinstance(value, str):
        raise ValueError(f"{path}={value!r} is not a SourceMods set")
    s = value.strip()
    if s.lower() in ("none", "null"):
        return "none"
    case_solver.sourcemods(path.split(".")[0], s)
    return s


def _coupled_sourcemods(opt):
    """coupled.sourcemods; absent (an older saved file) or null is none."""
    return opt["coupled"].get("sourcemods") or "none"


def _leaf_value(path, spec, value, run_dir=None):
    """The validated value of one option leaf, typed by spec["type"] or else the default's type."""
    if value is None:
        return None
    key = path.rsplit(".", 1)[-1]
    dflt = spec.get("default")
    kind = spec.get("type")
    if kind in PATH_TYPES:
        return _path_value(path, kind, spec, value, run_dir)
    if kind == "sourcemods":
        return _sourcemods_value(path, value)
    if kind == "name":
        if not isinstance(value, str) or not _NAME_RE.match(value.strip()):
            raise ValueError(f"{path}={value!r} may hold only letters, digits, '.', '_' and '-'")
        value = value.strip()
    elif kind == "str_list" or (kind is None and isinstance(dflt, list)):
        if not isinstance(value, list):
            raise ValueError(f"{path}={value!r} must be a list")
        want = str if kind == "str_list" else (type(dflt[0]) if dflt else None)
        bad = [v for v in value if want is not None and (not isinstance(v, want) or isinstance(v, bool))]
        if bad:
            raise ValueError(f"{path}: {bad[0]!r} is not a {want.__name__}")
        if spec.get("items") in PATH_TYPES:
            value = [_path_value(f"{path}[{i}]", spec["items"], {}, v, run_dir) for i, v in enumerate(value)]
    elif isinstance(dflt, bool):
        b = _bool_of(value)
        if b is None:
            raise ValueError(f"{path}={value!r} is not a yes/no value (true/false, Y/N, yes/no)")
        value = b
    elif kind in ("int", "float") or isinstance(dflt, (int, float)):
        typed = value
        value = _number(path, value, kind == "int" or (kind is None and isinstance(dflt, int)))
        if kind == "float" and isinstance(typed, str) and float(value).is_integer():
            value = int(value)        # '5' -> 5, as engage hands the value over
        err = _range_error(path, spec, value)
        if err:
            raise ValueError(err)
    valids = spec.get("valids")
    if valids and str(value) not in [str(v) for v in valids]:
        raise ValueError(f"{path}={value!r} is not a valid value (valid {sorted(str(v) for v in valids)})")
    if key == "walltime" and not _WALLTIME_RE.match(str(value)):
        raise ValueError(f"{path}={value!r} is not a walltime; expected HH:MM:SS (e.g. 12:00:00)")
    if key == "start_file_date" and not ic_solver.START_FILE_DATE_RE.fullmatch(str(value).strip()):
        raise ValueError(f"{path}={value!r} must be {ic_solver.START_FILE_DATE_FORMS} "
                         f"(e.g. 2004-03-06T00:00:00)")
    if key in _DATE_LEAVES:
        segment_solver.parse_run_datetime(value, path)
    return value


def _prompt_leaf(key, spec, opt_sec, path=None, shown=None, check=None, run_dir=None):
    """Prompt for one option; `check(value)` returning a message asks again."""
    default = opt_sec.get(key, spec.get("default"))
    valids = spec.get("valids")
    description = spec.get("description")
    warning = spec.get("warning")
    is_bool = isinstance(spec.get("default"), bool)
    if warning:
        print(f"{YELLOW}{warning}{RESET}")
    prompt = spec["prompt"]
    if valids:
        prompt += f" ({'|'.join(str(v) for v in valids)})"
    if spec.get("type") == "sourcemods":
        prompt += f" ({'|'.join(sourcemods_choices())} or a path)"
    if path == "simulation.project" and default in (None, ""):
        hpc = opt_sec.get("hpc_system", "derecho")
        choices = project_choices(hpc)
        if choices:
            prompt += f" (your projects: {'|'.join(choices)})"
        default = project_prompt_default(hpc)
    if default is not None:
        shown = "Y" if default is True else "N" if default is False else default
    if shown is not None:
        prompt += f" [{GREEN}{shown}{RESET}]"
    path = path or key
    while True:
        try:
            ans = input(f"{prompt}: ").strip()
        except EOFError:
            raise SystemExit(f"\n[waxrun] input ended at the prompt for {path}; nothing was written") from None
        if ans == "?":
            print(f"{YELLOW}{description or 'No description available.'}{RESET}")
            continue
        if ans == "":
            if path == "simulation.project":
                err = project_error(default, opt_sec.get("hpc_system", "derecho"))
                if err:
                    print(f"{RED}{err}: enter your project code{RESET}")
                    continue
            value = default
        else:
            try:
                value = _answer_value(ans, spec, path, valids, is_bool, run_dir)
            except ValueError as e:
                print(f"{RED}{e}{RESET}")
                continue
        err = check(value) if check else None
        if err:
            print(f"{RED}{err}{RESET}")
            continue
        opt_sec[key] = value
        return


def _answer_value(ans, spec, path, valids, is_bool, run_dir=None):
    """A non-blank answer as the leaf value; ValueError with the line to show."""
    if is_bool:
        b = _bool_of(ans)
        if b is None:
            raise ValueError(f"Please enter yes/no (got {ans!r}).")
        return b
    if valids and ans not in [str(v) for v in valids]:
        raise ValueError(f"{ans!r} is not a valid value. Choose from: {'|'.join(str(v) for v in valids)}")
    if spec.get("type") != "sourcemods" and not isinstance(spec.get("default"), (int, float)) \
            and ans.lower() in ("none", "null"):
        return None
    return _leaf_value(path, spec, ans, run_dir)


# leaves asked only when an earlier answer in their section calls for them
ASK_IF = {"nudging.nudging_dir": lambda sec: _bool_of(sec.get("generate")) is False,
          "coupled.cesm_exe": lambda sec: sec.get("build") == "verify_only",
          "coupled.sourcemods": lambda sec: sec.get("build") != "verify_only"}


def asked(path, opt_sec):
    cond = ASK_IF.get(path)
    return cond is None or cond(opt_sec)


def _run_dir_of(opt):
    """The absolute run folder; relative paths the user types are anchored to it."""
    return os.path.abspath(os.path.expanduser(str((opt.get("simulation") or {}).get("run_directory") or ".")))


def _prompt_section(desc_sec, opt_sec, mode, prefix="", skip=(), root=None):
    for key, spec in desc_sec.items():
        if key.startswith("_"):
            continue
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(spec, dict) and "LEVEL" in spec:
            if path in skip or not spec.get("prompt") or not _at_or_below(spec["LEVEL"], mode):
                continue
            if not asked(path, opt_sec):
                continue
            _prompt_leaf(key, spec, opt_sec, path, run_dir=_run_dir_of(root) if root is not None else None)
        elif isinstance(spec, dict):
            opt_sec.setdefault(key, {})
            _prompt_section(spec, opt_sec[key], mode, path, skip, root)


IC_DATE_LEAVES = ("ic.start_file_date", "ic.spinup_start")


def _hm(tod):
    return f" {tod // 3600:02d}:{tod % 3600 // 60:02d}" if tod else ""


def _spinup_dtime(opt):
    return opt["spinup"].get("dtime") or segment_solver.SPINUP_DTIME_DEFAULT


def default_spinup(pick, coupling_start, opt):
    """(start, note): the default spinup start for restart `pick`, kept within min..max hours."""
    sp = opt["spinup"]
    start, moved = ic_solver.default_spinup_start(pick, coupling_start, float(sp["min_hours"]),
                                                  float(sp["max_hours"]), _spinup_dtime(opt))
    r = ic_solver.redated(pick, coupling_start)
    redate = f"re-dated to {r:%Y-%m-%d}{_hm(r.hour * 3600 + r.minute * 60)}"
    hours = (coupling_start - start).total_seconds() / 3600
    if not moved:
        return start, f"{redate} ({hours:g} h of spinup)"
    raw = (coupling_start - r).total_seconds() / 3600
    return start, (f"{redate} would be {raw:g} h of spinup; default moved to "
                   f"{start:%Y-%m-%dT%H:%M:%S} ({hours:g} h)")


def start_file_choices(opt):
    """The restart times in ic.start_file_dir, or [] when it cannot be listed."""
    try:
        return ic_solver.series_stamps(ic_solver.enumerate_candidates(opt["ic"]["start_file_dir"]))
    except (OSError, ValueError, KeyError, TypeError):
        return []


def auto_pick(opt, start_file_date=None):
    """(day, spinup_start, note): the defaults of the two IC date prompts."""
    w = _window(opt)
    target = w["coupling_start"] - timedelta(seconds=int(opt["spinup"]["target_seconds"]))
    cands = ic_solver.enumerate_candidates(opt["ic"]["start_file_dir"])
    pick, _, _ = _native_pick(cands, w["coupling_start"], target, float(opt["spinup"]["min_hours"]),
                              float(opt["spinup"]["max_hours"]), start_file_date,
                              opt["simulation"]["res"], warn=False)
    start, redate = default_spinup(pick, w["coupling_start"], opt)
    day = ic_solver.stamp_iso(pick)
    what = "nearest restart" if start_file_date is None else "restart"
    return day, start.strftime(segment_solver.RUN_DATETIME_FORMAT), f"{what}: {day}, {redate}"


def resolve_start_file_date(opt, value):
    """The one restart time a start file date answer names, else None."""
    if value is None:
        return None
    try:
        cands = ic_solver.enumerate_candidates(opt["ic"]["start_file_dir"])
        times = ic_solver.series_stamps([c for c in cands if ic_solver.matches_date(c, value)])
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return times[0] if len(times) == 1 else None


def auto_pick_failure(e):
    text = str(getattr(e, "code", None) or e).strip()
    return f"no restart pick yet: {text.splitlines()[0] if text else type(e).__name__}"


def ic_date_error(opt, path, value, spinup_default=None):
    """None, or why an IC date answer is refused."""
    if path == "ic.start_file_date":
        if value is None:
            return None
        try:
            auto_pick(opt, value)
        except SystemExit as e:
            text = str(e.code or "").strip()
            if text.startswith(("Start file date", "[waxrun] ic.start_file_date")):
                return text.splitlines()[0]
        except Exception:
            pass
        return None
    start = value if value is not None else spinup_default
    if start is None:
        return None
    try:
        cs_ = _run_dt(opt["simulation"]["coupling_start"], "simulation.coupling_start")
    except (SystemExit, ValueError):
        return None
    sp = opt["spinup"]
    return ic_solver.spinup_start_error(
        start, cs_, float(sp["min_hours"]), float(sp["max_hours"]),
        _spinup_dtime(opt) if value is not None else None)


def _prompt_ic_dates(opt, desc, mode, filled):
    """Prompt for the start file date, then the spinup start date."""
    ic = opt["ic"]
    for path in IC_DATE_LEAVES:
        key = path.split(".")[1]
        spec = desc["ic"][key]
        if path in filled or not _at_or_below(spec["LEVEL"], mode):
            continue
        sfd, derived = ic.get("start_file_date"), None
        if ic.get(key) is None:
            try:
                day, start, note = auto_pick(opt, sfd if key == "spinup_start" else None)
                derived = start if key == "spinup_start" else day
            except (SystemExit, Exception) as e:
                note = auto_pick_failure(e)
            if key == "start_file_date" or sfd is not None:
                print(f"{YELLOW}{note}{RESET}")
        if key == "start_file_date":
            choices = start_file_choices(opt)
            if choices:
                print(f"{YELLOW}available: {ic_solver.series_span_line(choices)}{RESET}")
        _prompt_leaf(key, spec, ic, path, shown=derived,
                     check=lambda v, p=path, d=derived: ic_date_error(opt, p, v, d))


def prompt_options(opt, mode, desc=None, asked=()):
    """Prompt for every leaf at or below `mode` not filled by engage or already asked."""
    desc = desc or load_descriptions()
    skip = set(opt.get("_engage_filled", ())) | set(IC_DATE_LEAVES) | set(asked)
    print("\n")
    print("Instructions:")
    print(f"-> Default Selected input parameter is given in {GREEN}GREEN{RESET}")
    print(f"-> Warnings and Information are given in {YELLOW}YELLOW{RESET}")
    print(f"-> Errors are given in {RED}RED{RESET}")
    print("-> Valid values (if any) are given in brackets eg. (value1 | value2 | value3)")
    print("-> Enter '?' for any input parameter to get a detailed description")
    filled = set(opt.get("_engage_filled", ()))
    if filled:
        print(f"-> {YELLOW}{len(filled)} option(s) already set by engage are not asked again{RESET}")
    print("\n")
    print("Run Options:")
    for section in desc:
        if section.startswith("_") or section == "version":
            continue
        opt.setdefault(section, {})
        _prompt_section(desc[section], opt[section], mode, section, skip, opt)
    _prompt_ic_dates(opt, desc, mode, filled)
    return opt


def apply_answers(opt, answers):
    """Set nested prompt answers into `opt` as the prompts would (None included)."""
    for path in _leaf_paths(answers):
        src, dst = answers, opt
        parts = path.split(".")
        for part in parts[:-1]:
            src = src[part]
            dst = dst.setdefault(part, {})
        dst[parts[-1]] = copy.deepcopy(src[parts[-1]])
    return opt


def _wizard(opt, mode, desc, hidden=()):
    """The wizard's answers, or None when it cannot run here; Ctrl-C exits 130."""
    import tui_form_waxrun
    try:
        return tui_form_waxrun.run_waxrun_wizard(opt, mode, desc, hidden)
    except tui_form_waxrun._FormCancelled:
        print("[waxrun] wizard cancelled; nothing was written", file=sys.stderr)
        raise SystemExit(130) from None
    except tui_form_waxrun._FormUnavailable as e:
        print(f"waxrun wizard unavailable ({e}); using the linear prompts.")
        return None


def wizard_answers(engage, mode="BASIC", hidden=()):
    """Run the waxrun wizard for engage before it writes anything; the answers, or None."""
    desc = load_descriptions()
    opt = getparams(mode, engage=engage)
    return _wizard(opt, mode, desc, hidden)


def getparams(mode="BASIC", options_path=None, overrides=None, engage=None, interactive=False,
              tui=False, answers=None, asked=None):
    """The options dict: defaults, options JSON, engage payload, overrides, then prompts if interactive."""
    if mode not in LEVELS:
        raise ValueError(f"mode must be one of {LEVELS}, got {mode!r}")
    desc = load_descriptions()
    opt = _defaults(desc)
    opt["api_version"] = API_VERSION
    opt["mode"] = mode
    if options_path:
        loaded = json.loads(Path(options_path).read_text()) if not isinstance(options_path, dict) else options_path
        if "options" in loaded and isinstance(loaded["options"], dict):
            loaded = loaded["options"]
        loaded = {k: v for k, v in loaded.items() if not k.startswith("_")}
        loaded.pop("api_version", None)
        loaded.pop("mode", None)
        _deep_update(opt, loaded, keep_none=True)
    else:
        loaded = {}
    if engage:
        apply_engage_options(opt, engage)
        _note_engage_owned(opt, loaded)
    if overrides:
        _deep_update(opt, overrides)
        opt["_engage_filled"] = sorted(set(opt.get("_engage_filled", ())) | set(_leaf_paths(overrides)))
    apply_waxhome(opt)
    if asked:
        apply_answers(opt, asked)
    if answers is not None:
        apply_answers(opt, answers)
    elif interactive:
        got = _wizard(opt, mode, desc) if tui else None
        if got is None:
            prompt_options(opt, mode, desc, asked=_leaf_paths(asked or {}))
        else:
            apply_answers(opt, got)
    validate_options(opt, desc)
    hpc = opt["simulation"].get("hpc_system", "derecho")
    for section in ("job", "stacks"):
        if hpc not in opt[section]:
            opt[section][hpc] = copy.deepcopy(opt[section]["derecho"])
    if not opt["model"]["data"].get("kaijuhome"):
        opt["model"]["data"]["kaijuhome"] = os.environ.get("KAIJUHOME")
    return opt


def _path_get(d, path):
    for part in path.split("."):
        if not isinstance(d, dict) or part not in d:
            return None
        d = d[part]
    return d


def _note_engage_owned(opt, loaded):
    """NOTE the engage-filled leaves whose -wo value engage replaced; returns them."""
    changed = []
    for path in opt.get("_engage_filled", ()):
        if path == "simulation.run_directory":
            continue
        old = _path_get(loaded, path)
        new = _path_get(opt, path)
        if old in (None, "") or str(old) == str(new):
            continue
        changed.append(f"{path} {old!r} -> {new!r}")
    if changed:
        _note_once("[waxrun] NOTE: engage owns these waxrun options and replaced the -wo JSON's values "
                   "(change them in the engage options): " + ", ".join(changed))
    return changed


WAXHOME_ENV = tree_solver.WAXHOME_ENV
_NOTED = set()      # getparams can run twice per process


def _note_once(line):
    if line in _NOTED:
        return
    _NOTED.add(line)
    _warn(line)


def apply_waxhome(opt):
    """Default simulation.srcroot to $WAXHOME when nothing else set it."""
    sim = opt["simulation"]
    given = str(sim.get("srcroot") or "").strip()
    tree = tree_solver.waxhome()
    if given:
        if tree and os.path.normpath(os.path.abspath(os.path.expanduser(given))) != tree:
            _note_once(f"[waxrun] NOTE: simulation.srcroot is {given} (set explicitly), not "
                       f"${WAXHOME_ENV} ({tree}).")
        return opt
    if tree:
        sim["srcroot"] = tree
    return opt


def _refuse(msg, path, opt):
    hint = " (set by the engage payload; change it there)" if path in set(opt.get("_engage_filled", ())) else ""
    raise SystemExit(f"[waxrun] {msg}{hint}")


def _validate_leaves(desc_sec, opt_sec, prefix, opt, run_dir=None):
    for key, spec in desc_sec.items():
        if key.startswith("_") or not isinstance(spec, dict):
            continue
        path = f"{prefix}.{key}" if prefix else key
        if "LEVEL" in spec:
            if key in opt_sec:
                try:
                    opt_sec[key] = _leaf_value(path, spec, opt_sec[key], run_dir)
                except ValueError as e:
                    _refuse(str(e), path, opt)
        elif isinstance(opt_sec.get(key), dict):
            _validate_leaves(spec, opt_sec[key], path, opt, run_dir)


def _unknown_leaves(desc_sec, opt_sec, prefix=""):
    """Option keys with no options_description.json leaf; the cam passthrough is not checked."""
    out = []
    for key, v in opt_sec.items():
        path = f"{prefix}.{key}" if prefix else key
        if key.startswith("_") or path in INTERNAL_KEYS or path == "cam":
            continue
        spec = desc_sec.get(key)
        if not isinstance(spec, dict):
            out.append(path)
        elif "LEVEL" not in spec and isinstance(v, dict):
            out += _unknown_leaves(spec, v, path)
    return out


def _all_leaf_paths(desc_sec, prefix=""):
    out = []
    for key, spec in desc_sec.items():
        if key.startswith("_") or not isinstance(spec, dict):
            continue
        path = f"{prefix}.{key}" if prefix else key
        out += [path] if "LEVEL" in spec else _all_leaf_paths(spec, path)
    return out


def validate_options(opt, desc=None):
    """Validate and normalise in place the options that arrive from JSON; SystemExit on a bad leaf."""
    desc = desc or load_descriptions()
    hpc = opt["simulation"].get("hpc_system", "derecho")
    try:
        _leaf_value("simulation.hpc_system", desc["simulation"]["hpc_system"], hpc)
    except ValueError as e:
        _refuse(str(e), "simulation.hpc_system", opt)
    queue = ((opt.get("job") or {}).get(hpc) or {}).get("queue")
    if queue is not None:
        machine = opt["simulation"].get("cime_machine") or hpc
        err = queue_error(queue, opt["simulation"].get("srcroot"), machine, hpc)
        if err:
            _refuse(err, f"job.{hpc}.queue", opt)
    run_dir = _run_dir_of(opt)
    for section, spec in desc.items():
        if section.startswith("_") or section == "version" or not isinstance(spec, dict):
            continue
        if isinstance(opt.get(section), dict):
            _validate_leaves(spec, opt[section], section, opt, run_dir)
    sp = opt.get("spinup") or {}
    if None not in (sp.get("min_hours"), sp.get("max_hours")) and sp["min_hours"] > sp["max_hours"]:
        _refuse(f"spinup.min_hours={sp['min_hours']!r} is above spinup.max_hours={sp['max_hours']!r}",
                "spinup.min_hours", opt)
    if sp.get("wall_margin") is not None and sp["wall_margin"] < 1:
        _note_once(f"[waxrun] NOTE: spinup.wall_margin={sp['wall_margin']} is below 1: the spinup walltime "
                   f"estimate has no headroom.")
    unknown = _unknown_leaves({k: v for k, v in desc.items() if k != "version"}, opt)
    if unknown:
        import difflib
        known = _all_leaf_paths(desc)
        lines = []
        for p in unknown:
            close = difflib.get_close_matches(p, known, n=1)
            lines.append(f"{p}" + (f" (did you mean {close[0]}?)" if close else ""))
        raise SystemExit(f"[waxrun] not waxrun option(s) (config/options_description.json): "
                         + "; ".join(lines))
    n = opt.get("nudging") or {}
    if n.get("nudging_dir") and _bool_of(n.get("generate")) is not False:
        _refuse(f"nudging.nudging_dir={n['nudging_dir']!r} is used only when nudging.generate is no",
                "nudging.nudging_dir", opt)
    cp = opt.get("coupled") or {}
    if (cp.get("sourcemods") or "none") != "none" and cp.get("build") == "verify_only":
        _refuse(f"coupled.sourcemods={cp['sourcemods']!r} is used only when coupled.build is build",
                "coupled.sourcemods", opt)
    for section in ("spinup", "coupled"):
        nodes = opt[section].get("pe_nodes")
        whole = isinstance(nodes, int) or (isinstance(nodes, float) and nodes.is_integer())
        if nodes is not None and (not whole or nodes < 1):
            _refuse(f"{section}.pe_nodes={nodes!r} must be a whole number of nodes >= 1", f"{section}.pe_nodes", opt)
    res = opt["simulation"].get("res")
    try:
        ic_solver.expected_grid(res)
    except KeyError as e:
        _refuse(f"simulation.res={res!r}: {e.args[0]}", "simulation.res", opt)
    return opt


def apply_engage_options(opt, engage):
    """Map the engage payload onto the options; values taken are recorded in opt['_engage_filled']."""
    sim, d = opt["simulation"], opt["model"]["data"]
    esim = engage.get("simulation") or {}
    epbs = engage.get("pbs") or {}
    ecpl = engage.get("coupling") or {}
    evolt = engage.get("voltron") or {}
    # deep merge: "waxrun" must not drop engage's own "waccmx" values
    ewax = copy.deepcopy(engage.get("waccmx") or {})
    _deep_update(ewax, copy.deepcopy(engage.get("waxrun") or {}))
    filled = set()

    def put(section, key, value, path):
        section[key] = value
        filled.add(path)

    if esim.get("job_name"):
        put(sim, "job_name", esim["job_name"], "simulation.job_name")
    if esim.get("hpc_system"):
        put(sim, "hpc_system", esim["hpc_system"], "simulation.hpc_system")
    hpc = sim.get("hpc_system", "derecho")
    for section in ("job", "stacks"):
        if hpc not in opt[section]:
            opt[section][hpc] = copy.deepcopy(opt[section]["derecho"])

    run_directory = epbs.get("run_directory") or engage.get("run_directory")
    if run_directory:
        put(sim, "run_directory", run_directory, "simulation.run_directory")
    if esim.get("start_date"):
        put(sim, "coupling_start", _iso(esim["start_date"], "simulation.start_date"), "simulation.coupling_start")
    if esim.get("stop_date"):
        put(sim, "stop_time", _iso(esim["stop_date"], "simulation.stop_date"), "simulation.stop_time")
    if esim.get("segment_duration"):
        put(opt["coupled"], "stephours", _hours(esim["segment_duration"], "segment_duration"),
            "coupled.stephours")
    if ecpl.get("gr_warm_up_time") not in (None, ""):
        put(opt["mage"], "delaycouplinghrs", _hours(ecpl["gr_warm_up_time"], "gr_warm_up_time"),
            "mage.delaycouplinghrs")
    if ecpl.get("gcm_spin_up_time") not in (None, ""):
        put(opt["spinup"], "target_seconds", int(round(float(ecpl["gcm_spin_up_time"]))),
            "spinup.target_seconds")
    if ecpl.get("conda_env"):
        put(d, "conda_env", ecpl["conda_env"], "model.data.conda_env")
    dtc = (evolt.get("coupling") or {}).get("dtCouple")
    if dtc not in (None, ""):
        v = float(dtc)
        put(opt["mage"], "dtcouple", int(v) if v.is_integer() else v, "mage.dtcouple")
    dto = (evolt.get("output") or {}).get("dtOut")
    if dto not in (None, ""):
        v = float(dto)
        put(opt["mage"], "dtout", int(v) if v.is_integer() else v, "mage.dtout")

    proj = epbs.get("account_name") or epbs.get("project_code")
    if proj:
        put(sim, "project", proj, "simulation.project")
    for k in ("queue", "walltime", "job_priority"):
        if epbs.get(k):
            put(opt["job"][hpc], k, epbs[k], f"job.{hpc}.{k}")
    if epbs.get("kaiju_install_directory"):
        put(d, "kaijuhome", os.path.abspath(epbs["kaiju_install_directory"]), "model.data.kaijuhome")
    if epbs.get("kaiju_build_directory") and not d.get("voltron_exe"):
        put(d, "voltron_exe", os.path.join(os.path.abspath(epbs["kaiju_build_directory"]), "bin", "voltron_mpi.x"),
            "model.data.voltron_exe")
    for k in ("cesm_exe", "kaijuhome", "voltron_exe"):
        if ewax.get(k) and not isinstance(ewax[k], dict):
            target, path = (opt["coupled"], "coupled.cesm_exe") if k == "cesm_exe" else (d, f"model.data.{k}")
            put(target, k, ewax[k], path)
    nested = {k: v for k, v in ewax.items() if isinstance(v, dict)}
    if nested:
        _deep_update(opt, nested)
        filled |= set(_leaf_paths(nested))
    egrx = engage.get("grx") or {}
    if egrx.get("coupling_start"):
        put(sim, "coupling_start", _iso(egrx["coupling_start"], "grx.coupling_start"), "simulation.coupling_start")
    if egrx.get("stop"):
        put(sim, "stop_time", _iso(egrx["stop"], "grx.stop"), "simulation.stop_time")
    if egrx.get("delaycouplinghrs") not in (None, ""):
        put(opt["mage"], "delaycouplinghrs", int(egrx["delaycouplinghrs"]), "mage.delaycouplinghrs")
    if egrx.get("stephours") not in (None, ""):
        put(opt["coupled"], "stephours", int(egrx["stephours"]), "coupled.stephours")
    if egrx.get("target_start") and sim.get("coupling_start"):
        secs = (_dt(sim["coupling_start"]) - _dt(egrx["target_start"])).total_seconds()
        put(opt["spinup"], "target_seconds", int(round(secs)), "spinup.target_seconds")
    if egrx.get("project"):
        put(sim, "project", egrx["project"], "simulation.project")
    if egrx.get("job_priority"):
        put(opt["job"][hpc], "job_priority", egrx["job_priority"], f"job.{hpc}.job_priority")
    if egrx.get("run_directory"):
        put(sim, "run_directory", egrx["run_directory"], "simulation.run_directory")
    opt["_engage"] = engage
    opt["_engage_filled"] = sorted(set(opt.get("_engage_filled", ())) | filled)
    return opt


def _layout_of(opt):
    """layout() of the options; always recomputed, never stored."""
    return layout(opt["simulation"]["job_name"], opt["simulation"].get("run_directory") or ".")


def _window(opt):
    """The coupled window (start, stop, delay, stephours, MAGE t0, segment counts), checked."""
    sim, hpc = opt["simulation"], opt["simulation"].get("hpc_system", "derecho")
    cs_ = _run_dt(sim["coupling_start"], "simulation.coupling_start")
    stop = _run_dt(sim["stop_time"], "simulation.stop_time")
    if stop <= cs_:
        raise SystemExit(f"stop_time {stop} is not after coupling_start {cs_}")
    err = segment_solver.coupling_start_error(
        cs_, (opt.get("mage") or {}).get("dtcouple"),
        (opt.get("spinup") or {}).get("dtime") or segment_solver.SPINUP_DTIME_DEFAULT)
    if err:
        raise SystemExit(f"[waxrun] {err}")
    delay_h = _whole_hours(opt["mage"]["delaycouplinghrs"], "mage.delaycouplinghrs")
    step_h = _whole_hours(opt["coupled"]["stephours"], "coupled.stephours")
    window_s = (stop - cs_).total_seconds()
    wall_h = _secs(opt["job"][hpc]["walltime"]) / 3600.0
    rate = grx_rate(hpc)
    if rate.get("note"):
        _warn(f"[waxrun] {rate['note']}")
    if step_h < 1:
        raise SystemExit(f"coupled.stephours={step_h} (engage simulation.segment_duration / 3600) must "
                         f"be at least 1 h")
    if delay_h < 1:
        raise SystemExit(f"mage.delaycouplinghrs={delay_h} (engage coupling.gr_warm_up_time / 3600) must be "
                         f"at least 1 h and at most coupled.stephours={step_h} h")
    if window_s % (step_h * 3600):
        hint = ""
        if window_s % 3600:
            hint = (f"; the window is not a whole number of hours: stop_time must share coupling_start's "
                    f"minute and second ({cs_:%H:%M:%S})")
        raise SystemExit(f"coupled.stephours={step_h} does not divide the coupled window "
                         f"({window_s / 3600:g} h from {cs_:%Y-%m-%d %H:%M} to {stop:%Y-%m-%d %H:%M}); "
                         f"valid values (hours): {segment_solver.valid_stephours(window_s / 3600, delay_h, wall_h, rate=rate['rate'], margin=rate['margin'])}"
                         + hint)
    t0 = cs_ - timedelta(hours=delay_h)
    ic, iu, maxhours = segment_solver.compute_segments(t0, stop, step_h, delay_h)
    if ic != 1:
        raise SystemExit(f"gr_warm_up_time must be <= segment_duration (delay {delay_h} h, "
                         f"segment {step_h} h gives ic={ic})")
    if not segment_solver.segment_fits_walltime(step_h, wall_h, rate["rate"], rate["margin"]):
        _warn(f"[waxrun] WARNING: a {step_h} h coupled segment is expected to need "
              f"{_hms(segment_solver.segment_walltime_hours(step_h, rate['rate'], rate['margin']) * 3600)}; "
              f"walltime is {opt['job'][hpc]['walltime']}")
    return {"coupling_start": cs_, "stop": stop, "delay_h": delay_h, "step_h": step_h, "t0": t0,
            "ic": ic, "iu": iu, "maxhours": maxhours, "window_h": window_s / 3600.0, "wall_h": wall_h}


def _srcroot(opt, kind):
    """The tree a phase is built from: simulation.srcroot, or spinup.srcroot for a spinup if set."""
    base = (opt["simulation"].get("srcroot") or "").strip()
    override = (opt["spinup"].get("srcroot") or "").strip() if kind == "spinup" else ""
    src = override or base
    if not src:
        where = f"in {os.path.basename(_OPTIONS_FILE)}" if _OPTIONS_FILE else "with --srcroot"
        raise SystemExit(f"ERROR: simulation.srcroot is required: set it {where} or export ${WAXHOME_ENV} "
                         f"(clone steps: {CLONE_STEPS}).")
    return os.path.normpath(os.path.abspath(os.path.expanduser(src)))


def _tree(opt, srcroot, kind, sourcemods, single_tree):
    """tree_solver.tree_spec for one phase, cached; refuses on any FAIL."""
    if opt.get("_skip_tree_check"):
        _warn(f"[waxrun] WARNING: --skip-tree-check: {srcroot} is not checked")
        return None
    pins = opt["simulation"].get("tree_pins")
    pins = pins if not pins or os.path.isabs(pins) else str(HERE / pins)
    key = (srcroot, kind, sourcemods, single_tree, pins)
    if key not in _TREE_CACHE:
        _TREE_CACHE[key] = tree_solver.tree_spec(srcroot, kind, pins=tree_solver.load_pins(pins),
                                                 sourcemods=sourcemods, single_tree=single_tree)
    spec = copy.deepcopy(_TREE_CACHE[key])
    for f in spec["findings"]:
        if f["level"] != "OK":
            _warn(f"[waxrun] tree {f['level']:4s} {f['check']}: {f['detail'].splitlines()[0]}")
    if not spec["ok"]:
        key = "spinup.srcroot" if kind == "spinup" and (opt["spinup"].get("srcroot") or "").strip() \
            else "simulation.srcroot"
        origin = (f"${WAXHOME_ENV}" if tree_solver.waxhome() == srcroot
                  else os.path.basename(_OPTIONS_FILE) if _OPTIONS_FILE else None)
        fails = "; ".join(f"{f['check']}: {f['detail'].splitlines()[0]}"
                          for f in spec["findings"] if f["level"] == "FAIL")
        raise SystemExit(f"ERROR: {key} {srcroot}{f' ({origin})' if origin else ''} cannot build the {kind} "
                         f"case: {fails} (clone steps: {CLONE_STEPS}).")
    return spec


def _tree_module_script(srcroot, machine):
    """The module script case.setup derives from the tree's config_machines.xml, or None."""
    path = os.path.join(srcroot, "ccs_config", "machines", machine, "config_machines.xml")
    try:
        return case_solver.module_script(case_solver.module_commands(path))
    except (OSError, ET.ParseError, ValueError) as e:
        _warn(f"[waxrun] WARNING: cannot read the {machine} module_system in {path} ({e}); "
              "the ldd checks run without module loads")
        return None


def _tree_module_system_xml(srcroot, machine):
    """The raw <module_system> text of the tree's config_machines.xml (dry-run stand-in)."""
    path = os.path.join(srcroot, "ccs_config", "machines", machine, "config_machines.xml")
    try:
        m = case_solver._MODSYS.search(Path(path).read_text())
    except OSError:
        m = None
    return m.group(0).strip() if m else '<module_system type="module"></module_system>'


# namelist numbers and logicals, written unquoted
_NL_NUMBER = re.compile(r"^(\d+\*)?[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")
_NL_LOGICAL = re.compile(r"^(\d+\*)?(\.(true|false|t|f)\.|t|f|true|false)$", re.I)


# one quoted character constant, or a comma list of them, as given
_NL_QUOTED_LIST = re.compile(r"""^\s*('(?:[^']|'')*'|"(?:[^"]|"")*")(\s*,\s*('(?:[^']|'')*'|"(?:[^"]|"")*"))*\s*$""")
_NL_INTEGER = re.compile(r"^(\d+\*)?[+-]?\d+$")


def _nl_logical(item):
    """'.true.'/'.false.' for a logical spelling, else None (Build::Namelist reads only these)."""
    m = _NL_LOGICAL.match(str(item).strip())
    if not m:
        return None
    word = m.group(2).lower().strip(".")
    return (m.group(1) or "") + (".true." if word in ("t", "true") else ".false.")


def _nl_quote(text):
    """Quote as a Fortran string; "..." when the text holds a ' (CAM splits a doubled '')."""
    if "'" in text:
        return '"' + text.replace('"', '""') + '"'
    return "'" + text + "'"


def _render_override(v):
    """A cam.* value as namelist text by its shape, used when namelist_definition.xml is unreadable."""
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, (list, tuple)):
        return ", ".join(_render_override(x) for x in v)
    if not isinstance(v, str):
        return str(v)
    if _QUOTED_RHS.match(v):
        return v
    items = [i.strip() for i in v.split(",")]
    if all(_NL_NUMBER.match(i) or _NL_LOGICAL.match(i) for i in items):
        if not any(_NL_LOGICAL.match(i) for i in items):
            return v
        return ", ".join(_nl_logical(i) or i for i in items)
    return _nl_quote(v.strip())


_NL_TYPES_CACHE = {}
_NL_ENTRY = re.compile(r"<entry\s([^>]*)>", re.S)
_NL_ATTR = re.compile(r'([A-Za-z_]+)\s*=\s*"([^"]*)"')
_NL_TYPE = re.compile(r"^(logical|integer|real|char)(?:\*(\d+))?(\([^)]*\))?$")


def namelist_types(srcroot):
    """{variable: {type, len, array, valids}} from the tree's namelist_definition.xml ({} if unreadable)."""
    if not srcroot:
        return {}
    path = os.path.join(str(srcroot), "bld", "namelist_files", "namelist_definition.xml")
    if path not in _NL_TYPES_CACHE:
        try:
            text = Path(path).read_text()
        except OSError:
            text = ""
        out = {}
        for m in _NL_ENTRY.finditer(text):
            attrs = dict(_NL_ATTR.findall(m.group(1)))
            t = _NL_TYPE.match(attrs.get("type", "").strip())
            if not attrs.get("id") or not t:
                continue
            out[attrs["id"].lower()] = {
                "type": t.group(1), "len": int(t.group(2)) if t.group(2) else None,
                "array": bool(t.group(3)),
                "valids": [v.strip() for v in attrs.get("valid_values", "").split(",") if v.strip()]}
        _NL_TYPES_CACHE[path] = out
    return _NL_TYPES_CACHE[path]


def _nl_valid(leaf, text, valids):
    """`text` in the XML's valid values (any case), in the XML's own spelling."""
    if not valids:
        return text
    for v in valids:
        if v.lower() == text.lower():
            return v
    raise ValueError(f"{leaf}={text!r} is not valid ({'|'.join(valids)})")


def _nl_element(leaf, v, defn):
    kind, valids = defn["type"], defn.get("valids")
    if kind == "logical":
        if isinstance(v, str) and _nl_logical(v):
            return _nl_logical(v)
        b = _bool_of(v)
        if b is None:
            raise ValueError(f"{leaf}={v!r} is not a logical (true/false)")
        return ".true." if b else ".false."
    if kind in ("integer", "real"):
        if isinstance(v, bool) or not isinstance(v, (int, float, str)):
            raise ValueError(f"{leaf}={v!r} is not {'an integer' if kind == 'integer' else 'a number'}")
        if kind == "integer":
            if isinstance(v, float) or (isinstance(v, str) and not _NL_INTEGER.match(v.strip())):
                raise ValueError(f"{leaf}={v!r} is not an integer")
        elif isinstance(v, float) and not math.isfinite(v):
            raise ValueError(f"{leaf}={v!r} is not a finite number")
        elif isinstance(v, str) and not _NL_NUMBER.match(v.strip()):
            raise ValueError(f"{leaf}={v!r} is not a number")
        text = str(v).strip()
        return _nl_valid(leaf, text, valids) if valids else text
    if isinstance(v, bool) or not isinstance(v, (str, int, float)):
        raise ValueError(f"{leaf}={v!r} is not a character value")
    text = str(v).strip() if not isinstance(v, str) else v.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"" and _NL_QUOTED_LIST.match(text) \
            and len(_NL_STRING.findall(text)) == 1:
        inner = text[1:-1].replace(text[0] * 2, text[0])
        return _nl_quote(_nl_valid(leaf, inner, valids)) if valids else text
    return _nl_quote(_nl_valid(leaf, text, valids))


def format_nl_value(name, value, defn):
    """The user_nl_cam right-hand side of cam.<name>, typed by its namelist_definition.xml entry."""
    leaf = f"cam.{name}"
    if isinstance(value, dict) or value is None:
        raise ValueError(f"{leaf}={value!r} is not a namelist value")
    if isinstance(value, (list, tuple)):
        if not defn["array"]:
            raise ValueError(f"{leaf} takes one value, got the list {list(value)!r}")
        if not value:
            raise ValueError(f"{leaf} is an empty list")
        items = list(value)
    elif defn["type"] == "char" and isinstance(value, str) and _NL_QUOTED_LIST.match(value):
        items = [m.group(0) for m in _NL_STRING.finditer(value)]
        if len(items) > 1 and not defn["array"]:
            raise ValueError(f"{leaf} takes one value, got {value!r}")
    elif defn["type"] != "char" and isinstance(value, str) and defn["array"] and "," in value:
        items = [i.strip() for i in value.split(",")]
    else:
        items = [value]
    if any(isinstance(i, (list, tuple, dict)) for i in items):
        raise ValueError(f"{leaf}={value!r} is not a list of values")
    return ", ".join(_nl_element(leaf, i, defn) for i in items)


def _cam_overrides(opt, srcroot=None):
    """{namelist key: rendered value} of the cam.* overrides; an unknown key is refused."""
    types = namelist_types(srcroot)
    out = {}
    for k, v in (opt.get("cam") or {}).items():
        if k.startswith("_") or v is None:
            continue
        if k.lower() in DERIVED_NL_KEYS:
            raise SystemExit(f"cam.{k} is derived by waxrun and cannot be overridden")
        if not types:
            out[k] = _render_override(v)
            continue
        defn = types.get(k.lower())
        if defn is None:
            raise SystemExit(f"[waxrun] cam.{k} is not a CAM namelist variable (namelist_definition.xml)")
        try:
            out[k] = format_nl_value(k, v, defn)
        except ValueError as e:
            raise SystemExit(f"[waxrun] {e}") from None
    return out


INITHIST_VALIDS = ("NONE", "6-HOURLY", "DAILY", "MONTHLY", "YEARLY", "CAMIOP", "ENDOFRUN")


def _inithist(value, srcroot):
    """coupled.inithist in the spelling the tree's namelist_definition.xml accepts."""
    valids = (namelist_types(srcroot).get("inithist") or {}).get("valids") or list(INITHIST_VALIDS)
    try:
        return _nl_valid("coupled.inithist", str(value).strip().strip("'\""), valids)
    except ValueError as e:
        raise SystemExit(f"[waxrun] {e}") from None


def _nl_provenance(nl_patches):
    """{namelist key: the cam.* leaf that set it}, for the forcing gate's messages."""
    return {k.lower(): f"cam.{k} (options JSON)" for k in nl_patches}


def spinup_walltime(stop_seconds, wall_s_per_simday, margin, queue_max):
    """(walltime estimate rounded up to 30 min, number of chunks) for a spinup of stop_seconds."""
    est_s = float(wall_s_per_simday) * (stop_seconds / 86400.0) * float(margin)
    if est_s <= 0:
        raise SystemExit(f"[waxrun] the spinup walltime estimate is {est_s:g} s: spinup.wall_s_per_simday "
                         f"and spinup.wall_margin must be > 0")
    est_s = int(math.ceil(est_s / 1800.0) * 1800)
    max_s = _secs(queue_max)
    k = max(1, int(math.ceil(est_s / max_s)))
    return _hms(est_s), k


def spinup_chunks(stop_n, k):
    """Split stop_n into at most k chunks of ceil(stop_n / k), the remainder (>= 1) last."""
    stop_n = int(stop_n)
    c = int(math.ceil(stop_n / k)) if k > 1 else stop_n
    k = int(math.ceil(stop_n / c)) if c > 0 else 1
    if k <= 1:
        return [{"stop_n": stop_n, "remainder": False}]
    chunks = [{"stop_n": c, "remainder": False} for _ in range(k - 1)]
    chunks.append({"stop_n": int(stop_n) - c * (k - 1), "remainder": True})
    return chunks


def _forcing_precheck(nl_text, window, source, strict=True, provenance=None):
    """Check that the forcing files named in user_nl_cam text exist and cover the window."""
    start, end = forcing_solver._dt(window[0]), forcing_solver._dt(window[1])
    if end <= start:
        raise ValueError(f"forcing window end {end} is not after start {start}")
    files = forcing_solver.namelist_files(nl_text)
    provenance = provenance or {}
    report = {"ok": True, "source": str(source), "sst": None,
              "window": [start.strftime("%Y-%m-%d %H:%M:%S"), end.strftime("%Y-%m-%d %H:%M:%S")],
              "findings": []}
    for key in forcing_solver.FORCING_KEYS:
        short, path = forcing_solver.SHORT[key], files[key]
        where = provenance.get(key, str(source))
        if path is None or (not path and key in forcing_solver.OPTIONAL_KEYS):
            report[short] = None
            # present but unparsed means unquoted: buildnml would reject it
            if path is None and _nl_key_present(nl_text, key):
                report["findings"].append({"level": "error", "key": key, "file": None,
                                           "message": f"{key} is set in {where} but its value is not a "
                                                      f"quoted Fortran string, so it cannot be checked "
                                                      f"and ./preview_namelists will reject the line"})
                continue
            report["findings"].append({"level": "info" if key in forcing_solver.OPTIONAL_KEYS else "warning",
                                       "key": key, "file": None,
                                       "message": f"{key} is not set in {where} (a compset default may apply; "
                                                  f"the CaseDocs/atm_in gate in genCASE decides)"})
            continue
        if "$" in path:
            expanded = os.path.expandvars(path)
            if "$" in expanded:
                report[short] = {"file": path, "first": None, "last": None, "n": 0}
                report["findings"].append({"level": "info", "key": key, "file": path,
                                           "message": f"{key} = {path} (from {where}) names an unset "
                                                      f"variable; checked in genCASE (CaseDocs/atm_in)"})
                continue
            path = expanded
        if not os.path.isfile(path):
            report[short] = {"file": path, "first": None, "last": None, "n": 0}
            report["findings"].append({"level": "error", "key": key, "file": path,
                                       "message": f"{key} = {path} (from {where}) does not exist"})
            continue
        first, last, n = forcing_solver.file_span(path)
        report[short] = {"file": path, "first": first.strftime("%Y-%m-%d %H:%M:%S"),
                         "last": last.strftime("%Y-%m-%d %H:%M:%S"), "n": n}
        if start < first or end > last:
            report["findings"].append({
                "level": "error", "key": key, "file": path,
                "message": f"{key} = {os.path.basename(path)} covers {first:%Y-%m-%d %H:%M} .. "
                           f"{last:%Y-%m-%d %H:%M}; the run needs {start:%Y-%m-%d %H:%M} .. "
                           f"{end:%Y-%m-%d %H:%M} (last date on file {last:%Y-%m-%d}, "
                           f"requested end {end:%Y-%m-%d %H:%M})"})
        else:
            report["findings"].append({"level": "ok", "key": key, "file": path,
                                       "message": f"{key} covers {first:%Y-%m-%d} .. {last:%Y-%m-%d}"})
    errors = [f for f in report["findings"] if f["level"] == "error"]
    report["ok"] = not errors
    if errors and strict:
        raise SystemExit("forcing coverage gate FAILED for window "
                         f"{report['window'][0]} .. {report['window'][1]}:\n  "
                         + "\n  ".join(e["message"] for e in errors)
                         + "\nCAM does not extrapolate these files: use ones that cover the window.")
    return report


def _nl_key_present(nl_text, key):
    """True when `key` is assigned on an uncommented line, whatever the value."""
    pat = re.compile(r"^\s*" + re.escape(key) + r"\s*=", re.I)
    return any(pat.match(l) for l in nl_text.splitlines() if not l.lstrip().startswith("!"))


def _check_generator(n):
    """Refuse a nudging.gen_workers that does not fit the node's cores or memory."""
    workers = n["gen_workers"]
    if isinstance(workers, bool) or not isinstance(workers, int):
        raise SystemExit(f"[waxrun] nudging.gen_workers={workers!r} must be a whole number of workers")
    try:
        nudge_solver.gen_workers(workers)
    except ValueError as e:
        raise SystemExit(f"[waxrun] {e}") from None
    cap = float(nudge_solver.load_config()["pbs"].get("max_mem_gb") or 0)
    mem = nudge_solver.gen_mem_gb(workers)
    if cap and mem > cap:
        raise SystemExit(f"[waxrun] nudging.gen_workers={workers} needs {mem}GB, above the node limit of "
                         f"{int(cap)}GB (config/nudging_gen.json)")


def nudging_mode(n):
    """(enabled, generate, dir); nudging is off when generate is no and no directory is given."""
    generate = _as_bool(n["generate"], "nudging.generate")
    ndir = n.get("nudging_dir") or None
    return bool(generate or ndir), generate, (None if generate else ndir)


def _nudging(opt, lay, beg_dt, stop, dtime):
    """(nudging block of a spec, Nudge_* namelist values); Nudge_Path is set even when nudging is off."""
    n = opt["nudging"]
    enabled, generate, ndir = nudging_mode(n)
    _check_generator(n)
    tpd, template = nudge_solver.TIMES_PER_DAY, nudge_solver.FILE_TEMPLATE
    beg, end = nudge_solver.window(beg_dt, stop)
    assembled = _expand(n["assembled_dir"], lay)
    search = [ndir] if ndir else _expand(list(n["search_dirs"]), lay)
    outdir = os.path.join(lay["shared_dir"], "met_data", "gen", "{YYYY}")
    if generate and outdir not in search:
        search.append(outdir)
    block = {"enabled": enabled, "path": assembled.rstrip("/") + "/", "beg": beg.isoformat(),
             "end": end.isoformat(), "generate": generate, "nudging_dir": ndir,
             "gen_workers": int(n["gen_workers"]), "gen_dir": lay["nudging_gen_dir"],
             "gen_outdir": outdir, "rda_input": n["rda_input"], "gen_job": None, "gen_spans": [],
             "assembled_dir": assembled, "search_dirs": search, "times_per_day": tpd, "file_template": template,
             "vwin_lindex": int(n["vwin_lindex"]), "vwin_hindex": int(n["vwin_hindex"]),
             "slots": 0, "missing": 0, "coverage": {}}
    if enabled:
        slots = nudge_solver.slots(beg, end, tpd)
        cov = nudge_solver.coverage(search, slots, template)
        spans = nudge_solver.missing_spans(cov)
        block["slots"] = len(slots)
        block["missing"] = len(nudge_solver.missing(cov))
        block["coverage"] = {s.isoformat(): p for s, p in cov.items()}
        if spans and not generate:
            gaps = nudge_solver.missing(cov)
            raise SystemExit(f"[waxrun] nudging: {ndir} does not cover the run ({beg} .. {end}, "
                             f"{tpd} files a day): {len(gaps)} of {len(slots)} slot(s) missing or "
                             f"incomplete, first {gaps[0]:%Y-%m-%d %H:%M}, last {gaps[-1]:%Y-%m-%d %H:%M}. "
                             "Add the files, or answer yes to 'Generate nudging data?'")
        if spans:
            block["gen_spans"] = [[a.isoformat(), b.isoformat()] for a, b in spans]
    nl = nudge_solver.namelist(dtime, beg, end, enabled, block["path"], False, tpd, template)
    nl["Nudge_Vwin_Lindex"] = int(n["vwin_lindex"])
    nl["Nudge_Vwin_Hindex"] = int(n["vwin_hindex"])
    return block, nl


def _nudge_vars(nl):
    """Template variables from a nudge_solver.namelist dict (paths unquoted)."""
    out = {"nudge_model": nl["Nudge_Model"], "nudge_path": nl["Nudge_Path"].strip("'"),
           "model_times_per_day": int(nl["Model_Times_Per_day"]),
           "nudge_vwin_lindex": int(nl["Nudge_Vwin_Lindex"]),
           "nudge_vwin_hindex": int(nl["Nudge_Vwin_Hindex"])}
    for edge in ("Beg", "End"):
        for part in ("Year", "Month", "Day"):
            out[f"nudge_{edge.lower()}_{part.lower()}"] = int(nl[f"Nudge_{edge}_{part}"])
    return out


def _ic(opt, lay, w):
    """The IC block of the spinup spec (restart pick, re-date, checks, planned strip); reads only."""
    ic, sp, res = opt["ic"], opt["spinup"], opt["simulation"]["res"]
    target = w["coupling_start"] - timedelta(seconds=int(sp["target_seconds"]))
    min_h, max_h = float(sp["min_hours"]), float(sp["max_hours"])
    source = ic["start_file_dir"]
    sfd, start = ic.get("start_file_date"), ic.get("spinup_start")
    rejected = []
    if ic.get("pick"):
        if not os.path.isfile(ic["pick"]):
            raise SystemExit(f"[waxrun] ic.pick {ic['pick']} (the replay pin) does not exist: restore the "
                             f"file, or point options.ic.pick at the restart to use and replay with --clobber")
        pinned = ic_solver.candidate(ic["pick"])
        try:
            cands = ic_solver.enumerate_candidates(source)
        except FileNotFoundError:
            cands = [pinned]
        if not any(os.path.realpath(c.path) == os.path.realpath(pinned.path) for c in cands):
            cands = sorted(cands + [pinned], key=lambda c: (c.year, c.month, c.day, c.tod))
        if sfd and not ic_solver.matches_date(pinned, sfd):
            raise SystemExit(f"[waxrun] ic.pick {pinned.path} is not a restart of start file date {sfd}: "
                             f"edit or clear options.ic.pick or options.ic.start_file_date")
        # with a date set, the bounds apply to the spinup start instead (below)
        pick = (pinned if sfd or start else
                ic_solver.select_restart([pinned], w["coupling_start"], target, min_h, max_h))
        info = _inspect(pick.path)
        try:
            ic_solver.assert_grid(info, res)
        except ValueError as e:
            raise SystemExit(f"[waxrun] ic.pick: {e}") from None
        # reproduce the rejected list when an unchanged replay still lands on the pin
        try:
            again, _, again_rejected = _native_pick(cands, w["coupling_start"], target, min_h, max_h,
                                                    sfd, res, warn=False)
            if os.path.realpath(again.path) == os.path.realpath(pinned.path):
                rejected = again_rejected
        except SystemExit:
            pass
    else:
        cands = _candidates(source)
        pick, info, rejected = _native_pick(cands, w["coupling_start"], target, min_h, max_h,
                                            sfd, res)
    try:
        ic_solver.assert_required(info)
    except ValueError as e:
        raise SystemExit(f"[waxrun] ic: {e}") from None
    if start:
        spinup_from = _run_dt(start, "ic.spinup_start")
    else:
        spinup_from, redate = default_spinup(pick, w["coupling_start"], opt)
        if "default moved" in redate:
            _note_once(f"[waxrun] ic: start file {ic_solver.stamp_iso(pick)} {redate}")
    err = ic_solver.spinup_start_error(spinup_from, w["coupling_start"], min_h, max_h,
                                       _spinup_dtime(opt) if start else None)
    if err:
        raise SystemExit(f"[waxrun] {err}")
    run_startdate, start_tod, stop_seconds = ic_solver.redate(pick, w["coupling_start"], spinup_from)
    redated = datetime.strptime(run_startdate, "%Y-%m-%d") + timedelta(seconds=start_tod)
    phys_ok = ic.get("phys_ok") or ic_solver.PHYS_OK
    constituents = _constituents(_srcroot(opt, "spinup"), opt["simulation"]["compset"],
                                 sp["cam_config_opts_append"])
    offenders = sorted(ic_solver.physgrid_offenders(info, phys_ok, constituents))
    extra = sorted(set(offenders) - ic_solver.physgrid_offenders(info, phys_ok))
    if extra:
        _warn(f"[waxrun] ic: {', '.join(extra)} on ncol in {os.path.basename(pick.path)} are advected "
              f"by this tree's chemistry; stripped too")
    shared_dir = _expand(ic["shared_dir"], lay)
    # from this run's shared_dir, not options ic.stripped (a replay may use another run dir)
    stripped = os.path.join(shared_dir, ic_solver.stripped_name(pick.path)) if offenders else None
    if stripped and os.path.abspath(stripped) == os.path.abspath(pick.path):
        raise SystemExit(f"refusing to strip {pick.path} onto itself; set ic.shared_dir elsewhere")
    reference_diff = None
    if ic.get("reference"):
        try:
            ref = ic_solver.inspect_ic(ic["reference"]).varset
            final = info.varset - set(offenders)
            reference_diff = {"+": sorted(final - ref), "-": sorted(ref - final)}
        except OSError as e:
            reference_diff = {"error": str(e)}
    block = {
        "start_file_dir": source,
        "candidates": [{"stamp": ic_solver.stamp_of(c), "path": c.path} for c in cands],
        "pick": pick.path, "pick_stamp": ic_solver.stamp_of(pick), "pick_redated": redated.isoformat(),
        "target_start": target.isoformat(), "distance_hours": abs((redated - target).total_seconds()) / 3600.0,
        "start_file_date": sfd, "spinup_start": start, "min_hours": min_h, "max_hours": max_h,
        "stripped": stripped, "offenders": offenders, "shared_dir": shared_dir,
        "grid": {k: info.dims.get(k) for k in ("ncol", "ncol_d", "lev", "ilev")},
        "nvars": len(info.varset) - len(offenders), "date": info.date, "datesec": info.datesec,
        "ncdata": stripped or pick.path, "reference": ic.get("reference"), "reference_diff": reference_diff,
        "strip_cmd": (_planned_strip(["-O", "-x", "-v", ",".join(offenders), pick.path, stripped])
                      if offenders else None),
        "run_startdate": run_startdate, "start_tod": start_tod, "stop_seconds": stop_seconds,
        "rejected": rejected,
    }
    return block


def _planned_strip(args):
    try:
        return ic_solver._ncks_cmd(args)
    except ic_solver.NcksNotFound as e:
        raise SystemExit(f"[waxrun] {e}") from None


_CONSTITUENT_WARNED = set()


def _constituents(srcroot, compset, opts_append):
    """The spinup tree's transported constituents, or None (warned once)."""
    cons, why = ic_solver.constituents_for(srcroot, compset, opts_append)
    if cons is None and why not in _CONSTITUENT_WARNED:
        _CONSTITUENT_WARNED.add(why)
        _warn(f"[waxrun] WARNING: no constituent list ({why}); stripping every 3-D ncol variable "
              f"not in ic.phys_ok")
    return cons


def _assert_nudge_vertical(nudging, ic_path, what):
    """Refuse nudging files whose vertical grid differs from the IC's (CAM does not interpolate)."""
    if not nudging.get("enabled"):
        return
    slot = next((p for _, p in sorted(nudging.get("coverage", {}).items()) if p), None)
    if slot is None:
        return          # all slots generated later; the generator job validates them
    bad = nudge_solver.vertical_mismatch(slot, ic_path)
    if bad:
        raise SystemExit(f"[waxrun] nudging: {slot} and the {what} {ic_path} have different vertical "
                         f"grids ({'; '.join(bad)}); the nudging files need the IC's levels")


_IC_INFO = {}


def _inspect(path):
    key = os.path.realpath(path)
    if key not in _IC_INFO:
        _IC_INFO[key] = ic_solver.inspect_ic(path)
    return _IC_INFO[key]


def _candidates(source):
    try:
        return ic_solver.enumerate_candidates(source)
    except FileNotFoundError as e:
        raise SystemExit(f"[waxrun] {e}") from None


def _native_pick(cands, coupling_start, target, min_h, max_h, start_file_date, res, warn=True):
    """(pick, ICInfo, rejected): the restart pick, skipping candidates not on the native `res` grid."""
    try:
        ic_solver.expected_grid(res)
    except KeyError as e:
        raise SystemExit(f"[waxrun] simulation.res={res!r}: {e.args[0]}") from None
    if start_file_date is not None and not ic_solver.START_FILE_DATE_RE.fullmatch(str(start_file_date).strip()):
        raise SystemExit(f"[waxrun] ic.start_file_date={start_file_date!r} must be "
                         f"{ic_solver.START_FILE_DATE_FORMS}")
    pool, rejected = list(cands), []

    def note():
        return ("" if not rejected else
                f"\nNot native {res} (rejected before the pick):\n  "
                + "\n  ".join(f"{r['stamp']}  {r['path']}: {r['why']}" for r in rejected))

    while True:
        try:
            pick = ic_solver.select_restart(pool, coupling_start, target, min_h, max_h, start_file_date)
        except SystemExit as e:
            raise SystemExit(f"{e.code}{note()}") from None
        info = _inspect(pick.path)
        try:
            ic_solver.assert_grid(info, res)
            if rejected and warn:
                _warn(f"[waxrun] ic: {len(rejected)} non-native candidate(s) skipped before the pick "
                      f"({', '.join(os.path.basename(r['path']) for r in rejected)})")
            return pick, info, rejected
        except ValueError as e:
            why = "; ".join(l.strip() for l in str(e).splitlines()[1:]) or str(e)
            rejected.append({"stamp": ic_solver.stamp_of(pick), "path": pick.path, "why": why})
            pool = [c for c in pool if c.path != pick.path]
            if not pool:
                raise SystemExit(f"[waxrun] no native {res} restart in ic.start_file_dir{note()}") from None


def _build_policy(section, kind):
    """The spec 'build' dict; the spinup always builds."""
    if kind == "spinup":
        return {"mode": "build", "clean_all": False, "exeroot": None}
    mode = str(section["build"])
    if mode == "build":
        return {"mode": "build", "clean_all": True, "ldd_must_be_zero": True}
    if mode == "verify_only":
        exe = section.get("cesm_exe")
        if not exe:
            raise SystemExit("coupled.build=verify_only needs coupled.cesm_exe")
        exe = os.path.abspath(exe)
        if not os.path.exists(exe):
            raise SystemExit(f"[waxrun] coupled.cesm_exe={exe!r} does not exist (coupled.build=verify_only)")
        if os.path.isdir(exe) and os.path.exists(os.path.join(exe, "env_case.xml")):
            return {"mode": "verify_only", "clean_all": False, "ldd_must_be_zero": True,
                    "cesm_exe": os.path.join(case_solver.xmlquery(exe, "EXEROOT")["EXEROOT"], "cesm.exe"),
                    "source_caseroot": exe}
        _warn(f"[waxrun] WARNING: coupled.cesm_exe={exe} is not a caseroot; the build step needs the "
              "built case's caseroot and will refuse")
        return {"mode": "verify_only", "clean_all": False, "ldd_must_be_zero": True, "cesm_exe": exe,
                "source_caseroot": None}
    raise SystemExit(f"coupled.build must be 'build' or 'verify_only', got {mode!r}")


def _identity(opt, srcroot, caseroot, cor, account):
    """The create_newcase identity fields of a spec."""
    sim = opt["simulation"]
    return {
        "srcroot": srcroot, "caseroot": caseroot, "cime_output_root": cor, "project": account,
        "compset": sim["compset"], "res": sim["res"],
        "machine": sim.get("cime_machine") or sim.get("hpc_system", "derecho"),
        "cime_machine": sim.get("cime_machine"),
        "compiler": sim["compiler"], "mpilib": sim["mpilib"],
        "run_unsupported": _as_bool(sim["run_unsupported"], "simulation.run_unsupported"),
    }


def _pe_spec(nodes, leaf):
    nodes = int(nodes)
    if nodes < 1:
        raise SystemExit(f"[waxrun] {leaf}={nodes} must be >= 1 node")
    return {"ntasks_atm": nodes * TASKS_PER_NODE, "nodes": nodes, "tasks_per_node": TASKS_PER_NODE}


def _dtime(value, leaf, coupled=False):
    """A CAM timestep (s) dividing 3600 and 86400; coupled=True also applies the dtCouple rule."""
    dtime = int(value)
    if dtime <= 0 or 3600 % dtime or 86400 % dtime:
        raise SystemExit(f"{leaf}={value} s must divide 3600 s and 86400 s (ATM_NCPL and "
                         f"Model_Times_Per_day are whole numbers)")
    err = segment_solver.dtcouple_error(dtime, leaf) if coupled else None
    if err:
        raise SystemExit(f"[waxrun] {err}")
    return dtime


_NLEV = re.compile(r"-nlev\s+(\d+)")


def _assert_nlev(opt, append, leaf):
    """The `-nlev N` of `append`, checked against the grid's lev and the Nudge_Vwin_* window."""
    if append is not None and not isinstance(append, str):
        raise SystemExit(f"{leaf}={append!r} must be a string of CAM configure options")
    found = _NLEV.findall(append or "")
    if not found:
        raise SystemExit(f"{leaf}={append!r} must carry '-nlev <N>': the L145 IC, the nudging "
                         f"window and CAM_CPPDEFS -DPLEV all depend on it")
    if len(found) > 1:
        raise SystemExit(f"{leaf}={append!r} gives -nlev {len(found)} times; give it once")
    nlev = int(found[0])
    res = opt["simulation"]["res"]
    try:
        lev = ic_solver.expected_grid(res).get("lev")
    except KeyError:
        lev = None       # reported by ic_solver.assert_grid
    if lev and nlev != lev:
        raise SystemExit(f"{leaf} says -nlev {nlev} but config/grids.json gives lev={lev} for "
                         f"simulation.res={res}: CAM would build {nlev} levels while "
                         f"ic_solver.assert_grid demands an L{lev} IC and the nudging files are L{lev}")
    n = opt["nudging"]
    lo, hi = int(n["vwin_lindex"]), int(n["vwin_hindex"])
    if not 1 <= lo <= hi <= nlev + 1:
        raise SystemExit(f"nudging.vwin_lindex={lo}, vwin_hindex={hi} is not a window on {nlev} "
                         f"levels (1 <= lindex <= hindex <= nlev + 1 = {nlev + 1}); {leaf} says "
                         f"-nlev {nlev}")
    return nlev


def _sourcemods(kind, policy):
    try:
        return {"policy": policy, **case_solver.sourcemods(kind, policy)}
    except ValueError as e:
        raise SystemExit(f"[waxrun] {e}")


def _sourcemods_changed(payload, specs):
    """Warn once when a replayed SourceMods set's files differ from the md5s the payload saved."""
    changed = set()
    for spec in specs:
        old = ((payload.get(spec["kind"]) or {}).get("sourcemods") or {}).get("md5")
        new = spec["sourcemods"].get("md5")
        if isinstance(old, dict) and new is not None:
            changed |= {r for r in set(old) | set(new) if old.get(r) != new.get(r)}
    changed = sorted(changed)
    if changed:
        _warn(f"[waxrun] WARNING: your SourceMods changed since {_OPTIONS_FILE} was saved "
              f"({', '.join(changed)}): this replay uses the files as they are now")


def _sourcemods_bytes(kind, policy):
    """{<comp dir>/<name>: file bytes} of a SourceMods set."""
    return {rel: Path(p).read_bytes() for rel, p in _sourcemods(kind, policy)["files"].items()}


def share_exe_plan(opt):
    """(True, None) when the spinup's cesm.exe can serve the coupled case too, else (False, reason)."""
    if not _as_bool((opt.get("build") or {}).get("share_exe", True), "build.share_exe"):
        return False, None
    sp, cp = opt["spinup"], opt["coupled"]
    why = []
    if os.path.realpath(_srcroot(opt, "spinup")) != os.path.realpath(_srcroot(opt, "coupled")):
        why.append("spinup.srcroot is a different CAM tree")
    if _sourcemods_bytes("spinup", sp["sourcemods"]) != _sourcemods_bytes("coupled", _coupled_sourcemods(opt)):
        why.append(f"spinup.sourcemods={sp['sourcemods']} and coupled.sourcemods={_coupled_sourcemods(opt)} differ")
    if str(sp["cam_config_opts_append"]).split() != str(cp["cam_config_opts_append"]).split():
        why.append("spinup and coupled cam_config_opts_append differ")
    if str(cp["build"]) != "build":
        why.append(f"coupled.build={cp['build']} already names the binary")
    if why:
        _note_once(f"[waxrun] NOTE: {'; '.join(why)}: building the spinup and coupled cases separately "
                   f"(build.share_exe).")
        return False, "; ".join(why)
    return True, None


def _stack_modules(opt, hpc):
    """(module_system_file, runtime_modules) of stacks.<hpc>, refused when they disagree."""
    sim, stacks = opt["simulation"], opt["stacks"][hpc]
    msf = stacks["module_system_file"]
    msf = msf if os.path.isabs(msf) else str(HERE / msf)
    settings = {"COMPILER": sim["compiler"], "MPILIB": sim["mpilib"]}
    block_mods = case_solver.runtime_modules(msf, settings)
    runtime_modules = list(stacks["runtime_modules"])
    # in order: Lmod loads them in sequence
    if list(block_mods) != runtime_modules:
        raise SystemExit(f"stacks.{hpc}.runtime_modules {runtime_modules} differ from the spliced block "
                         f"{msf} ({list(block_mods)}), in content or order")
    return msf, runtime_modules


def standalone_spinup_options(options, engage=None):
    """The spinup spec; reads only, writes nothing."""
    opt = options
    lay = _layout_of(opt)
    w = _window(opt)
    sim, sp, hpc = opt["simulation"], opt["spinup"], opt["simulation"].get("hpc_system", "derecho")
    _check_project(opt)
    job, account = sim["job_name"], sim["project"]
    case_name = f"{job}_spinup"
    caseroot = lay["caseroot_spinup"]
    cor = lay["cime_output_root"]
    srcroot = _srcroot(opt, "spinup")
    single_tree = os.path.realpath(srcroot) == os.path.realpath(_srcroot(opt, "coupled"))
    sourcemods = _sourcemods("spinup", sp["sourcemods"])
    _assert_nlev(opt, sp["cam_config_opts_append"], "spinup.cam_config_opts_append")
    tree = _tree(opt, srcroot, "spinup", sourcemods["policy"], single_tree)
    dtime = _dtime(sp["dtime"], "spinup.dtime")
    if dtime > 600:
        raise SystemExit(f"spinup.dtime {dtime} s is above the 600 s limit of a standalone WACCM-X spinup")
    ic = _ic(opt, lay, w)
    stop_option, stop_n = ic_solver.stop_n(ic["stop_seconds"])
    queue = opt["job"][hpc]["queue"]
    err = queue_error(queue, srcroot, sim.get("cime_machine") or hpc, hpc, nodes=sp["pe_nodes"],
                      what=f"the spinup case (spinup.pe_nodes={sp['pe_nodes']})")
    if err:
        raise SystemExit(f"[waxrun] {err}")
    qmax = queue_max_walltime(hpc)
    if queue not in qmax:
        raise SystemExit(f"job.{hpc}.queue={queue!r} (engage pbs.queue) has no walltime_max in engage "
                         f"config/machines.yaml queues; valid queues: {', '.join(sorted(qmax))}")
    wall_s = float(sp["wall_s_per_simday"]) * (REFERENCE_DTIME / dtime)
    if dtime != REFERENCE_DTIME:
        _warn(f"[waxrun] WARNING: spinup.wall_s_per_simday was measured at dtime={REFERENCE_DTIME} s; "
              f"scaled to {wall_s:.0f} s/simday for dtime={dtime} s (re-measure it before a production "
              f"spinup)")
    walltime_est, k = spinup_walltime(ic["stop_seconds"], wall_s, sp["wall_margin"], qmax[queue])
    chunks = spinup_chunks(stop_n, k)
    if k > 1:
        # per-chunk rounding can exceed the queue cap: add chunks until one fits
        unit = 3600 if stop_option == "nhours" else 1
        while True:
            walltime_est, _ = spinup_walltime(chunks[0]["stop_n"] * unit, wall_s, sp["wall_margin"],
                                              qmax[queue])
            if _secs(walltime_est) <= _secs(qmax[queue]) or chunks[0]["stop_n"] <= 1:
                break
            k += 1
            chunks = spinup_chunks(stop_n, k)
        if _secs(walltime_est) > _secs(qmax[queue]):
            raise SystemExit(f"spinup: even a {chunks[0]['stop_n']} {stop_option[1:]} chunk needs "
                             f"{walltime_est} (> the {queue} queue's {qmax[queue]}); pick another "
                             f"job.{hpc}.queue (engage pbs.queue)")
    nudging, nl = _nudging(opt, lay, ic["run_startdate"], w["stop"], dtime)
    _assert_nudge_vertical(nudging, ic["pick"], "spinup IC")
    nl_vars = {
        "ncdata": ic["ncdata"],
        "nudge_t": bool(nudging["enabled"]),
        "inithist": "ENDOFRUN",
        "solar_irrad_data_file": templates.REFERENCE_SOLAR["spinup"]["solar_irrad_data_file"],
        "solar_parms_data_file": templates.REFERENCE_SOLAR["spinup"]["solar_parms_data_file"],
        **_nudge_vars(nl),
    }
    tod = _tod(w["coupling_start"])
    product = os.path.join(cor, case_name, "run", f"{case_name}.cam.i.{w['coupling_start']:%Y-%m-%d}-{tod:05d}.nc")
    spec = {
        "kind": "spinup", "api_version": API_VERSION, "job_name": job, "case_name": case_name,
        **_identity(opt, srcroot, caseroot, cor, account),
        "rundir": os.path.join(cor, case_name, "run"),
        "exeroot": os.path.join(cor, case_name, "bld"), "dryrun_root": lay["dryrun_dir"], "run_dir": lay["run_dir"],
        "tree": tree, "single_tree": single_tree, "srcroot_describe": (tree or {}).get("describe"),
        "sourcemods": sourcemods, "ic": ic,
        "coupling_start": w["coupling_start"].isoformat(), "stop_time": w["stop"].isoformat(),
        "run_startdate": ic["run_startdate"], "start_tod": ic["start_tod"], "stop_option": stop_option,
        "stop_n": stop_n, "stop_seconds": ic["stop_seconds"], "resubmit": 0, "chunks": chunks,
        "dtime": dtime, "model_times_per_day": 86400 // dtime,
        "pe": _pe_spec(sp["pe_nodes"], "spinup.pe_nodes"),
        "config_xml": {"cam_config_opts_append": sp["cam_config_opts_append"],
                       "pio_netcdf_format": sp["pio_netcdf_format"], "ncpl_base_period": "hour",
                       "atm_ncpl": case_solver.ncpl_from_dtime(dtime),
                       "rof_ncpl": case_solver.ncpl_from_dtime(dtime), "pe_nodes": int(sp["pe_nodes"])},
        "user_nl_template": "spinup", "nl_vars": nl_vars, "nl_patches": _cam_overrides(opt, srcroot),
        "nudging": nudging, "forcing": None,
        "product_path": product, "walltime_est": walltime_est, "job_wallclock": walltime_est,
        "queue": queue, "job_priority": opt["job"][hpc].get("job_priority") or None,
        "charge_account": account,
        "build": _build_policy(sp, "spinup"), "runtime_modules": None, "hpc_system": hpc,
    }
    spec["create_newcase"] = case_solver.create_newcase_cmd(spec)
    if share_exe_plan(opt)[0]:
        # this case builds the coupled binary, so it gets the coupled case's splice and checks
        msf, mods = _stack_modules(opt, hpc)
        spec.update({"share_exe": True, "module_system_file": msf, "runtime_modules": mods,
                     "cesm_exe": os.path.join(spec["exeroot"], "cesm.exe")})
        spec["build"] = {"mode": "build", "clean_all": True, "ldd_must_be_zero": True, "exeroot": None,
                         "shared_with": "coupled"}
    else:
        spec["runtime_modules"] = _tree_module_script(srcroot, spec["machine"])
    spec["xmlchanges"] = case_solver.spinup_xmlchanges(spec)
    spec["forcing"] = _forcing_precheck(render_user_nl_cam(spec),
                                        (ic["run_startdate"], w["stop"] + timedelta(days=1)),
                                        templates.TEMPLATES["spinup"],
                                        provenance=_nl_provenance(spec["nl_patches"]))
    return spec


def harmonize_forcing(spinup_spec):
    """The solar forcing files of the spinup render, reused by the coupled case."""
    src = forcing_solver.namelist_files(render_user_nl_cam(spinup_spec), forcing_solver.SOLAR_KEYS)
    return {k: v for k, v in src.items() if v}


def coupled_options(options, engage, spinup_spec):
    """The coupled spec; reads only, writes nothing."""
    opt = options
    lay = _layout_of(opt)
    w = _window(opt)
    sim, cp, hpc = opt["simulation"], opt["coupled"], opt["simulation"].get("hpc_system", "derecho")
    _check_project(opt)
    err = queue_error(opt["job"][hpc]["queue"], None, sim.get("cime_machine") or hpc, hpc,
                      nodes=cp["pe_nodes"], what=f"the coupled GRX job (coupled.pe_nodes={cp['pe_nodes']} "
                                                 f"WACCM-X nodes alone)")
    if err:
        raise SystemExit(f"[waxrun] {err}")
    job, account = sim["job_name"], sim["project"]
    caseroot = lay["caseroot_coupled"]
    cor = lay["cime_output_root"]
    srcroot = _srcroot(opt, "coupled")
    single_tree = os.path.realpath(srcroot) == os.path.realpath(spinup_spec["srcroot"])
    sourcemods = _sourcemods("coupled", _coupled_sourcemods(opt))
    _assert_nlev(opt, cp["cam_config_opts_append"], "coupled.cam_config_opts_append")
    tree = _tree(opt, srcroot, "coupled", sourcemods["policy"], False)
    dtc = opt["mage"].get("dtcouple")
    if dtc in (None, ""):
        raise SystemExit("mage.dtcouple (makeitso voltron.coupling.dtCouple) is required: it is the "
                         "coupled physics timestep -> ATM_NCPL = 3600 / dtCouple")
    if float(dtc) != int(float(dtc)):
        raise SystemExit(f"voltron.coupling.dtCouple={dtc} s must be a whole number of seconds")
    dtime = _dtime(int(float(dtc)), "mage.dtcouple", coupled=True)
    segment_solver.warn_dtcouple_once(dtime, emit=_warn)
    run_startdate, start_tod = segment_solver.coupled_start(w["t0"], w["delay_h"])
    msf, runtime_modules = _stack_modules(opt, hpc)
    nudging, nl = _nudging(opt, lay, run_startdate, w["stop"], dtime)
    _assert_nudge_vertical(nudging, spinup_spec["ic"]["pick"], "spinup IC (whose coefficients the product carries)")
    nudging["gen_spans"], nudging["gen_job"] = list(spinup_spec["nudging"]["gen_spans"]), None
    nudging["coverage"] = {}
    nl_vars = {
        "ncdata": spinup_spec["product_path"],
        "inithist": _inithist(cp["inithist"], srcroot),
        "solar_wind_data_file": "",
        **harmonize_forcing(spinup_spec),
        **_nudge_vars(nl),
    }
    nl_patches = _cam_overrides(opt, srcroot)
    if spinup_spec.get("share_exe"):
        # the spinup case built cesm.exe; verify it instead of building again
        build = {"mode": "verify_only", "clean_all": False, "ldd_must_be_zero": True,
                 "cesm_exe": spinup_spec["cesm_exe"], "source_caseroot": spinup_spec["caseroot"],
                 "no_newer_sources": True, "shared_from": "spinup"}
    else:
        build = _build_policy(cp, "coupled")
    exeroot = os.path.join(cor, job, "bld")
    cesm_exe = build.get("cesm_exe") or os.path.join(exeroot, "cesm.exe")
    d = opt["model"]["data"]
    voltron_exe = d.get("voltron_exe")
    spec = {
        "kind": "coupled", "api_version": API_VERSION, "job_name": job, "case_name": job,
        **_identity(opt, srcroot, caseroot, cor, account),
        "rundir": lay["rundir"], "exeroot": exeroot,
        "dryrun_root": lay["dryrun_dir"], "run_dir": lay["run_dir"], "cesm_exe": cesm_exe,
        "tree": tree, "single_tree": single_tree, "srcroot_describe": (tree or {}).get("describe"),
        "spinup_srcroot": spinup_spec["srcroot"], "sourcemods": sourcemods,
        "ncdata": spinup_spec["product_path"],
        "coupling_start": w["coupling_start"].isoformat(), "stop_time": w["stop"].isoformat(),
        "t0": w["t0"].isoformat(), "delaycouplinghrs": w["delay_h"],
        "run_startdate": run_startdate, "start_tod": start_tod, "stephours": w["step_h"],
        "segments": {"ic": w["ic"], "iu": w["iu"], "maxhours": w["maxhours"], "n_coupled": w["iu"] - w["ic"] - 1,
                     "continue_flip": w["ic"] + 2},
        "dtime": dtime, "model_times_per_day": 86400 // dtime, "dtcouple": dtc,
        "pe": _pe_spec(cp["pe_nodes"], "coupled.pe_nodes"),
        "config_xml": {"cam_config_opts_append": cp["cam_config_opts_append"], "pio_netcdf_format": None,
                       "ncpl_base_period": "hour", "atm_ncpl": case_solver.ncpl_from_dtime(dtime),
                       "rof_ncpl": case_solver.ncpl_from_dtime(dtime), "pe_nodes": int(cp["pe_nodes"])},
        "user_nl_template": "coupled", "nl_vars": nl_vars, "nl_patches": nl_patches,
        "forcing": None, "module_system_file": msf, "runtime_modules": runtime_modules,
        "build": build, "nudging": nudging, "voltron_exe": voltron_exe, "kaijuhome": d.get("kaijuhome"),
        "conda_env": d.get("conda_env"), "queue": opt["job"][hpc]["queue"],
        "walltime": opt["job"][hpc]["walltime"], "job_priority": opt["job"][hpc].get("job_priority"),
        "charge_account": account, "hpc_system": hpc,
        "spinup_product": spinup_spec["product_path"], "spinup_caseroot": spinup_spec["caseroot"],
    }
    spec["create_newcase"] = case_solver.create_newcase_cmd(spec)
    spec["xmlchanges"] = case_solver.coupled_xmlchanges(spec)
    spec["forcing"] = _forcing_precheck(render_user_nl_cam(spec),
                                        (run_startdate, w["stop"] + timedelta(days=1)),
                                        templates.TEMPLATES["coupled"],
                                        provenance=_nl_provenance(nl_patches))
    return spec


def render_user_nl_cam(spec):
    """The case's user_nl_cam text: the template rendered from spec, with the cam.* overrides."""
    kind = spec["user_nl_template"]
    nl_vars = dict(spec["nl_vars"], **{k: spec[k] for k in templates.SPEC_VARIABLES.get(kind, ())})
    text = templates.render(kind, nl_vars)
    text = case_solver.patch_user_nl(text, spec["nl_patches"])
    _assert_nl_lengths(spec, text)
    return text


_STR_LEN_CACHE = {}
_NL_STRING = re.compile(r"'((?:[^']|'')*)'|\"((?:[^\"]|\"\")*)\"")
DEFAULT_STR_LEN = 256


def namelist_str_lens(srcroot):
    """{variable: N} for every char*N entry of the tree's namelist_definition.xml."""
    path = os.path.join(str(srcroot or ""), "bld", "namelist_files", "namelist_definition.xml")
    if path not in _STR_LEN_CACHE:
        try:
            text = Path(path).read_text()
        except OSError:
            text = ""
        _STR_LEN_CACHE[path] = {m.group(1).lower(): int(m.group(2)) for m in
                                re.finditer(r'<entry\s+id="([^"]+)"\s+type="char\*(\d+)(?:\([^)"]*\))?"', text)}
    return _STR_LEN_CACHE[path]


def nl_strings(text, str_lens, default=DEFAULT_STR_LEN):
    """[(variable, limit, element)] for every quoted string in the user_nl_cam `text`."""
    lines = text.split("\n")
    out = []
    for name, a, b in case_solver._statements(lines):
        rhs = "\n".join(case_solver._split_comment(l)[0] for l in lines[a:b + 1]).split("=", 1)[1]
        limit = str_lens.get(name.lower(), default)
        for m in _NL_STRING.finditer(rhs):
            out.append((name, limit, m.group(1) if m.group(1) is not None else m.group(2)))
    return out


def nl_length_errors(text, str_lens, default=DEFAULT_STR_LEN):
    """[(variable, length, limit, element)] for every quoted string longer than its char*N limit."""
    return [(n, len(el), lim, el) for n, lim, el in nl_strings(text, str_lens, default) if len(el) > lim]


def _assert_nl_lengths(spec, text):
    """Refuse an over-long namelist string before create_newcase rather than at preview_namelists."""
    bad = nl_length_errors(text, namelist_str_lens(spec.get("srcroot")))
    if bad:
        raise SystemExit(f"[waxrun] {spec.get('kind', '')} user_nl_cam: "
                         + "; ".join(f"{n} is {ln} characters, longer than its namelist limit of {lim}: {el}"
                                     for n, ln, lim, el in bad)
                         + "\nShorten the run directory (--run-dir / pbs.run_directory) or job_name.")


def _run(cmd, cwd=None):
    print("+ " + (cmd if isinstance(cmd, str) else " ".join(map(str, cmd))))
    return subprocess.run(cmd, check=True, cwd=cwd)


class _Recorder:
    """Accumulates commands.sh; with dry_run it records without executing."""

    def __init__(self, dry_run):
        self.dry_run = dry_run
        self.lines = ["#!/usr/bin/env bash", "# generated by waxrun genCASE" + (" --dry-run" if dry_run else ""),
                      "set -euo pipefail", ""]

    def note(self, text):
        self.lines.append(f"# {text}")

    def run(self, argv, cwd=None):
        line = case_solver.shell_line(argv)
        self.lines.append(f"cd {shlex.quote(str(cwd))} && {line}" if cwd else line)
        if not self.dry_run:
            _run(argv, cwd=cwd)

    def text(self):
        return "\n".join(self.lines) + "\n"


def genCASE(spec, clobber=False, dry_run=False, build=True):
    """Create, configure and build one case from its spec; dry_run only writes files under dryrun/."""
    kind = spec["kind"]
    coupled = kind == "coupled"
    spliced = coupled or bool(spec.get("share_exe"))
    caseroot = spec["caseroot"]
    out = os.path.join(spec["dryrun_root"], spec["case_name"]) if dry_run else caseroot
    rec = _Recorder(dry_run)
    if dry_run:
        os.makedirs(out, exist_ok=True)
        for stale in ("xmlchanges.txt", "create_newcase.sh", "sourcemods.txt",
                      "env_mach_specific.xml", "env_mach_specific.xml.bak_tree", case_solver.SPLICE_GUARD):
            p = os.path.join(out, stale)
            if os.path.exists(p):
                os.remove(p)
    else:
        if os.path.exists(caseroot):
            if os.path.exists(os.path.join(caseroot, case_solver.SPLICE_GUARD)) and not clobber:
                raise FileExistsError(f"{caseroot} is a spliced case ({case_solver.SPLICE_GUARD}); "
                                      "never case.setup --reset/--clean it; pass --clobber to recreate")
            if not clobber:
                raise FileExistsError(f"caseroot {caseroot} exists (pass --clobber to recreate)")
            shutil.rmtree(caseroot)
        os.makedirs(os.path.dirname(caseroot), exist_ok=True)
        os.makedirs(spec["cime_output_root"], exist_ok=True)
    xml_log = os.path.join(out, "xmlchanges.txt")
    result = {"kind": kind, "caseroot": caseroot, "dry_run_dir": out if dry_run else None,
              "user_nl_cam": os.path.join(out, "user_nl_cam"), "commands": os.path.join(out, "commands.sh"),
              "xmlchanges": xml_log, "spec": os.path.join(out, "spec.json"),
              "create_newcase": os.path.join(out, "create_newcase.sh"),
              "sourcemods_file": os.path.join(out, "sourcemods.txt"), "tree": spec.get("tree")}

    argv = spec.get("create_newcase") or case_solver.create_newcase_cmd(spec)
    rec.note(f"{kind} case {spec['case_name']}: create_newcase on {spec['srcroot']} "
             f"({spec['srcroot_describe']}); never --pecount — config_pes.xml gives the reference "
             f"{spec['pe']['nodes']} x {spec['pe']['tasks_per_node']} layout")
    rec.run(argv)
    Path(out, "create_newcase.sh").write_text(case_solver.shell_line(argv) + "\n")
    if not dry_run:
        # the pe_nodes override is applied in pre_setup, so the tree default is expected here
        result["pe"] = case_solver.pe_layout(caseroot, expect=case_solver.TREE_NTASKS_ATM)

    def _apply(phase):
        rec.lines += [f"cd {shlex.quote(caseroot)} && " + l for l in
                      case_solver.run_xmlchanges(caseroot, spec["xmlchanges"][phase], dry_run=True, log=xml_log)]
        if dry_run:
            return
        try:
            case_solver.run_xmlchanges(caseroot, spec["xmlchanges"][phase])
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"{kind} case {spec['case_name']}: the {phase} xmlchange {e.cmd} failed (rc={e.returncode}); "
                f"the keys before it are already written. Fix the key and recreate the case with "
                f"--clobber.") from e

    rec.note("CAM_CONFIG_OPTS is appended first: CAM's buildnml re-runs buildcpp only when "
             "Buildconf/camconf/config_cache.xml is missing or env_build.xml moved, so a preview "
             "before this append leaves CAM_CPPDEFS at -DPLEV=130 -DPCOLS=16 silently")
    _apply("pre_setup")
    if not dry_run:
        rec.note(f"case_solver.pe_layout(caseroot, expect={spec['pe']['ntasks_atm']}): the pe_nodes "
                 "override has been applied now, so the layout the PBS deck assumes is verified")
        result["pe"] = case_solver.pe_layout(caseroot, expect=spec["pe"])
    rec.run(["./case.setup"], cwd=caseroot)
    rec.note("case_solver.assert_calendar(caseroot): CALENDAR == GREGORIAN (the IC re-date, STOP_N and "
             "the nudging slots use Gregorian date arithmetic)")
    if not dry_run:
        result["calendar"] = case_solver.assert_calendar(caseroot)
    _apply("post_setup")

    if spliced:
        if not coupled:
            rec.note("build.share_exe: this spinup builds the cesm.exe the coupled case runs, so it gets "
                     "the coupled case's 23.09 <module_system> splice")
        rec.note("env_mach_specific.xml.bak_tree = the pre-splice copy, a true tree snapshot: "
                 "case.setup generated this file from the tree's ccs_config (25.10)")
        rec.run(["cp", "env_mach_specific.xml", "env_mach_specific.xml.bak_tree"], cwd=caseroot)
        rec.note(f"case_solver.splice_module_system(caseroot, {spec['module_system_file']}): <module_system> only; "
                 f"writes {case_solver.SPLICE_GUARD} (genCASE refuses case.setup --reset/--clean afterwards)")
        if dry_run:
            Path(out, "env_mach_specific.xml").write_text(_pristine_env_mach_specific(spec))
            case_solver.splice_module_system(out, spec["module_system_file"],
                                             path=os.path.join(out, "env_mach_specific.xml"))
        else:
            case_solver.splice_module_system(caseroot, spec["module_system_file"])
        result["runtime_modules"] = spec["runtime_modules"]

    policy = spec["sourcemods"]["policy"]
    rec.note(f"case_solver.copy_sourcemods({kind}, {policy!r}): "
             + (", ".join(f"SourceMods/{r} <- {v}" for r, v in sorted(spec['sourcemods']['files'].items()))
                or "nothing (the pinned tree's own sources)"))
    try:
        copied = case_solver.copy_sourcemods(caseroot if not dry_run else out, kind, policy, dry_run=dry_run,
                                             expect=spec["sourcemods"].get("md5"))
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        raise SystemExit(f"[waxrun] {e}")
    Path(out, "sourcemods.txt").write_text("".join(f"{r}\n" for r in copied))
    result["sourcemods"] = copied

    text = render_user_nl_cam(spec)
    Path(out, "user_nl_cam").write_text(text)
    rec.note(f"user_nl_cam := templates.render({spec['user_nl_template']!r}, nl_vars) from "
             f"{templates.TEMPLATES[spec['user_nl_template']]}"
             + (f" + case_solver.patch_user_nl({sorted(spec['nl_patches'])})" if spec["nl_patches"] else "")
             + (", solar forcing taken from the spinup render" if coupled else ""))
    rec.note(f"rm -f {case_solver.CAM_CONFIG_CACHE} (defensive: a cache written before the "
             f"CAM_CONFIG_OPTS append survives case.setup --reset and case.build --clean-all)")
    if not dry_run:
        case_solver.clear_config_cache(caseroot)
    rec.run(["./preview_namelists"], cwd=caseroot)
    rec.note("case_solver.assert_cam_config(caseroot): CAM_CONFIG_OPTS carries "
             f"{spec['config_xml']['cam_config_opts_append']!r} exactly once and CAM_CPPDEFS the "
             "matching -DPLEV/-DPCOLS (FAIL otherwise); "
             f"case_solver.dtime_from_case(caseroot) == {spec['dtime']} s")
    if not dry_run:
        result["cam_config"] = case_solver.assert_cam_config(
            caseroot, spec["config_xml"]["cam_config_opts_append"])
        if case_solver.has_fail(result["cam_config"]):
            raise RuntimeError(
                f"{caseroot} would build the wrong CAM: "
                + "; ".join(f["detail"].splitlines()[0] for f in result["cam_config"] if f["level"] == "FAIL"))
        result["dtime"] = case_solver.dtime_from_case(caseroot, expect=spec["dtime"])

    window = (spec["run_startdate"], _dt(spec["stop_time"]) + timedelta(days=1))
    rec.note(f"forcing_solver.coverage(caseroot, {window[0]} .. {window[1]:%Y-%m-%d %H:%M}): refuses on a gap")
    if dry_run:
        spec["forcing"] = _forcing_precheck(Path(out, "user_nl_cam").read_text(), window,
                                            os.path.join(out, "user_nl_cam"))
    else:
        spec["forcing"] = forcing_solver.coverage(caseroot, window, strict=True)
        _validate_atm_in(spec)
    result["forcing"] = spec["forcing"]

    if not coupled:
        ic = spec["ic"]
        if ic["offenders"]:
            rec.note(f"ic_solver.strip_offenders: {ic['offenders']} on ncol are physgrid offenders "
                     "(READ_DYN_FIELD_3D: Could not find BR field)")
            rec.lines.append(case_solver.shell_line(ic["strip_cmd"]))
            if not dry_run:
                pick = ic_solver.validate_ic(ic["pick"], spec["res"], spec.get("phys_ok") or ic_solver.PHYS_OK,
                                             ic["shared_dir"], ic.get("reference"),
                                             offenders=ic["offenders"])
                if pick["path"] != ic["ncdata"]:
                    raise RuntimeError(f"stripped IC landed at {pick['path']}, spec expects {ic['ncdata']}")
        else:
            rec.note(f"IC {ic['pick']} has no physgrid offenders: used as-is ({ic['nvars']} vars)")
        n = spec["nudging"]
        if n["enabled"]:
            cov = {datetime.fromisoformat(k): v for k, v in n["coverage"].items()}
            if n["gen_spans"]:
                rec.note(f"nudging: {n['missing']} slot(s) missing -> generator job(s) {n.get('gen_job')} "
                         f"write {n['gen_outdir']}; links to those files dangle until the job runs")
                cov = nudge_solver.plan_generated(cov, n["gen_outdir"])
            rec.note(f"nudge_solver.assemble({n['assembled_dir']}): one symlink dir of {len(cov)} slots == Nudge_Path")
            if not dry_run:
                nudge_solver.assemble(n["assembled_dir"], cov, run_root=spec["run_dir"])
                slots = nudge_solver.slots(n["beg"], n["end"])
                gaps = nudge_solver.missing(nudge_solver.coverage([n["assembled_dir"]], slots))
                if gaps and not n["gen_spans"]:
                    raise RuntimeError(f"{len(gaps)} nudging slot(s) unresolved in {n['assembled_dir']}")
                if gaps:
                    _warn(f"[waxrun] nudging: {len(gaps)} slot(s) wait for the generator job (checked again "
                          "before the spinup starts)")
            rec.note(f"PRERUN_SCRIPT={os.path.join(caseroot, 'waxrun_prerun.sh')} re-checks every slot of "
                     f"{n['assembled_dir']} inside case.run (after the generator prereq completed); the submit "
                     "wrapper does the same before ./case.submit")

    if build:
        if spliced and spec["build"]["mode"] == "build":
            rec.note("case_solver.f7_preflight(caseroot, srcroot): module dry load (FAIL), set(SFC ifort) (FAIL)")
            if not dry_run:
                findings = case_solver.f7_preflight(caseroot, spec["srcroot"])
                result["f7"] = findings
                if case_solver.has_fail(findings):
                    raise RuntimeError("module preflight failed: "
                                       + "; ".join(f["detail"].splitlines()[0] for f in findings if f["level"] == "FAIL"))
        rec.lines += case_solver.build(spec, dry_run=True)
        if not dry_run:
            case_solver.build(spec)
            if spliced:
                foreign = case_solver.ldd_foreign_stack(spec["cesm_exe"], spec["runtime_modules"])
                if foreign:
                    raise RuntimeError(f"{spec['cesm_exe']} resolves libraries from another stack: {foreign}")
                # voltron runs beside cesm.exe on the same stack
                volt = spec.get("voltron_exe")
                if volt and os.path.exists(volt):
                    foreign = case_solver.ldd_foreign_stack(volt, spec["runtime_modules"])
                    if foreign:
                        raise RuntimeError(f"{volt} resolves libraries from another stack: {foreign}")
    else:
        rec.note("build skipped (--skip-build / build=False)")

    Path(out, "commands.sh").write_text(rec.text())
    Path(out, "spec.json").write_text(json.dumps(spec, indent=2, default=str) + "\n")
    if not dry_run:
        Path(caseroot, "waxrun_commands.sh").write_text(rec.text())
        Path(caseroot, "waxrun_spec.json").write_text(json.dumps(spec, indent=2, default=str) + "\n")
    result["commands_text"] = rec.text()
    return result


def _pristine_env_mach_specific(spec):
    """A minimal env_mach_specific.xml with the tree's <module_system>, for --dry-run."""
    block = _tree_module_system_xml(spec["srcroot"], spec["machine"])
    indented = "\n".join(("  " + l) if l.strip() else l for l in block.split("\n"))
    return ('<?xml version="1.0"?>\n<file id="env_mach_specific.xml" version="2.0">\n'
            + indented + "\n</file>\n")


def _validate_atm_in(spec):
    """Check that the keys the case depends on reached CaseDocs/atm_in."""
    atm_in = Path(spec["caseroot"], "CaseDocs", "atm_in")
    if not atm_in.exists():
        raise RuntimeError(f"{atm_in} missing after preview_namelists")
    text = atm_in.read_text().lower()
    want = ["ncdata", "model_times_per_day", "nudge_model"]
    if spec["kind"] == "coupled":
        want.append("ionos_epotential_mage")
    missing = [k for k in want if f"{k} " not in text and f"{k}=" not in text and f"{k}\t" not in text]
    if missing:
        raise RuntimeError(f"{atm_in} lacks {missing}")


def _postrun_script(spec):
    """CIME POSTRUN_SCRIPT for a chunked spinup (CIME does not set CONTINUE_RUN with RESUBMIT=0)."""
    k = len(spec["chunks"])
    remainder = spec["chunks"][-1]["stop_n"]
    return "\n".join([
        "#!/usr/bin/env bash",
        f"# CIME POSTRUN_SCRIPT for {spec['case_name']}: continue run after each chunk, shorter last chunk.",
        "set -euo pipefail",
        f"cd \"${{1:-{spec['caseroot']}}}\"",
        "n=$(( $(cat .waxrun_chunk 2>/dev/null || echo 0) + 1 ))",
        "echo \"$n\" > .waxrun_chunk",
        "./xmlchange CONTINUE_RUN=TRUE",
        f"if [ \"$n\" -eq {k - 1} ]; then ./xmlchange STOP_N={remainder}; fi",
        "",
    ])


def _nudging_recheck_lines(n, assembled=None):
    """Shell lines that count missing slots in the assembled nudging dir (`assembled`: its spelling)."""
    a = n["assembled_dir"] if assembled is None else assembled
    return [
        f"# nudging check: every slot in {a} must exist",
        "missing=0",
        f"for f in {shlex.quote(a)}/*.nc; do [ -e \"$f\" ] || {{ echo \"MISSING $f\" >&2; missing=$((missing+1)); }}; done",
        f"[ \"$(ls {shlex.quote(a)}/*.nc 2>/dev/null | wc -l)\" -ge {n['slots']} ] || {{ echo \"fewer than {n['slots']} nudging slots in {a}\" >&2; exit 3; }}",
    ]


def _prerun_script(spec):
    """CIME PRERUN_SCRIPT: refuse to start while a nudging slot is missing (the DAG skips the wrapper)."""
    n = spec["nudging"]
    return "\n".join([
        "#!/usr/bin/env bash",
        f"# CIME PRERUN_SCRIPT for {spec['case_name']}: refuse to start while a nudging slot is missing.",
        "set -euo pipefail",
        *_nudging_recheck_lines(n),
        "[ \"$missing\" -eq 0 ] || { echo \"REFUSING to start: $missing nudging slot(s) missing\" >&2; exit 3; }",
        "",
    ])


def _submit_wrapper(spec, nn, k, root=None):
    """The text of one spinup submit wrapper; paths are relative to its own directory `root`."""
    n = spec["nudging"]
    root = root or spec["run_dir"]
    run_root = spec["run_dir"]
    lines = [
        "#!/usr/bin/env bash",
        f"# Submit spinup chunk {nn}/{k} of {spec['case_name']}: id=$(bash "
        f"{spec['job_name']}_spinup_submit_{nn:02d}.sh [prereq job id])",
        "set -euo pipefail",
        "prereq=\"${1:-}\"",
        "unset CDPATH",
        'cd "$(dirname "$0")"',
    ]
    if n["enabled"]:
        lines += _nudging_recheck_lines(n, nudge_solver.rel_or_abs(n["assembled_dir"], root, run_root)) + [
            "[ \"$missing\" -eq 0 ] || { echo \"REFUSING case.submit: $missing nudging slot(s) missing\" >&2; exit 3; }",
        ]
    if nn == 1:
        ncdata = nudge_solver.rel_or_abs(spec["ic"]["ncdata"], root, run_root)
        lines += [f"[ -e {shlex.quote(ncdata)} ] || {{ echo \"ncdata {ncdata} missing\" >&2; exit 4; }}"]
    lines += [f"cd {shlex.quote(nudge_solver.rel_or_abs(spec['caseroot'], root, run_root))}"]
    lines += [
        "if [ -n \"$prereq\" ]; then out=$(./case.submit --prereq \"$prereq\" 2>&1); else out=$(./case.submit 2>&1); fi",
        "echo \"$out\" >&2",
        "echo \"$out\" | grep -oP 'Submitted job case\\.run with id \\K\\d+'",
        "",
    ]
    return "\n".join(lines)


def genPBS(spec, kind, dry_run=False):
    """Write the spinup's nudging generator jobs or submit wrappers (by `kind`); returns their paths."""
    if spec["kind"] != "spinup":
        raise ValueError(f"genPBS({kind}) takes the spinup spec, got {spec['kind']}")
    root = spec["dryrun_root"] if dry_run else spec["run_dir"]
    if kind == "nudging":
        n = spec["nudging"]
        gen_dir = os.path.join(root, "nudging_gen")
        if os.path.isdir(gen_dir):
            # drop stale wrappers and markers; keep TMP/ (ESMF weight cache, maybe a live job's cwd)
            for stale in sorted(glob.glob(os.path.join(gen_dir, "nudge_gen*.pbs"))
                                + glob.glob(os.path.join(gen_dir, "nudge_gen*.pbs.OK"))):
                os.remove(stale)
            if not os.listdir(gen_dir):
                os.rmdir(gen_dir)
        if not (n["enabled"] and n["gen_spans"]):
            n["gen_job"] = None
            return []
        spans = [(datetime.fromisoformat(a).date(), datetime.fromisoformat(b).date()) for a, b in n["gen_spans"]]
        paths = nudge_solver.gen_job(spans, n["gen_outdir"], spec["project"],
                                     f"{spec['job_name']}_nudge", gen_dir, n["rda_input"],
                                     workers=n["gen_workers"],
                                     lindex=n["vwin_lindex"], hindex=n["vwin_hindex"],
                                     run_root=spec["run_dir"], job_priority=spec.get("job_priority"))
        n["gen_job"] = paths
        n["gen_dir"] = gen_dir
        n["submit"] = [nudge_solver.submit_command(p, base=spec["run_dir"]) for p in paths]
        return paths
    if kind == "spinup_submit":
        os.makedirs(root, exist_ok=True)
        k = len(spec["chunks"])
        paths = []
        for nn in range(1, k + 1):
            p = os.path.join(root, f"{spec['job_name']}_spinup_submit_{nn:02d}.sh")
            Path(p).write_text(_submit_wrapper(spec, nn, k, root))
            os.chmod(p, 0o755)
            paths.append(p)
        hooks = os.path.join(root, spec["case_name"]) if dry_run else spec["caseroot"]
        if k > 1:
            pr = os.path.join(hooks, "waxrun_postrun.sh")
            os.makedirs(os.path.dirname(pr), exist_ok=True)
            Path(pr).write_text(_postrun_script(spec))
            os.chmod(pr, 0o755)
            spec["postrun_script"] = pr
        if spec["nudging"]["enabled"]:
            pre = os.path.join(hooks, "waxrun_prerun.sh")
            os.makedirs(os.path.dirname(pre), exist_ok=True)
            Path(pre).write_text(_prerun_script(spec))
            os.chmod(pre, 0o755)
            spec["prerun_script"] = pre
        spec["submit_scripts"] = paths
        return paths
    raise ValueError(f"unknown genPBS kind {kind!r} (nudging | spinup_submit)")


def spec_diff(old, new, prefix=""):
    """Key paths whose values differ between two JSON-shaped specs (both directions)."""
    out = []
    if isinstance(old, dict) and isinstance(new, dict):
        for k in sorted(set(old) | set(new)):
            p = f"{prefix}.{k}" if prefix else str(k)
            if k not in old or k not in new:
                out.append(p)
            else:
                out += spec_diff(old[k], new[k], p)
        return out
    if json.dumps(old, sort_keys=True, default=str) != json.dumps(new, sort_keys=True, default=str):
        out.append(prefix or "<root>")
    return out


RUN_DIRECTORY_TOKEN = "<run_directory>"


def _subst_run_dir(obj, run_dir):
    """`obj` with the path `run_dir` replaced by <run_directory> in strings, list items and keys."""
    pat = re.compile(re.escape(run_dir) + r"(?=$|[/\s'\",;:)\]])")
    def sub(v):
        if isinstance(v, str):
            return pat.sub(RUN_DIRECTORY_TOKEN, v)
        if isinstance(v, list):
            return [sub(x) for x in v]
        if isinstance(v, dict):
            return {sub(k): sub(x) for k, x in v.items()}
        return v
    return sub(obj)


def _spec_view(spec):
    """A spec for comparison: tree finding details dropped, run-dir paths as <run_directory>."""
    view = dict(spec)
    tree = view.get("tree")
    if isinstance(tree, dict) and isinstance(tree.get("findings"), list):
        view["tree"] = dict(tree, findings=[{k: v for k, v in f.items() if k != "detail"}
                                            if isinstance(f, dict) else f for f in tree["findings"]])
    run_dir = str(spec.get("run_dir") or "").rstrip("/")
    if run_dir and os.path.isabs(run_dir):
        view = _subst_run_dir(view, run_dir)
    return view


def same_specs(payload, pure):
    """True when a payload's stored spinup and coupled specs equal the recomputed ones."""
    return not any(spec_diff(_spec_view(payload[kind]), _spec_view(pure[kind])) for kind in ("spinup", "coupled"))


def _replay_pin(payload):
    """The IC a -wo replay pins, or None when an edited start file date asks for a re-pick."""
    sic, oic = payload["spinup"]["ic"], payload["options"].get("ic") or {}
    pin = oic.get("pick") or sic["pick"]
    if oic.get("start_file_date") != sic.get("start_file_date") and pin == sic["pick"]:
        _warn(f"[waxrun] start file date changed ({sic.get('start_file_date') or 'automatic'} -> "
              f"{oic.get('start_file_date') or 'automatic'}): the IC is re-picked")
        return None
    return pin


# payload paths stored relative to the payload file; the specs keep absolute paths
PAYLOAD_PATHS = (("options", "simulation", "run_directory"), ("options", "ic", "stripped"),
                 ("options", "ic", "shared_dir"), ("options", "nudging", "assembled_dir"),
                 ("generated", "nudging", "gen_job"), ("generated", "nudging", "gen_dir"),
                 ("generated", "nudging_scripts"), ("generated", "spinup_submit_scripts"),
                 ("generated", "postrun_script"), ("generated", "prerun_script"))
CASE_RESULT_PATHS = ("caseroot", "dry_run_dir", "user_nl_cam", "commands", "xmlchanges", "spec",
                     "create_newcase", "sourcemods_file")


def _payload_path_slots(payload):
    """(container, key) of every path leaf of a payload; placeholder spellings are skipped."""
    def slot(container, key):
        v = container.get(key) if isinstance(container, dict) else None
        if isinstance(v, str) and v and "<" not in v:
            yield container, key
        elif isinstance(v, list):
            for i, item in enumerate(v):
                if isinstance(item, str) and item and "<" not in item:
                    yield v, i
    for path in PAYLOAD_PATHS:
        d = payload
        for k in path[:-1]:
            d = d.get(k) if isinstance(d, dict) else None
        yield from slot(d, path[-1])
    cases = (payload.get("generated") or {}).get("cases")
    if isinstance(cases, dict):
        for result in cases.values():
            for key in CASE_RESULT_PATHS:
                yield from slot(result, key)


def relocate_payload(payload, from_dir, to_dir=None, run_root=None):
    """A copy of `payload` with paths read relative to `from_dir`, re-spelled relative to `to_dir`."""
    # not deepcopy: two payload keys share one list in memory and would be relocated twice
    out = json.loads(json.dumps(payload, default=str))
    opt = out.get("options")
    from_dir = os.path.abspath(str(from_dir))

    def absolute(v):
        return os.path.normpath(os.path.join(from_dir, os.path.expanduser(v)))

    rd = ((opt or {}).get("simulation") or {}).get("run_directory") if isinstance(opt, dict) else None
    run_root = os.path.abspath(str(run_root)) if run_root else (absolute(str(rd)) if rd else None)
    to = os.path.abspath(str(to_dir)) if to_dir is not None else None
    related = bool(to and run_root and os.path.commonpath([run_root, to]) in (run_root, to))
    for container, key in list(_payload_path_slots(out)):
        a = absolute(container[key])
        if related and os.path.commonpath([run_root, a]) == run_root:
            rel = os.path.relpath(a, to)
            container[key] = rel if os.path.realpath(os.path.join(to, rel)) == os.path.realpath(a) else a
        else:
            container[key] = a
    return out


def options_shape_problem(data):
    """None, or the first section of a -wo file that is not a JSON object."""
    opts, prefix = data, ""
    if "options" in data:
        for key in ("options", "spinup", "coupled", "generated"):
            if key in data and not isinstance(data[key], dict):
                return f"{key} must be a JSON object"
        opts, prefix = data["options"], "options."

    def walk(desc, node, pre):
        for key, spec in desc.items():
            if key.startswith("_") or not isinstance(spec, dict) or "LEVEL" in spec or key not in node:
                continue
            if not isinstance(node[key], dict):
                return f"{pre}{key} must be a JSON object"
            found = walk(spec, node[key], f"{pre}{key}.")
            if found:
                return found
        return None
    desc = json.loads(Path(OPTIONS_DESCRIPTION).read_text())
    return walk({k: v for k, v in desc.items() if k != "version"}, opts, prefix)


def load_options_file(path):
    """A -wo file as ("payload", replay payload) or ("options", options dict)."""
    try:
        data = json.loads(Path(path).read_text())
    except FileNotFoundError:
        raise SystemExit(f"[waxrun] {path} does not exist") from None
    except (OSError, UnicodeDecodeError) as e:
        raise SystemExit(f"[waxrun] {path} cannot be read ({getattr(e, 'strerror', None) or e})") from None
    except json.JSONDecodeError as e:
        raise SystemExit(f"[waxrun] {path} is not valid JSON (line {e.lineno} column {e.colno}: {e.msg})") from None
    if not isinstance(data, dict):
        raise SystemExit(f"[waxrun] {path} must be a JSON object")
    problem = options_shape_problem(data)
    if problem:
        raise SystemExit(f"[waxrun] {path}: {problem}")
    if isinstance(data.get("options"), dict):
        specs = [k for k in ("spinup", "coupled") if k in data]
        if len(specs) == 2:
            return "payload", _load_payload(path)
        if specs:
            raise SystemExit(f"[waxrun] {path}: a saved waxrun_parameters.json has 'spinup' and 'coupled'; "
                             f"an options file has neither")
        return "options", data["options"]
    return "options", data


def _omitted_note(opt, loaded, desc):
    """NOTE the BASIC leaves an options file omits (they take their defaults); returns them."""
    filled = set(opt.get("_engage_filled", ())) | set(_leaf_paths(loaded or {}))

    def walk(desc_sec, opt_sec, prefix):
        out = []
        for key, spec in desc_sec.items():
            if key.startswith("_") or not isinstance(spec, dict):
                continue
            path = f"{prefix}.{key}" if prefix else key
            if "LEVEL" in spec:
                if (spec.get("prompt") and _at_or_below(spec["LEVEL"], "BASIC") and path not in filled
                        and asked(path, opt_sec)):
                    out.append(path)
            else:
                out += walk(spec, (opt_sec or {}).get(key) or {}, path)
        return out
    names = walk({k: v for k, v in desc.items() if k != "version"}, opt, "")
    if names:
        if len(names) <= 8:
            text = f"{', '.join(names)}: using the default{'' if len(names) == 1 else 's'}."
        else:
            text = (f"{len(names)} keys ({', '.join(names[:3])}, ...): using the defaults; the full set is "
                    f"saved in {PARAMETERS_JSON}.")
        _note_once(f"[waxrun] NOTE: the waxrun options JSON omits {text}")
    return names


def _load_payload(path):
    """A replay payload with its relative paths resolved against the file's directory."""
    payload = json.loads(Path(path).read_text())
    for key in ("options", "spinup", "coupled"):
        if key not in payload:
            raise SystemExit(f"{path} is not a waxrun_parameters.json payload (no {key!r})")
    if not payload["spinup"].get("ic", {}).get("pick"):
        raise SystemExit(f"{path}: spinup.ic.pick missing; nothing to pin")
    return relocate_payload(payload, os.path.dirname(os.path.abspath(str(path))))


def save_parameters(path, options, spinup_spec, coupled_spec, generated=None):
    """Write the replay payload (options, both specs, generated files); run-dir paths relative to it."""
    payload = {"api_version": API_VERSION,
               "options": {k: v for k, v in options.items() if k != "_engage"},
               "spinup": spinup_spec, "coupled": coupled_spec, "generated": generated or {}}
    payload = relocate_payload(payload, os.getcwd(), os.path.dirname(os.path.abspath(str(path))))
    Path(path).write_text(json.dumps(payload, indent=2, default=str) + "\n")
    return path


def use_tui(no_tui):
    """True unless --no-tui or there is no terminal."""
    return not no_tui and sys.stdin.isatty() and sys.stdout.isatty()


def _build_parser():
    p = argparse.ArgumentParser(description="waxrun — configure the WACCM-X spinup + coupled cases of a GRX run")
    p.add_argument("--mode", default="BASIC", choices=LEVELS)
    p.add_argument("--interactive", action="store_true", help="prompt for the options not supplied by engage")
    p.add_argument("--no-tui", action="store_true",
                   help="with --interactive: ask them as linear prompts instead of the wizard")
    p.add_argument("--tui", action="store_true",
                   help="the wizard (the default; kept for old command lines)")
    p.add_argument("--engage", help="engage options payload as a JSON string")
    p.add_argument("--engage-file", help="engage options payload as a JSON file")
    p.add_argument("-wo", "--waxrun_options",
                   help="path to waxrun JSON file of options")
    p.add_argument("--dry-run", action="store_true",
                   help="write the case specs under <run_dir>/dryrun/; run no CIME script")
    p.add_argument("--skip-build", action="store_true", help="generate the cases but run no case.build")
    p.add_argument("--phase", default="both", choices=PHASES)
    p.add_argument("--clobber", action="store_true",
                   help="recreate existing caseroots and overwrite an existing parameters file")
    p.add_argument("--options", "-o", help="replay file to write (default <run_dir>/waxrun_parameters.json)")
    p.add_argument("--job-name")
    p.add_argument("--start", help="coupled start (ISO)")
    p.add_argument("--stop", help="coupled stop (ISO)")
    p.add_argument("--run-dir")
    p.add_argument("--account")
    p.add_argument("--start-file-dir", help="the restart series (ic.start_file_dir)")
    p.add_argument("--srcroot", help=f"the CAM tree both cases are created from "
                                     f"(default: ${WAXHOME_ENV} when set)")
    p.add_argument("--spinup-srcroot",
                   help="a second CAM tree for the spinup case (default: --srcroot)")
    p.add_argument("--skip-tree-check", action="store_true",
                   help="skip the CAM tree check (both specs record tree = null)")
    return p


def _cli_overrides(args):
    ov = {"simulation": {}, "ic": {}, "spinup": {}}
    if args.job_name:
        ov["simulation"]["job_name"] = args.job_name
    if args.start:
        ov["simulation"]["coupling_start"] = _run_dt(args.start, "--start").strftime(segment_solver.RUN_DATETIME_FORMAT)
    if args.stop:
        ov["simulation"]["stop_time"] = _run_dt(args.stop, "--stop").strftime(segment_solver.RUN_DATETIME_FORMAT)
    if args.run_dir:
        ov["simulation"]["run_directory"] = os.path.abspath(args.run_dir)
    if args.account:
        ov["simulation"]["project"] = args.account
    if args.start_file_dir:
        ov["ic"]["start_file_dir"] = args.start_file_dir
    if args.srcroot:
        ov["simulation"]["srcroot"] = os.path.abspath(args.srcroot)
    if args.spinup_srcroot:
        ov["spinup"]["srcroot"] = os.path.abspath(args.spinup_srcroot)
    return {k: v for k, v in ov.items() if v}


def _plan(args):
    """Options and both specs, writing nothing: (engage, payload, opt, spinup_spec, coupled_spec, pure)."""
    global _OPTIONS_FILE
    _OPTIONS_FILE = args.waxrun_options
    engage = None
    if args.engage_file:
        engage = json.loads(Path(args.engage_file).read_text())
    elif args.engage:
        engage = json.loads(args.engage)
    kind, data = load_options_file(args.waxrun_options) if args.waxrun_options else (None, None)
    payload = data if kind == "payload" else None
    if kind == "options":
        # an options file: nothing is asked, omitted leaves take their defaults
        desc = load_descriptions()
        opt = getparams(args.mode, options_path=copy.deepcopy(data), overrides=_cli_overrides(args),
                        engage=engage, interactive=False)
        _omitted_note(opt, {k: v for k, v in data.items() if not k.startswith("_")}, desc)
    else:
        opt = getparams(args.mode, options_path=payload, overrides=_cli_overrides(args), engage=engage,
                        interactive=args.interactive and not payload, tui=use_tui(args.no_tui),
                        answers=None if payload else (engage or {}).get("prompt_answers"),
                        asked=None if payload else (engage or {}).get("asked_answers"))
    if payload and "sourcemods" not in ((payload.get("options") or {}).get("coupled") or {}) \
            and "coupled.sourcemods" not in opt.get("_engage_filled", ()):
        # saved before coupled.sourcemods existed: keep it absent so the file is rewritten as it was
        opt["coupled"].pop("sourcemods", None)
    if args.skip_tree_check:
        opt["_skip_tree_check"] = True
    if payload:
        opt["ic"]["pick"] = _replay_pin(payload)
    spinup_spec = standalone_spinup_options(opt, engage)
    coupled_spec = coupled_options(opt, engage, spinup_spec)
    if payload:
        _sourcemods_changed(payload, (spinup_spec, coupled_spec))
    err = skip_build_error(args, coupled_spec)
    if err:
        raise SystemExit(err)
    opt["ic"]["pick"], opt["ic"]["stripped"] = spinup_spec["ic"]["pick"], spinup_spec["ic"]["stripped"]
    pure = {"spinup": copy.deepcopy(spinup_spec), "coupled": copy.deepcopy(coupled_spec)}
    return engage, payload, opt, spinup_spec, coupled_spec, pure


def skip_build_error(args, coupled_spec):
    """None, or why --skip-build is refused: the coupled cesm.exe does not exist and nothing builds it."""
    if not getattr(args, "skip_build", False) or getattr(args, "dry_run", False) \
            or getattr(args, "phase", "both") == "spinup":
        return None
    exe = coupled_spec.get("cesm_exe")
    if exe and os.path.isfile(exe):
        return None
    return (f"--skip-build but the coupled cesm.exe does not exist: {exe} "
            "(build the coupled case, or drop --skip-build)")


def _phases_of(generated):
    phase = (generated or {}).get("phase")
    return {"spinup", "coupled"} if phase == "both" else ({phase} if phase in PHASES else set())


def _earlier_generated(path, args, payload, pure):
    """The `generated` block of the same run's existing payload at `path`, None to overwrite it.

    Any other existing file raises FileExistsError."""
    if not os.path.exists(path):
        return None
    try:
        prev = json.loads(Path(path).read_text())
        if isinstance(prev, dict):
            prev = relocate_payload(prev, os.path.dirname(os.path.abspath(str(path))))
    except (OSError, ValueError):
        prev = None
    if (not args.dry_run and isinstance(prev, dict) and all(k in prev for k in ("options", "spinup", "coupled"))
            and same_specs(prev, pure)):
        return prev.get("generated") or {}
    if args.clobber or (payload and os.path.samefile(path, args.waxrun_options)) \
            or (args.dry_run and not args.options):
        return None
    raise FileExistsError(f"{path} exists (pass --clobber to overwrite; -o/--options is where waxrun "
                          f"writes its replay payload, -wo replays one)")


def plan_check(argv=None):
    """Run waxrun(argv)'s checks without writing anything; raises what waxrun() would."""
    args = _build_parser().parse_args(argv)
    engage, payload, opt, spinup_spec, coupled_spec, pure = _plan(args)
    lay = _layout_of(opt)
    root = lay["dryrun_dir"] if args.dry_run else lay["run_dir"]
    _earlier_generated(args.options or os.path.join(root, PARAMETERS_JSON), args, payload, pure)
    return None


def waxrun(argv=None):
    """Generate both cases and their jobs and save the replay payload.

    Returns (coupled_spec or options, spinup_submit_scripts, coupled_case_files)."""
    args = _build_parser().parse_args(argv)
    ic_solver.require_ncks()
    engage, payload, opt, spinup_spec, coupled_spec, pure = _plan(args)
    lay = _layout_of(opt)
    phases = {"spinup", "coupled"} if args.phase == "both" else {args.phase}
    root = lay["dryrun_dir"] if args.dry_run else lay["run_dir"]
    path = args.options or os.path.join(root, PARAMETERS_JSON)
    earlier = _earlier_generated(path, args, payload, pure)
    build = not args.skip_build
    submit_scripts, nudging_scripts = [], []
    generated = {"dry_run": args.dry_run, "phase": args.phase, "build": build, "cases": {}}
    if "spinup" in phases:
        nudging_scripts = genPBS(spinup_spec, "nudging", dry_run=args.dry_run)
        generated["cases"]["spinup"] = genCASE(spinup_spec, clobber=args.clobber, dry_run=args.dry_run, build=build)
        submit_scripts = genPBS(spinup_spec, "spinup_submit", dry_run=args.dry_run)
        generated["nudging"] = {k: spinup_spec["nudging"].get(k) for k in ("gen_job", "gen_dir", "submit")}
        generated["postrun_script"] = spinup_spec.get("postrun_script")
        generated["prerun_script"] = spinup_spec.get("prerun_script")
    if "coupled" in phases:
        generated["cases"]["coupled"] = genCASE(coupled_spec, clobber=args.clobber, dry_run=args.dry_run, build=build)
    for r in generated["cases"].values():
        r.pop("commands_text", None)
    if earlier is not None:
        # the same run's other phase: keep what it generated
        generated["cases"] = {**(earlier.get("cases") or {}), **generated["cases"]}
        done = _phases_of(earlier) | phases
        generated["phase"] = "both" if done == {"spinup", "coupled"} else args.phase
        if "spinup" not in phases:
            for key in ("nudging", "postrun_script", "prerun_script"):
                if key in earlier:
                    generated[key] = earlier[key]
            nudging_scripts = list(earlier.get("nudging_scripts") or [])
            submit_scripts = list(earlier.get("spinup_submit_scripts") or [])
    coupled_case_files = [coupled_spec["caseroot"], os.path.join(coupled_spec["caseroot"], "user_nl_cam")]
    coupled_spec["nudging_scripts"] = nudging_scripts
    coupled_spec["nudging"]["gen_job"] = nudging_scripts[0] if nudging_scripts else None
    coupled_spec["nudging"]["submit"] = spinup_spec["nudging"].get("submit")
    coupled_spec["spinup_submit_scripts"] = submit_scripts
    coupled_spec["spinup_product"] = spinup_spec["product_path"]
    generated["nudging_scripts"] = nudging_scripts
    generated["spinup_submit_scripts"] = submit_scripts
    os.makedirs(root, exist_ok=True)
    save_parameters(path, opt, pure["spinup"], pure["coupled"], generated)
    coupled_spec["waxrun_parameters"] = path
    LAST.clear()
    LAST.update({"options": opt, "spinup": spinup_spec, "coupled": coupled_spec, "path": path, "generated": generated})
    print(f"waxrun: {'dry-run ' if args.dry_run else ''}parameters saved to {path}")
    if engage:
        return coupled_spec, submit_scripts, coupled_case_files
    return opt, submit_scripts, coupled_case_files


def main(argv=None):
    first, submit_scripts, coupled_case_files = waxrun(argv)
    spinup, coupled = LAST["spinup"], LAST["coupled"]
    print(f"waxrun: spinup {spinup['case_name']}: {spinup['run_startdate']} tod {spinup['start_tod']} "
          f"-> {spinup['coupling_start']} ({spinup['stop_option']} {spinup['stop_n']}, "
          f"{len(spinup['chunks'])} chunk(s), walltime {spinup['walltime_est']}); product {spinup['product_path']}")
    print(f"waxrun: coupled {coupled['case_name']}: {coupled['run_startdate']} tod {coupled['start_tod']} "
          f"STOP_N {coupled['stephours']} h x {coupled['segments']['n_coupled']} segments; cesm.exe {coupled['cesm_exe']}"
          + ("" if "coupled" in LAST["generated"]["cases"] else " (spec only: --phase spinup)"))
    if coupled["nudging_scripts"]:
        print(f"waxrun: nudging generator job(s): {coupled['nudging_scripts']}")
    print(f"waxrun: spinup submit wrappers: {submit_scripts}")
    print(f"waxrun: coupled case files: {coupled_case_files}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
