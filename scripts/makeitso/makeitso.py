#!/usr/bin/env python


"""makeitso for the MAGE magnetosphere software.

Interactively prepares a MAGE magnetosphere run (decks and job scripts). Can also be imported
and driven through makeitso(args).
"""


import argparse
import contextlib
import copy
import datetime
import filecmp
import glob
import io
import json
import os
import shlex
import shutil
import subprocess
import sys
import re
import math
import tempfile
import time

# netCDF4 before h5py: both bundle HDF5 and the first loaded wins; the other order breaks
# netCDF4.Dataset() later in the process ("NetCDF: HDF error").
try:
    import netCDF4  # noqa: F401  (load order only)
except ImportError:
    pass
import h5py
from jinja2 import StrictUndefined, Template, UndefinedError
import yaml

from jobgen import JobGen, account_value, user_projects


DESCRIPTION = "Interactive script to prepare a MAGE magnetosphere model run"

# The GR lead-in (spin-up and warm-up) runs in <run>/gr_spinup/; engage imports this constant.
GR_SPINUP_DIR = "gr_spinup"

PARAMETERS_JSON = "makeitso_parameters.json"

LAYOUT_HELP = (
    f"The run is written to pbs.run_directory; a segmented run's GR spin-up runs in {GR_SPINUP_DIR}/. "
    f"Layout and replay: docs/source/makeitso/makeitso.rst.")

RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[33m'
RESET = '\033[0m'

MODES = ["BASIC", "INTERMEDIATE", "EXPERT"]

# The run date/time spelling shared by makeitso, engage, waxrun, tiegcmrun and cda2wind.
RUN_DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S"

KAIJU_ENV_HINT = "source kaiju-private/scripts/setupEnvironment.sh"

# Fall back to this script's directory so `makeitso.py -h` works without KAIJUHOME.
SUPPORT_FILES_DIRECTORY = (
    os.path.join(os.environ["KAIJUHOME"], "scripts", "makeitso") if os.environ.get("KAIJUHOME")
    else os.path.dirname(os.path.abspath(__file__)))

OPTION_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
)

INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.ini")

CONFIG_DIR = os.path.join(SUPPORT_FILES_DIRECTORY, "config")

JSON_INDENT = 4

# Segments must start on a multiple of voltron.restart.dtRes; plan_segments may lower dtRes by
# at most this factor to make them do so.
MAX_DTRES_REDUCTION = 4

# cda2wind refuses a fetch window shorter than this.
CDA2WIND_MIN_WINDOW_SECONDS = 7200

# Every deck names this bare file; a user wind file is linked in under this name.
CDA2WIND_OUTPUT = "bcwind.h5"

# Free-text deck strings. kaiju reads them list-directed, so a value containing '/', ',' or a
# blank is silently truncated (or ignored).
DECK_STRING_OPTIONS = (
    ("gamera", "sim", "runid"),
    ("gamera", "sim", "H5Grid"),
    ("gamera", "ring", "gid"),
    ("gamera", "wind", "tsfile"),
)
_DECK_STRING_UNREADABLE = re.compile(r"[/,\s]")

# Ghost cells on each side of an LFM grid file.
LFM_GRID_GHOSTS = 4

GRID_RADIUS_RTOL = 1e-3

# GAMERA grid types allowed per run kind. H stops voltron unless the deck sets sabotage/isKareem.
RUN_GRID_TYPES = {"gr": ("D", "Q", "O", "H"), "grt": ("D", "Q", "O"), "grx": ("Q",)}

# Fixed values every GRX deck writes (not user options).
GRX_DECK_KEYS = {"tIO": "0.0", "doDeep": "T", "doExcessMap": "T", "nBounce": "1.0"}

# Defaults for arguments not supplied (e.g. when makeitso() is called from code).
args_default = {
    "clobber": False,
    "debug": False,
    "mode": "BASIC",
    "options_path": None,
    "verbose": False,
    "tui": False,
    # engage passes False: it writes PARAMETERS_JSON itself.
    "save_parameters": True,
    "rederive": None,
    # Called as before_write(options) just before anything is written.
    "before_write": None,
}


def require_env(name: str, hint: str = KAIJU_ENV_HINT):
    """Return os.environ[name], or exit with a one-line message if it is unset."""
    value = os.environ.get(name)
    if not value:
        raise SystemExit(f"ERROR: environment variable {name} is not set; {hint}.")
    return value


def as_list(value, tokens: bool = False):
    """Normalise a list-valued option to a list; [None] means no entries.

    Text splits on newlines and ';' (entries are whole shell commands); tokens=True also splits
    on commas and whitespace.
    """
    import ast
    if isinstance(value, (list, tuple)):
        return list(value)
    if value is None:
        return [None]
    s = str(value).strip()
    if s.lower() in ("", "none", "null", "[none]", "[null]"):
        return [None]
    if s.startswith("[") and s.endswith("]"):
        for parse in (ast.literal_eval, json.loads):
            try:
                parsed = parse(s)
            except (ValueError, SyntaxError, TypeError):
                continue
            if isinstance(parsed, (list, tuple)):
                return list(parsed) or [None]
    parts = []
    for line in s.replace(";", "\n").splitlines():
        line = line.strip()
        # strip quotes around one whole entry only ('a' b 'c' keeps its quotes)
        if len(line) >= 2 and line[0] == line[-1] and line[0] in "'\"" and line[0] not in line[1:-1]:
            line = line[1:-1].strip()
        if tokens:
            parts += [w.strip("'\"") for w in line.replace(",", " ").split()]
        elif line:
            parts.append(line)
    return [p for p in parts if p] or [None]


# List options whose entries are single words (as_list(tokens=True)).
TOKEN_LIST_FIELDS = ()


def parse_run_datetime(value, field: str = "date"):
    """Parse RUN_DATETIME_FORMAT to a naive datetime; ValueError naming the field otherwise."""
    if isinstance(value, datetime.datetime):
        return value
    try:
        return datetime.datetime.strptime(str(value).strip(), RUN_DATETIME_FORMAT)
    except (TypeError, ValueError):
        raise ValueError(f"{field} must be given as yyyy-mm-ddThh:mm:ss "
                         f"(e.g. 2013-03-17T00:00:00); got {value!r}") from None


_YES, _NO = ("y", "yes", "true"), ("n", "no", "false")
_TRUE, _FALSE = ("t", "true", "yes"), ("f", "false", "no")


def yn(value):
    """'Y' / 'N' for a yes/no spelling (Y/y/N/n/yes/no/true/false, any case), else None."""
    if isinstance(value, bool):
        return "Y" if value else "N"
    s = str(value).strip().lower() if value is not None else ""
    return "Y" if s in _YES else "N" if s in _NO else None


def tf(value):
    """'T' / 'F' for a Fortran-flag spelling (T/t/F/f/true/false/yes/no, any case), else None."""
    if isinstance(value, bool):
        return "T" if value else "F"
    s = str(value).strip().lower() if value is not None else ""
    return "T" if s in _TRUE else "F" if s in _FALSE else None


_NAME_RE = re.compile(r"^[A-Za-z0-9_.+-]+$")
REQUIRED_NAMES = ("simulation.job_name", "gamera.sim.runid")
# valids checked by the replay rules instead of the type check
VALIDS_BY_RULES = ("pbs.queue",)
_HMS_RE = re.compile(r"^\d+:\d{2}:\d{2}$")
_INT_RE = re.compile(r"^[+-]?\d+$")


def _anchor(value, run_dir):
    """Expand ~ and join a relative path onto run_dir."""
    s = os.path.expanduser(str(value).strip())
    if run_dir is None or os.path.isabs(s) or s.startswith("$"):
        return s
    return os.path.normpath(os.path.join(str(run_dir), s))


def _coerce_option(path, desc, value, run_dir=None, check_files=True):
    """Validate and canonicalise an option value against its description; ValueError if invalid."""
    kind = (desc or {}).get("type")
    if value is None or kind is None:
        return value
    if kind == "lines":
        return as_list(value)
    if isinstance(value, (list, dict)):
        raise ValueError(f"{path}={value!r} is not a single value")
    text = str(value).strip()
    if kind in ("yn", "tf"):
        out = yn(value) if kind == "yn" else tf(value)
        if out is None:
            raise ValueError(f"{path}={value!r} is not {'Y or N' if kind == 'yn' else 'T or F'}")
        value = out
    elif kind in ("int", "float"):
        if isinstance(value, bool):
            raise ValueError(f"{path}={value!r} is not a number")
        try:
            num = float(text.replace("d", "e").replace("D", "e"))
        except ValueError:
            raise ValueError(f"{path}={value!r} is not a number") from None
        if not math.isfinite(num):
            raise ValueError(f"{path}={value!r} is not a finite number")
        if kind == "int":
            if not num.is_integer():
                raise ValueError(f"{path}={value!r} is not a whole number")
            if isinstance(value, str) and not _INT_RE.match(text):
                value = str(int(num))
            elif isinstance(value, float):
                value = int(num)
        elif isinstance(value, str):
            value = text
        lo, hi = desc.get("min"), desc.get("max")
        if lo is not None and (num <= lo if desc.get("exclusive_min") else num < lo):
            raise ValueError(f"{path}={value!r} must be {'>' if desc.get('exclusive_min') else '>='} {lo}")
        if hi is not None and num > hi:
            raise ValueError(f"{path}={value!r} must be <= {hi}")
    elif kind == "name":
        if text.lower() in ("", "none", "null"):
            if path in REQUIRED_NAMES:
                raise ValueError(f"{path} needs a value")
            return value
        if not _NAME_RE.match(text):
            raise ValueError(f"{path}={value!r} may hold only letters, digits and . _ + -")
        value = text
    elif kind == "iso":
        value = parse_run_datetime(value, path).strftime(RUN_DATETIME_FORMAT)
    elif kind == "hms":
        if not _HMS_RE.match(text):
            raise ValueError(f"{path}={value!r} is not HH:MM:SS")
        value = text
    elif kind in ("file", "dir", "path"):
        if text.lower() in ("", "none", "null"):
            return None
        value = _anchor(text, run_dir) if kind != "file" else os.path.expanduser(text)
        if check_files and kind == "file" and not os.path.isfile(value):
            raise ValueError(f"{path}={text!r} does not exist")
        if check_files and kind == "dir" and os.path.isabs(value) and not os.path.isdir(value):
            raise ValueError(f"{path}={text!r} is not a directory")
    else:
        value = text
    valids = desc.get("valids") if path not in VALIDS_BY_RULES else None
    if valids and str(value) not in [str(v) for v in valids]:
        raise ValueError(f"{path}={value!r} is not one of {'|'.join(str(v) for v in valids)}")
    return value


def option_value_error(path, desc, value, run_dir=None):
    """The error message for an invalid value, else None."""
    try:
        _coerce_option(path, desc, value, run_dir)
    except ValueError as e:
        return str(e)
    return None


def normalize_option(path, desc, value, run_dir=None, check_files=True):
    return _coerce_option(path, desc, value, run_dir, check_files)


def normalize_mode(mode):
    """Upper-case and validate a user mode."""
    m = str(mode).strip().upper() if mode is not None else "BASIC"
    if m not in MODES:
        raise ValueError(f"mode must be one of {'|'.join(MODES)} (case-insensitive); got {mode!r}")
    return m


def _import_own_tui_form():
    """Load this directory's tui_form.py by path (tiegcmrun ships one of the same name)."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "makeitso_tui_form", os.path.join(here, "tui_form.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def use_tui(no_tui):
    """True unless --no-tui or there is no terminal."""
    return not no_tui and sys.stdin.isatty() and sys.stdout.isatty()


def create_command_line_parser():
    """Return the makeitso argument parser."""
    parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=LAYOUT_HELP)
    parser.add_argument(
        "--clobber", action="store_true", default=args_default["clobber"],
        help="Overwrite existing options file (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", default=args_default["debug"],
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--mode", default=args_default["mode"], type=str.upper, choices=MODES, metavar="MODE",
        help="User mode: BASIC|INTERMEDIATE|EXPERT, case-insensitive (default: %(default)s)."
    )
    parser.add_argument(
        "--options_path", "-o", default=args_default["options_path"],
        help="Path to JSON file of options (default: %(default)s)."
    )
    parser.add_argument(
        "--rederive", default=args_default["rederive"], metavar="GROUP[,GROUP]",
        help="With -o: re-derive these groups; '--rederive list' names them."
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        default=args_default["verbose"],
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "--no-tui", action="store_true",
        help="Ask the options as linear prompts instead of the wizard."
    )
    parser.add_argument(
        "--tui", action="store_true",
        help="The wizard (the default; kept for old command lines)."
    )
    return parser

# Keys from older options JSONs, refused with what replaces them.
OLD_FORMAT_KEYS = {
    ("provenance",): "no longer written",
    ("voltron", "output", "tsOut"): "use voltron.output.dtCon",
    ("voltron", "coupling", "doDeep"): "model default; a GRX deck sets it itself",
    ("voltron", "spinup", "tIO"): "model default 0.0: restart k is at k*dtRes; a GRX deck sets it itself",
    ("raiju", "cpl", "nBounce"): "model default; a GRX deck sets it itself",
    ("raiju", "plasmasphere", "doExcessMap"): "model default; a GRX deck sets it itself",
    ("gamera", "restart", "resID"): "always gamera.sim.runid",
    ("gamera", "source", "nBounce"): "model default",
    ("pbs", "modules"): "from config/machines.yaml",
    ("pbs", "moduledir"): "from config/machines.yaml",
    ("pbs", "mpiexec_command"): "from config/machines.yaml",
    ("pbs", "ncpus"): "from config/machines.yaml",
    ("pbs", "other"): "use simulation.node_type",
    ("pbs", "environment_variables"): "from config/machines.yaml",
}


def old_format_problem(options: dict):
    """The error message if options carry OLD_FORMAT_KEYS, else None."""
    old = []
    for path, what in OLD_FORMAT_KEYS.items():
        node = options
        for k in path[:-1]:
            node = node.get(k) if isinstance(node, dict) else None
        if isinstance(node, dict) and path[-1] in node:
            old.append(f"{'.'.join(path)} ({what})")
    if not old:
        return None
    return (f"the options JSON was written by an older makeitso: {'; '.join(old)}. Remove "
            f"{'it' if len(old) == 1 else 'them'} from the JSON, or generate a new one.")


# Derived values a -o replay keeps as saved unless --rederive names their group.
REDERIVE_GROUPS = {
    "job": "the machine settings of the pbs block (queue, walltime, job_priority / group_list / "
           "local_modules, ranks and threads per node) and gamera.coupling.blockHalo / "
           "voltron.coupling.doAsyncCoupling, from config/machines.yaml for simulation.hpc_system "
           "(pbs.account_name is kept)",
    "dates": "voltron.time.tFin and pbs.num_segments, from start_date, stop_date and segment_duration",
    "grid": "gamera.sim.H5Grid, the GAMERA decomposition (gamera.iPdir/jPdir/kPdir.N), "
            "pbs.num_helpers, pbs.select2 and voltron.helpers.numHelpers / useHelpers, from "
            "gamera_grid_type",
    "names": "gamera.sim.runid, from job_name",
}


def parse_rederive(value):
    """Parse --rederive ('job,grid', 'all', 'list') to a tuple of group names."""
    if value in (None, "", (), []):
        return ()
    items = value if isinstance(value, (list, tuple)) else str(value).split(",")
    groups = [str(g).strip().lower() for g in items if str(g).strip()]
    if groups == ["list"]:
        return ("list",)
    if "all" in groups:
        return tuple(REDERIVE_GROUPS)
    unknown = [g for g in groups if g not in REDERIVE_GROUPS]
    if unknown or not groups:
        raise ValueError(f"unknown --rederive group {', '.join(unknown) or repr(value)} "
                         f"(groups: {', '.join(REDERIVE_GROUPS)}, all, list)")
    return tuple(dict.fromkeys(groups))


def rederive_groups_text():
    return "\n".join(f"{name:<6} {what}" for name, what in REDERIVE_GROUPS.items())


def rederive_pop(options: dict, groups):
    """Remove the keys of the --rederive groups so the fill step derives them again."""
    pbs = options.get("pbs") or {}
    for group in groups:
        if group == "job":
            for f in MACHINE_PBS_FIELDS:
                pbs.pop(f, None)
            options.get("gamera", {}).get("coupling", {}).pop("blockHalo", None)
            options.get("voltron", {}).get("coupling", {}).pop("doAsyncCoupling", None)
        elif group == "dates":
            options.get("voltron", {}).get("time", {}).pop("tFin", None)
            pbs.pop("num_segments", None)
        elif group == "grid":
            options.get("gamera", {}).get("sim", {}).pop("H5Grid", None)
            for axis in ("iPdir", "jPdir", "kPdir"):
                options.get("gamera", {}).get(axis, {}).pop("N", None)
            for k in ("num_helpers", "select2"):
                pbs.pop(k, None)
            for k in ("numHelpers", "useHelpers"):
                options.get("voltron", {}).get("helpers", {}).pop(k, None)
        elif group == "names":
            options.get("gamera", {}).get("sim", {}).pop("runid", None)
    return options


def load_option_descriptions(path: str = OPTION_DESCRIPTIONS_FILE):
    """Return the option descriptions with the machine facts of config/machines.yaml applied."""
    with open(path, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)

    return apply_system_config(option_descriptions)


def apply_system_config(option_descriptions: dict, hostname: str = None):
    """Set hpc_system valids/default and the account default from config/machines.yaml (in place)."""
    jg = JobGen(CONFIG_DIR, hostname=hostname)
    hpc = option_descriptions.get("simulation", {}).get("hpc_system")
    if hpc is not None:
        hpc["valids"] = list(jg.machines)
        hpc["default"] = jg.name
    account = option_descriptions.get("pbs", {}).get("_common", {}).get("account_name")
    if account is not None and account.get("prompt") is not None:
        account["default"] = account_value(jg.m.get("account_default"))
    return option_descriptions


def account_default(hpc_platform: str, config_dir: str = None):
    """The machine's default batch account from config/machines.yaml, or None."""
    machine = _machines_yaml(config_dir or CONFIG_DIR).get(str(hpc_platform)) or {}
    return account_value(machine.get("account_default"))


def account_choices(hpc_platform: str, config_dir: str = None):
    """The user's project codes where the machine needs a per-user account, else []."""
    machine = _machines_yaml(config_dir or CONFIG_DIR).get(str(hpc_platform)) or {}
    if (not machine.get("account_groups") or not account_required(hpc_platform, config_dir)
            or account_default(hpc_platform, config_dir)):
        return []
    return user_projects()


def account_prompt_default(hpc_platform: str, config_dir: str = None):
    """The machine's account default, else the user's only project code, else None."""
    default = account_default(hpc_platform, config_dir)
    if default:
        return default
    choices = account_choices(hpc_platform, config_dir)
    return choices[0] if len(choices) == 1 else None


def seed_account_prompt(description: dict, hpc_platform: str, config_dir: str = None):
    """Set the account prompt's default and hint choices for the machine."""
    description["default"] = account_prompt_default(hpc_platform, config_dir)
    choices = account_choices(hpc_platform, config_dir)
    if choices:
        description["choices"] = choices
    else:
        description.pop("choices", None)
    return description


def account_required(hpc_platform: str, config_dir: str = None) -> bool:
    """True when the machine's jobs are charged to an account."""
    machine = _machines_yaml(config_dir or CONFIG_DIR).get(str(hpc_platform))
    return bool((machine or {}).get("account_required"))


_MACHINES_YAML = {}


def _machines_yaml(config_dir: str) -> dict:
    """The `machines` mapping of machines.yaml, cached until the file changes."""
    path = os.path.join(config_dir, "machines.yaml")
    key = (path, os.stat(path).st_mtime_ns)
    if key not in _MACHINES_YAML:
        with open(path, encoding="utf-8") as f:
            _MACHINES_YAML[key] = yaml.safe_load(f)["machines"]
    return _MACHINES_YAML[key]


def _blank(value) -> bool:
    return value is None or str(value).strip() in ("", "None", "null")


def account_error(value, hpc_platform: str, where: str = "pbs.account_name",
                  source: str = "the options JSON", config_dir: str = None):
    """The error for a missing account where one is required and has no default, else None."""
    if not _blank(value) or not account_required(hpc_platform, config_dir):
        return None
    if account_default(hpc_platform, config_dir):
        return None
    return f"{where} is required on {hpc_platform} (no shared default): add it to {source}."


def account_prompt_check(hpc_platform: str):
    """Prompt check that refuses a blank account where one is required."""
    def check(value):
        if account_error(value, hpc_platform) is None:
            return None
        return (f"an account is required on {hpc_platform} (config/machines.yaml has no shared "
                f"default): enter your project code")
    return check


def settle_account(options: dict, label: str = "the options JSON", section: str = "pbs",
                   key: str = "account_name", config_dir: str = None):
    """Drop, default or refuse the account of loaded options according to the machine."""
    hpc = str(options["simulation"]["hpc_system"])
    blk = options.setdefault(section, {})
    value = blk.get(key)
    if not account_required(hpc, config_dir):
        if not _blank(value):
            _note(f"{label}: {section}.{key}={value!r} ignored: {hpc} jobs carry no account line "
                  f"(config/machines.yaml account_required is not set).")
        blk[key] = None
        return options
    err = account_error(value, hpc, where=f"{section}.{key}", source=label, config_dir=config_dir)
    if err:
        raise ValueError(err)
    if _blank(value):
        blk[key] = account_default(hpc, config_dir)
        _note(f"{label}: {section}.{key} is blank: the {hpc} default {blk[key]!r} (your login name, "
              f"config/machines.yaml account_default) is used.")
    return options


def queue_names(machine: dict):
    """The machine's queue names, or None if it lists none."""
    queues = machine.get("queues")
    if queues is None:
        return None
    return list(queues)


def queue_limits(machine: dict, queue) -> dict:
    """{ncpus_max, walltime_max} of a queue; {} if none are recorded."""
    queues = machine.get("queues")
    if isinstance(queues, dict):
        return dict(queues.get(queue) or {})
    return {}


def hms_seconds(value):
    """Seconds of an HH:MM[:SS] walltime, or None when it is not one."""
    try:
        parts = [int(p) for p in str(value).strip().split(":")]
    except (TypeError, ValueError):
        return None
    if not 2 <= len(parts) <= 3:
        return None
    parts += [0] * (3 - len(parts))
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


def queue_limit_errors(machine: dict, hpc: str, queue, nodes: int, cores_per_node: int, walltime,
                       what: str = "the job"):
    """Errors for a job exceeding the queue's ncpus_max or walltime_max."""
    ss = _import_segment_solver()
    limits = queue_limits(machine, queue)
    problems = []
    err = ss.ncpus_error(limits.get("ncpus_max"), nodes, cores_per_node, queue, what, hpc)
    if err:
        problems.append(err + ".")
    wall, wall_max = hms_seconds(walltime), hms_seconds(limits.get("walltime_max"))
    if wall is not None and wall_max is not None and wall > wall_max:
        problems.append(f"{what}: walltime {walltime} is more than queue {queue!r} on {hpc} allows "
                        f"({limits['walltime_max']}, config/machines.yaml queues): shorten pbs.walltime "
                        f"(and the segments) or use another queue.")
    return problems


def check_queue_limits(options: dict):
    """Raise ValueError if makeitso's jobs exceed the queue limits."""
    sim, pbs = options["simulation"], options["pbs"]
    hpc = str(sim["hpc_system"])
    jg = JobGen(CONFIG_DIR, machine=hpc, node_type=sim.get("node_type"))
    try:
        nodes = int(pbs["select"]) + int(pbs.get("select2") or 1 + int(pbs.get("num_helpers") or 0))
    except (KeyError, TypeError, ValueError):
        return options
    problems = queue_limit_errors(jg.m, hpc, pbs.get("queue"), nodes, int(jg._cores()),
                                  pbs.get("walltime"), what="the MAGE job (pbs.select + pbs.select2)")
    if problems:
        raise ValueError(" ".join(problems))
    return options


def gr_rate_info(hpc: str, grid) -> dict:
    """The standalone GR model rate on `grid` from config/machines.yaml."""
    rates = (_machines_yaml(CONFIG_DIR).get(str(hpc)) or {}).get("rates") or {}
    return _import_segment_solver().lookup_rate(rates, "gr", grid)


_WALLTIME_ADVISED = set()


def _warn_once(key, text):
    if key not in _WALLTIME_ADVISED:
        _WALLTIME_ADVISED.add(key)
        print(f"{YELLOW}WARNING: {text}{RESET}")


def standalone_walltime_error(hpc: str, queue, walltime, grid, segment_seconds=None, run_seconds=None,
                              announce: bool = False):
    """The error if a segment (or the whole run) does not fit the walltime, else None."""
    if hms_seconds(walltime) is None:
        return None
    info = gr_rate_info(hpc, grid)
    if announce and info.get("note") and (hpc, grid) not in _WALLTIME_ADVISED:
        _WALLTIME_ADVISED.add((hpc, grid))
        _note(info["note"].replace("NOTE: ", "") + f" ({hpc}: config/machines.yaml rates)")
    limits = queue_limits(_machines_yaml(CONFIG_DIR).get(str(hpc)) or {}, queue)
    whole = run_seconds is not None
    return _import_segment_solver().walltime_error(
        info, walltime, limits.get("walltime_max"), segment_seconds, queue=queue, hpc=hpc,
        where=f"grid {grid} on {hpc}", what="the whole run (use_segments=N)" if whole else "segment_duration",
        remedy="answer Y and pick a segment length" if whole else "shorten segment_duration",
        warn=(lambda text: _warn_once((hpc, grid, str(walltime), segment_seconds), text)) if announce else None)


def check_standalone_walltime(options: dict, args: dict):
    """Raise ValueError if a standalone run's segments do not fit the walltime."""
    if "coupling" in args or args.get("grx"):
        return options
    sim, pbs = options["simulation"], options["pbs"]
    try:
        run_s = (parse_run_datetime(sim["stop_date"], "simulation.stop_date")
                 - parse_run_datetime(sim["start_date"], "simulation.start_date")).total_seconds()
    except (KeyError, TypeError, ValueError):
        return options
    segmented = yn(sim.get("use_segments", "N")) == "Y"
    seg = float(sim["segment_duration"]) if segmented else run_s
    # a segment under 1 h only needs room for itself, so skip the 1 h check
    sub_hour = min(seg, run_s) < 3600
    err = ((None if sub_hour else
            standalone_walltime_error(sim["hpc_system"], pbs.get("queue"), pbs.get("walltime"),
                                      sim.get("gamera_grid_type"), announce=True))
           or standalone_walltime_error(sim["hpc_system"], pbs.get("queue"), pbs.get("walltime"),
                                        sim.get("gamera_grid_type"), min(seg, run_s),
                                        run_seconds=None if segmented else run_s, announce=True))
    if err:
        raise ValueError(err)
    return options


def standalone_segment_default(hpc: str, queue, walltime, grid, run_seconds) -> str:
    """Default segment_duration in seconds: the largest whole-hour length that fits the walltime."""
    ss = _import_segment_solver()
    run_h = float(run_seconds) / 3600
    limits = queue_limits(_machines_yaml(CONFIG_DIR).get(str(hpc)) or {}, queue)
    try:
        _, _, fit_h = ss.fit_hours(walltime, limits.get("walltime_max"))
    except (TypeError, ValueError):
        return str(float(run_seconds))
    if run_h < 1:
        return str(float(run_seconds))
    top = max(1, ss.max_segment_hours(fit_h, gr_rate_info(hpc, grid)["k"]))
    if top >= run_h:
        return str(float(run_seconds))
    return str(float(top * 3600))


def machine_names() -> list:
    return list(JobGen(CONFIG_DIR).machines)


def host_machine(hostname: str = None) -> str:
    return JobGen(CONFIG_DIR, hostname=hostname).name


def machine_keyed_default(table, hpc_platform: str):
    """The machine's entry of a per-machine default table, else its "generic" entry."""
    if not isinstance(table, dict):
        return table
    return table.get(str(hpc_platform), table.get("generic"))


_RAW_OPTION_DESCRIPTIONS = {}


def machine_science_default(path: str, hpc_platform: str):
    """machine_keyed_default of the option at 'comp.sect.key'."""
    if "od" not in _RAW_OPTION_DESCRIPTIONS:
        with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
            _RAW_OPTION_DESCRIPTIONS["od"] = json.load(f)
    node = _RAW_OPTION_DESCRIPTIONS["od"]
    for part in path.split("."):
        node = node[part]
    return machine_keyed_default(node.get("default"), hpc_platform)


def synthesize_pbs_descriptions(hpc_platform: str, option_descriptions: dict):
    """Build a machine's PBS option block from the "_machine" template and config/machines.yaml."""
    jg = JobGen(CONFIG_DIR, machine=hpc_platform)
    od = copy.deepcopy(option_descriptions["pbs"]["_machine"])
    od["queue"]["default"] = jg.m.get("queue_default", od["queue"]["default"])
    if jg.m.get("walltime_default"):
        od["walltime"]["default"] = jg.m["walltime_default"]
    # None disables validation; a stale list would make the prompt loop forever
    od["queue"]["valids"] = queue_names(jg.m)
    if "job_priority_default" in jg.m:
        od["job_priority"]["default"] = jg.m["job_priority_default"]
    if "job_priority" in od and jg.m.get("job_priorities"):
        od["job_priority"]["valids"] = list(jg.m["job_priorities"])
    if "job_priority" in od:
        jpv = od["job_priority"].get("valids")
        if jpv is not None and od["job_priority"]["default"] not in jpv:
            od["job_priority"]["valids"] = None
    extra = set(jg.m.get("pbs_prompt_extra", []))
    for f in OPTIONAL_PBS_FIELDS:
        if f not in extra:
            od.pop(f, None)
    return od


def node_types(hpc_platform: str):
    """The machine's node type names ([] when it has one kind)."""
    return list(JobGen(CONFIG_DIR, machine=hpc_platform).m.get("node_types") or [])


def node_type_description(option_descriptions: dict, hpc_platform: str):
    """The simulation.node_type description for a machine, or None if it has one node kind."""
    jg = JobGen(CONFIG_DIR, machine=hpc_platform)
    nts = jg.m.get("node_types")
    if not nts:
        return None
    desc = copy.deepcopy(option_descriptions["simulation"]["node_type"])
    desc["valids"] = list(nts)
    if desc.get("default") in (None, ""):
        desc["default"] = jg.node_type_name
    return desc


def node_type_name(hpc_platform: str, node_type: str = None):
    return JobGen(CONFIG_DIR, machine=hpc_platform, node_type=node_type).node_type_name


def node_cores(hpc_platform: str, node_type: str = None) -> int:
    return int(JobGen(CONFIG_DIR, machine=hpc_platform, node_type=node_type)._cores())


def node_model(hpc_platform: str, node_type: str = None):
    return JobGen(CONFIG_DIR, machine=hpc_platform, node_type=node_type)._node_model()


def launcher(hpc_platform: str) -> str:
    return JobGen(CONFIG_DIR, machine=hpc_platform).m["mpirun"]


# Threads-per-rank field -> the ranks-per-node field its default is derived from.
THREAD_FIELDS = {"ompthreads": "mpiprocs", "helper_ompthreads": "helper_mpiprocs"}


def thread_default(hpc_platform: str, node_type, ranks_per_node):
    """Default OMP threads per rank: cores per node // ranks per node."""
    try:
        ranks = int(ranks_per_node)
    except (TypeError, ValueError):
        return None
    if ranks <= 0:
        return None
    return str(max(1, node_cores(hpc_platform, node_type) // ranks))


def set_thread_default(platform_od: dict, name: str, chosen: dict, hpc_platform: str, node_type):
    """Set the default of a threads field from the ranks per node already chosen."""
    ranks_field = THREAD_FIELDS.get(name)
    if ranks_field is None or name not in platform_od or ranks_field not in chosen:
        return
    if platform_od[name].get("prompt") is None:
        return
    default = thread_default(hpc_platform, node_type, chosen[ranks_field])
    if default is not None:
        platform_od[name]["default"] = default


def check_node_type(options: dict):
    """Drop or refuse simulation.node_type according to the machine's node types."""
    sim = options.get("simulation", {})
    nt = sim.get("node_type")
    if nt in (None, ""):
        sim.pop("node_type", None)
        return options
    hpc = str(sim["hpc_system"])
    known = node_types(hpc)
    if not known:
        _note(f"simulation.node_type={nt!r} ignored: {hpc} has a single node type "
              f"(config/machines.yaml).")
        del sim["node_type"]
    elif nt not in known:
        raise ValueError(f"simulation.node_type={nt!r} is not a node type of {hpc}; "
                         f"known: {'|'.join(known)} (config/machines.yaml).")
    return options


def check_pe_layout(options: dict):
    """Refuse ranks x threads per node exceeding the node type's cores."""
    sim, pbs = options["simulation"], options["pbs"]
    nt = sim.get("node_type")
    cores = node_cores(str(sim["hpc_system"]), nt)
    for ranks, threads in (("mpiprocs", "ompthreads"), ("helper_mpiprocs", "helper_ompthreads")):
        if ranks not in pbs or threads not in pbs:
            continue
        if int(pbs[ranks]) * int(pbs[threads]) > cores:
            where = f"{sim['hpc_system']}" + (f" node type {nt!r}" if nt else "")
            raise ValueError(
                f"pbs.{ranks}={pbs[ranks]} x pbs.{threads}={pbs[threads]} exceeds the {cores} cores "
                f"per node of {where}: lower pbs.{threads} to at most {cores // int(pbs[ranks])}.")
    return options


def derive_select(options: dict):
    """Set pbs.select = iPdir*jPdir*kPdir / mpiprocs; refuse a fractional node count."""
    g = options["gamera"]
    try:
        ni, nj, nk = (int(g[axis]["N"]) for axis in ("iPdir", "jPdir", "kPdir"))
        ranks_per_node = int(options["pbs"]["mpiprocs"])
        select_nodes = ni*nj*nk/ranks_per_node
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        select_nodes = None
    if select_nodes is None or select_nodes <= 0 or int(select_nodes) != select_nodes:
        raise ValueError(f"gamera.iPdir.N x jPdir.N x kPdir.N / pbs.mpiprocs = "
                         f"{g['iPdir']['N']} x {g['jPdir']['N']} x {g['kPdir']['N']} / {options['pbs']['mpiprocs']} "
                         f"is not a whole number of nodes")
    options["pbs"]["select"] = str(int(select_nodes))
    return options


def get_pbs_descriptions(option_descriptions: dict, hpc_platform: str):
    """The machine's PBS option block, synthesized and cached on first use."""
    if hpc_platform not in option_descriptions["pbs"]:
        option_descriptions["pbs"][hpc_platform] = synthesize_pbs_descriptions(
            hpc_platform, option_descriptions)
    return option_descriptions["pbs"][hpc_platform]


# PBS directive fields only some machines use (machines.yaml `pbs_prompt_extra`).
OPTIONAL_PBS_FIELDS = ("job_priority", "group_list", "local_modules")

# Every leaf of these sections is rendered into template.ini.
MODEL_SECTIONS = ("gamera", "voltron", "chimp", "remix", "raiju")

# PBS fields render_pbs reads that have an option_descriptions default.
RENDERED_PBS_FIELDS = ("run_directory", "queue", "walltime", "mpiprocs", "ompthreads",
                       "num_helpers", "helper_mpiprocs", "helper_ompthreads")


def _note(msg: str):
    print(f"{YELLOW}NOTE: {msg}{RESET}")


FILL_NOTE_MAX = 8
FILL_NOTE_SHOWN = 3


def fill_note_text(tool: str, names: list, saved: str) -> str:
    """The note listing the keys an options JSON omits (counted beyond FILL_NOTE_MAX)."""
    if len(names) <= FILL_NOTE_MAX:
        return (f"the {tool} options JSON omits {', '.join(names)}: using the default"
                f"{'' if len(names) == 1 else 's'}.")
    return (f"the {tool} options JSON omits {len(names)} keys ({', '.join(names[:FILL_NOTE_SHOWN])}, ...): "
            f"using the defaults; the full set is saved in {saved}.")


def _keyed_default(default, *keys):
    """Look up a keyed default table by the first key found, else return the default."""
    if isinstance(default, dict):
        for k in keys:
            if k in default:
                return default[k]
        return default.get("generic")
    return default


def ensure_sections(options: dict, option_descriptions: dict = None):
    """Create the sections loaded options omit (in place); refuse a section that is not a dict."""
    if option_descriptions is None:
        if "od" not in _RAW_OPTION_DESCRIPTIONS:
            with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
                _RAW_OPTION_DESCRIPTIONS["od"] = json.load(f)
        option_descriptions = _RAW_OPTION_DESCRIPTIONS["od"]
    paths = [("simulation",), ("pbs",)]
    for comp in MODEL_SECTIONS:
        paths += [(comp,)] + [(comp, sect) for sect in option_descriptions[comp]]
    for path in paths:
        node = options
        for key in path[:-1]:
            node = node[key]
        if not isinstance(node.setdefault(path[-1], {}), dict):
            raise ValueError(f"{'.'.join(path)} must be a JSON object")
    return options


def _raw_option_descriptions() -> dict:
    """option_descriptions.json without machine facts, read once."""
    if "od" not in _RAW_OPTION_DESCRIPTIONS:
        with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
            _RAW_OPTION_DESCRIPTIONS["od"] = json.load(f)
    return _RAW_OPTION_DESCRIPTIONS["od"]


def read_options_json(path: str, shown: str = None) -> dict:
    """Read an options JSON object; ValueError naming the file if that fails."""
    shown = shown or path
    try:
        with open(path, "r", encoding="utf-8") as f:
            text = f.read()
    except FileNotFoundError:
        raise ValueError(f"{shown} does not exist") from None
    except UnicodeDecodeError:
        raise ValueError(f"{shown} is not valid JSON (not UTF-8 text)") from None
    except OSError as e:
        raise ValueError(f"{shown} cannot be read ({e.strerror or e})") from None
    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(f"{shown} is not valid JSON (line {e.lineno} column {e.colno}: {e.msg})") from None
    if not isinstance(data, dict):
        raise ValueError(f"{shown} must be a JSON object")
    return data


def section_problem(options: dict, option_descriptions: dict = None):
    """The error for the first section that is present but not an object, else None."""
    od = option_descriptions or _raw_option_descriptions()
    paths = [("simulation",), ("pbs",)]
    for comp in MODEL_SECTIONS:
        paths += [(comp,)] + [(comp, sect) for sect in od[comp]]
    for path in paths:
        node = options
        for key in path[:-1]:
            node = node.get(key) if isinstance(node, dict) else None
        if isinstance(node, dict) and path[-1] in node and not isinstance(node[path[-1]], dict):
            return f"{'.'.join(path)} must be a JSON object"
    return None


# Derived pbs keys with no option description.
DERIVED_PBS_KEYS = ("select", "select2", "num_segments")

# Obsolete keys that a saved options JSON may still hold and replay with.
PREVIOUS_KEYS = (
    "simulation.LFM_grid_type", "simulation.runid", "pbs.helper_select", "pbs.hpc_system",
    "pbs.num_jobs", "pbs.runid", "gamera.ringknobs.doVClean", "gamera.source.doBounceDT",
    "gamera.source.doWolfLim", "raiju.grid.ThetaU", "raiju.output.doLossExtras",
    "raiju.output.writeGhosts", "remix.conductance.doRamp", "remix.precipitation.alpha",
    "voltron.coupling.rTrc", "voltron.ebsquish.epsSquish", "simulation.bcwind_t0", "rcm")


def add_schema_paths(schema: dict, dotted) -> dict:
    """Add dotted key paths to a key schema (in place)."""
    for path in dotted:
        *parents, leaf = path.split(".")
        node = schema
        for key in parents:
            node = node.setdefault(key, {})
        if isinstance(node, dict):
            node.setdefault(leaf, None)
    return schema


def options_key_schema(option_descriptions: dict = None) -> dict:
    """The keys a makeitso options JSON may hold, as {key: sub-schema or None for a leaf}."""
    od = option_descriptions or _raw_option_descriptions()
    pbs = list(od["pbs"]["_common"]) + list(od["pbs"]["_machine"]) + list(DERIVED_PBS_KEYS)
    schema = {"simulation": dict.fromkeys(od["simulation"]), "pbs": dict.fromkeys(pbs)}
    for comp in MODEL_SECTIONS:
        schema[comp] = {sect: dict.fromkeys(od[comp][sect]) for sect in od[comp]}
    return add_schema_paths(schema, list(PREVIOUS_KEYS) + [".".join(p) for p in OLD_FORMAT_KEYS])


def unknown_keys_problem(options: dict, schema: dict):
    """The error naming keys not in `schema` (keys starting with '_' allowed), else None."""
    import difflib

    def paths(sch, prefix=""):
        out = []
        for k, v in sch.items():
            out.append(prefix + k)
            if isinstance(v, dict):
                out += paths(v, f"{prefix}{k}.")
        return out

    everywhere = paths(schema)
    items = []

    def walk(node, sch, prefix):
        for k, v in node.items():
            if str(k).startswith("_"):
                continue
            if k not in sch:
                close = difflib.get_close_matches(str(k), list(sch), n=1)
                hint = (prefix + close[0]) if close else next(
                    (p for p in everywhere if p.rsplit(".", 1)[-1] == k), None)
                items.append(f"{prefix}{k}" + (f" (did you mean {hint}?)" if hint else ""))
            elif isinstance(sch[k], dict) and isinstance(v, dict):
                walk(v, sch[k], f"{prefix}{k}.")

    walk(options, schema, "")
    if not items:
        return None
    return f"unknown key{'' if len(items) == 1 else 's'} " + "; ".join(items)


def fill_missing_loaded_options(options: dict, option_descriptions: dict, args: dict = None):
    """Fill the entries a loaded (-o) options JSON omits with the prompt defaults.

    Required entries with no default are refused; template.ini would render them as empty.
    """
    ensure_sections(options, option_descriptions)
    grx = bool((args or {}).get("grx"))
    sim = options["simulation"]
    required = ["start_date", "stop_date", "gamera_grid_type"] if not grx else ["gamera_grid_type"]
    missing = [f"simulation.{k}" for k in required if sim.get(k) in (None, "")]
    if missing:
        raise ValueError(f"The run options have no value for {', '.join(missing)}: add "
                         f"{'it' if len(missing) == 1 else 'them'} to the options JSON.")
    err = grid_type_error(sim["gamera_grid_type"], run_component(args))
    if err:
        raise ValueError(err)
    filled = []
    for k in ("job_name", "bcwind_available", "gamera_grid_inner_radius", "gamera_grid_outer_radius",
              "hpc_system", "use_segments"):
        if sim.get(k) in (None, ""):
            sim[k] = option_descriptions["simulation"][k].get("default")
            filled.append(f"simulation.{k}")
    grid = str(sim["gamera_grid_type"])
    hpc = str(sim["hpc_system"])

    pbs = options["pbs"]
    for k, default in (("kaiju_install_directory", lambda: os.environ["KAIJUHOME"]),
                       ("kaiju_build_directory", lambda: kaiju_build_default(run_component(args)))):
        if pbs.get(k) in (None, ""):
            pbs[k] = default()
            filled.append(f"pbs.{k}")
    if "account_name" not in pbs and account_required(hpc) and account_default(hpc):
        pbs["account_name"] = account_default(hpc)
        filled.append("pbs.account_name")
    machine_od = get_pbs_descriptions(option_descriptions, hpc)
    for k, desc in list(option_descriptions["pbs"]["_common"].items()) + list(machine_od.items()):
        if k not in pbs and (k in RENDERED_PBS_FIELDS or k in OPTIONAL_PBS_FIELDS):
            default = _keyed_default(desc.get("default"), grid)
            if k in THREAD_FIELDS:
                default = thread_default(hpc, sim.get("node_type"),
                                         pbs.get(THREAD_FIELDS[k])) or default
            if default is not None:
                pbs[k] = default
                filled.append(f"pbs.{k}")
    if "num_helpers" in pbs and "select2" not in pbs:
        pbs["select2"] = str(1 + int(pbs["num_helpers"]))
    # needs the queue and walltime filled above
    if not grx and sim.get("segment_duration") in (None, ""):
        if yn(sim["use_segments"]) == "Y":
            sim["segment_duration"] = standalone_segment_default(hpc, pbs.get("queue"), pbs.get("walltime"),
                                                                 grid, run_seconds(sim))
        else:
            sim["segment_duration"] = str(run_seconds(sim))
        filled.append("simulation.segment_duration")

    derived = {("gamera", "sim", "H5Grid"): f"lfm{grid}.h5",
               ("gamera", "sim", "runid"): sim["job_name"]}
    if "num_helpers" in pbs:
        derived[("voltron", "helpers", "numHelpers")] = int(pbs["num_helpers"])
        derived[("voltron", "helpers", "useHelpers")] = "T" if int(pbs["num_helpers"]) else "F"
    if not grx:
        derived[("voltron", "time", "tFin")] = run_seconds(sim)
    else:
        t0, t1 = (parse_run_datetime(args["grx"][k], f"grx.{k}") for k in ("t0", "t1"))
        derived[("voltron", "time", "tFin")] = (t1 - t0).total_seconds()
    missing = []
    for comp in MODEL_SECTIONS:
        for sect, fields in option_descriptions[comp].items():
            o = options.setdefault(comp, {}).setdefault(sect, {})
            for k, desc in fields.items():
                if k in o:
                    continue
                if (comp, sect, k) in derived:
                    default = derived[(comp, sect, k)]
                else:
                    default = _keyed_default(desc.get("default"), grid, hpc)
                if default is None:
                    missing.append(f"{comp}.{sect}.{k}")
                    continue
                o[k] = default
                filled.append(f"{comp}.{sect}.{k}")
    if missing:
        raise ValueError(f"The run options have no value for {', '.join(missing)} and there is "
                         f"no default: add {'it' if len(missing) == 1 else 'them'} to the "
                         f"options JSON.")
    if filled:
        _note(fill_note_text("makeitso", filled, PARAMETERS_JSON))
    return options


def _run_dir_of(options: dict, args: dict = None) -> str:
    """The absolute run directory that relative paths are anchored to."""
    grx = (args or {}).get("grx") or {}
    rd = grx.get("run_dir") or (options.get("pbs") or {}).get("run_directory") or "."
    return os.path.abspath(os.path.expanduser(str(rd)))


def validate_options(options: dict, option_descriptions: dict, args: dict = None):
    """Validate and canonicalise every option in place; all problems raise one ValueError."""
    problems = []
    run_dir = _run_dir_of(options, args)

    def walk(descs, values, prefix):
        if not isinstance(values, dict):
            return
        for k, desc in descs.items():
            if not isinstance(desc, dict) or "LEVEL" not in desc or k not in values:
                continue
            path = f"{prefix}.{k}"
            try:
                values[k] = _coerce_option(path, desc, values[k],
                                           None if path == "pbs.run_directory" else run_dir,
                                           check_files=path != "simulation.bcwind_file")
            except ValueError as e:
                problems.append(str(e))

    walk(option_descriptions["simulation"], options.get("simulation"), "simulation")
    hpc = str((options.get("simulation") or {}).get("hpc_system"))
    pbs_desc = dict(option_descriptions["pbs"]["_common"])
    if hpc in (option_descriptions["simulation"]["hpc_system"].get("valids") or [hpc]):
        pbs_desc.update(get_pbs_descriptions(option_descriptions, hpc))
    walk(pbs_desc, options.get("pbs"), "pbs")
    for comp in MODEL_SECTIONS:
        for sect, fields in option_descriptions[comp].items():
            walk(fields, (options.get(comp) or {}).get(sect), f"{comp}.{sect}")
    stride = ((options.get("voltron") or {}).get("coupling") or {}).get("qkSquishStride")
    try:
        n = int(stride)
        if n >= 1 and n & (n - 1):
            problems.append(f"voltron.coupling.qkSquishStride={stride!r} must be a power of two")
    except (TypeError, ValueError):
        pass
    sim = options.get("simulation") or {}
    try:
        rin, rout = float(sim["gamera_grid_inner_radius"]), float(sim["gamera_grid_outer_radius"])
        if rin >= rout:
            problems.append(f"simulation.gamera_grid_inner_radius={sim['gamera_grid_inner_radius']!r} must be "
                            f"below gamera_grid_outer_radius={sim['gamera_grid_outer_radius']!r}")
    except (KeyError, TypeError, ValueError):
        pass
    if problems:
        raise ValueError("\n".join(problems))
    return options


def check_helper_layout(options: dict):
    """Refuse select2 / numHelpers / useHelpers that do not follow pbs.num_helpers."""
    pbs, helpers = options.get("pbs") or {}, (options.get("voltron") or {}).get("helpers") or {}
    try:
        nh = int(pbs["num_helpers"])
    except (KeyError, TypeError, ValueError):
        return options
    problems = []
    if "select2" in pbs and str(pbs["select2"]) != str(1 + nh):
        problems.append(f"pbs.select2={pbs['select2']!r} is not 1 + pbs.num_helpers ({1 + nh})")
    if "numHelpers" in helpers and str(helpers["numHelpers"]) != str(nh):
        problems.append(f"voltron.helpers.numHelpers={helpers['numHelpers']!r} is not pbs.num_helpers ({nh})")
    if "useHelpers" in helpers and tf(helpers["useHelpers"]) != ("T" if nh else "F"):
        problems.append(f"voltron.helpers.useHelpers={helpers['useHelpers']!r} does not follow "
                        f"pbs.num_helpers={nh} ({'T' if nh else 'F'})")
    if problems:
        raise ValueError("; ".join(problems) + " (pass --rederive grid to follow pbs.num_helpers)")
    return options


def kaiju_build_default(component: str) -> str:
    return os.path.join(os.environ["KAIJUHOME"],
                        {"gr": "build_mpi", "grt": "build_grt", "grx": "build_grx"}[component])


# PBS fields re-derived from config/machines.yaml by --rederive job.
MACHINE_PBS_FIELDS = ("queue", "walltime", "job_priority", "group_list", "local_modules",
                      "mpiprocs", "ompthreads", "helper_mpiprocs", "helper_ompthreads")


def machine_block_evidence(pbs: dict, platform_od: dict, hpc: str):
    """Reasons a loaded PBS block looks made for another machine; a queue no machine has raises."""
    evidence = []
    valids = platform_od.get("queue", {}).get("valids")
    queue = pbs.get("queue")
    if valids and not _blank(queue) and queue not in valids:
        others = {q for name, m in _machines_yaml(CONFIG_DIR).items() if name != hpc
                  for q in (queue_names(m) or [])}
        if queue not in others:
            raise ValueError(f"pbs.queue={queue!r} is not a queue of {hpc} (config/machines.yaml queues: "
                             f"{'|'.join(valids)}): set pbs.queue to one of them.")
        evidence.append(f"pbs.queue {queue!r} is not a {hpc} queue")
    for f in OPTIONAL_PBS_FIELDS:
        if f in pbs and f not in platform_od:
            evidence.append(f"pbs.{f} is not a {hpc} field")
    return evidence


def apply_grid_layout(options: dict, option_descriptions: dict):
    """Set the values derived from the grid type and pbs.num_helpers."""
    grid = str(options["simulation"]["gamera_grid_type"])
    gamera, pbs = options["gamera"], options["pbs"]
    gamera["sim"]["H5Grid"] = f"lfm{grid}.h5"
    for axis in ("iPdir", "jPdir", "kPdir"):
        gamera[axis]["N"] = option_descriptions["gamera"][axis]["N"]["default"][grid]
    num_helpers = int(pbs["num_helpers"])
    pbs["select2"] = str(1 + num_helpers)
    helpers = options["voltron"]["helpers"]
    helpers["useHelpers"] = "T" if num_helpers else "F"
    helpers["numHelpers"] = num_helpers
    return options


def run_seconds(sim: dict) -> float:
    """stop_date - start_date in seconds; ValueError unless positive."""
    t0 = parse_run_datetime(sim["start_date"], "simulation.start_date")
    t1 = parse_run_datetime(sim["stop_date"], "simulation.stop_date")
    if t1 <= t0:
        raise ValueError(f"stop_date ({sim['stop_date']}) must be after start_date "
                         f"({sim['start_date']}).")
    return (t1 - t0).total_seconds()


def replay_hints(options: dict, option_descriptions: dict, args: dict) -> list:
    """Hints for replayed derived values that do not follow their inputs, naming the --rederive group."""
    sim, pbs = options["simulation"], options["pbs"]
    hpc, grid = str(sim["hpc_system"]), str(sim["gamera_grid_type"])
    hints = []
    evidence = machine_block_evidence(pbs, get_pbs_descriptions(option_descriptions, hpc), hpc)
    if evidence:
        hints.append(f"the pbs block looks made for another machine than simulation.hpc_system={hpc!r} "
                     f"({'; '.join(evidence)}): pass --rederive job to follow the machine")
    if not args.get("grx"):
        tfin = options.get("voltron", {}).get("time", {}).get("tFin")
        try:
            same = float(tfin) == run_seconds(sim)
        except (TypeError, ValueError):
            same = False
        if not same:
            hints.append(f"voltron.time.tFin={tfin!r} is not stop_date - start_date: pass --rederive "
                         f"dates to follow the dates")
    h5 = options.get("gamera", {}).get("sim", {}).get("H5Grid")
    helpers = options.get("voltron", {}).get("helpers", {})
    try:
        nh = int(pbs.get("num_helpers"))
        helpers_off = str(pbs.get("select2")) != str(1 + nh) or str(helpers.get("numHelpers")) != str(nh)
    except (TypeError, ValueError):
        helpers_off = False
    if h5 != f"lfm{grid}.h5" or helpers_off:
        hints.append(f"gamera.sim.H5Grid / pbs.select2 / voltron.helpers do not follow "
                     f"gamera_grid_type={grid!r} and pbs.num_helpers: pass --rederive grid to follow them")
    if options.get("gamera", {}).get("sim", {}).get("runid") != sim.get("job_name"):
        hints.append(f"gamera.sim.runid is not simulation.job_name={sim.get('job_name')!r}: pass "
                     f"--rederive names to follow it")
    return hints


def prepare_replay(options: dict, option_descriptions: dict, args: dict):
    """Fill, check and annotate a loaded (-o) options JSON for replay."""
    fill_missing_loaded_options(options, option_descriptions, args)
    check_node_type(options)
    settle_account(options)
    for hint in replay_hints(options, option_descriptions, args):
        _note(hint + ".")
    return options


def update_option_descriptions(option_descriptions: dict, args: dict):
    """Make each value a caller (engage) passes in `args` the option's default, not asked."""
    if "simulation" in args:
        simulation = args["simulation"]
        for k in simulation:
            od = option_descriptions["simulation"].get(k)
            if od is None:
                continue
            od["prompt"] = None
            od["default"] = simulation[k]

    # A coupled run starts gr_warm_up_time before engage's start date.
    if "coupling" in args:
        coupling = args["coupling"]
        simulation = args["simulation"]
        gr_warm_up_time = float(coupling["gr_warm_up_time"])
        dt = datetime.timedelta(seconds=gr_warm_up_time)
        start_date = simulation["start_date"]
        t0 = parse_run_datetime(start_date, "simulation.start_date")
        t0 -= dt
        start_date = datetime.datetime.isoformat(t0)
        option_descriptions["simulation"]["start_date"]["prompt"] = None
        option_descriptions["simulation"]["start_date"]["default"] = (
            start_date
        )
        if not args.get("grx"):
            t1 = parse_run_datetime(simulation["stop_date"], "simulation.stop_date")
            option_descriptions["voltron"]["time"]["tFin"].update(prompt=None, default=(t1 - t0).total_seconds())
            option_descriptions["simulation"]["use_segments"].update(prompt=None, default="Y")

    # pbs.num_segments is not taken: makeitso derives it over the warm-up-extended span.
    if "pbs" in args:
        pbs = args["pbs"]
        for k in ["account_name", "kaiju_install_directory",
                 "kaiju_build_directory", "run_directory"]:
            od = option_descriptions["pbs"]["_common"][k]
            od["prompt"] = None
            od["default"] = pbs.get(k)

        hpc_platform = args["simulation"]["hpc_system"]
        platform_od = get_pbs_descriptions(option_descriptions, hpc_platform)
        for k in args["pbs"]:
            if k in platform_od:
                platform_od[k]["prompt"] = None
                platform_od[k]["default"] = pbs[k]
        if "local_modules" in platform_od and "local_modules" in pbs:
            platform_od["local_modules"]["default"] = pbs["local_modules"]
    return option_descriptions


def get_run_option(name: str, description: dict, mode: str, check=None, path: str = None,
                   run_dir: str = None):
    """Prompt the user for a single run option; return the default if it is not asked.

    check(value) returns an error message for a typed value (asked again), or None.
    """
    level = description["LEVEL"]
    prompt = description.get("prompt", None)
    default = description.get("default", None)
    valids = description.get("valids", None)

    if prompt is None:
        return default

    if mode == "BASIC" and level in ["INTERMEDIATE", "EXPERT"]:
        return default
    if mode == "INTERMEDIATE" and level in ["EXPERT"]:
        return default

    if valids is not None:
        vs = "|".join(valids)
        prompt += f" ({vs})"
    elif description.get("choices"):
        # suggestions only: any other value is accepted
        prompt += f" (your projects: {'|'.join(description['choices'])})"

    if default is not None:
        prompt += f" [{GREEN}{default}{RESET}]"

    ok = False
    while not ok:

        option_value = input(f"{prompt}: ")

        if option_value == "":
            option_value = default
        # canonical Y/N or T/F before the valids check
        kind = description.get("type")
        if kind in ("yn", "tf") and option_value is not None:
            option_value = (yn(option_value) if kind == "yn" else tf(option_value)) or option_value

        if valids is not None and option_value not in valids:
            print(f"{YELLOW}Invalid value for option {name}: {option_value}!{RESET}")
            continue
        err = option_value_error(path or name, description, option_value, run_dir)
        if not err and check is not None:
            err = check(option_value)
        if err:
            print(f"{YELLOW}Invalid value for option {name}: {err}{RESET}")
            continue
        option_value = normalize_option(path or name, description, option_value, run_dir)

        ok = True

    # not str()-ed: render_pbs would iterate a str one character per line
    if option_value is None or isinstance(option_value, list):
        return option_value
    if isinstance(default, list):
        return as_list(option_value, tokens=name in TOKEN_LIST_FIELDS)
    return str(option_value)


def fetch_bcwind_time_range(bcwind_path: str):
    """(start_date, stop_date) of a bcwind file: its T=0 (gamera's run start) and its last UT."""
    info = bcwind_file_info(bcwind_path)
    start_date = info["t0"] or info["w0"].strftime(RUN_DATETIME_FORMAT)
    return start_date, info["w1"].strftime(RUN_DATETIME_FORMAT)


def bcwind_file_info(bcwind_path: str) -> dict:
    """First/last UT (w0, w1), median sample spacing in s, and T=0 (t0, or None) of a wind file."""
    with h5py.File(bcwind_path, "r") as f:
        if "UT" not in f or f["UT"].shape[0] == 0:
            raise ValueError(f"bcwind file {bcwind_path} has no usable 'UT' time dataset.")
        ut = [x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) for x in f["UT"][:]]
        t_first = float(f["T"][0]) if "T" in f and f["T"].shape[0] > 0 else None
    times = [datetime.datetime.fromisoformat(u.replace(" ", "T")) for u in ut]
    steps = sorted((b - a).total_seconds() for a, b in zip(times, times[1:]))
    spacing = steps[len(steps) // 2] if steps else 0.0
    t0 = None
    if t_first is not None:
        t0 = (times[0] - datetime.timedelta(seconds=round(t_first))).strftime(RUN_DATETIME_FORMAT)
    return {"w0": times[0], "w1": times[-1], "spacing": spacing, "t0": t0}


def prompt_user_for_run_options(option_descriptions: dict, args: dict):
    """Prompt for makeitso's options and return the nested options dict."""
    mode = args["mode"]

    print("\n")
    print("Instructions:")
    print(f"-> Default selected input parameter is given in {GREEN}GREEN{RESET}")
    print(f"-> Warnings and information are given in {YELLOW}YELLOW{RESET}")
    print(f"-> Errors are given in {RED}RED{RESET}")
    print("-> Valid values (if any) are given in brackets eg. (value1 | value2 | value3)")
    print("\n")
    print("Run Options:")

    options = {}

    o = options["simulation"] = {}
    od = option_descriptions["simulation"]

    for on in ["job_name"]:
        o[on] = get_run_option(on, od[on], mode, path=f"simulation.{on}")

    for on in ["bcwind_available"]:
        o[on] = get_run_option(on, od[on], mode, path=f"simulation.{on}")
    if yn(o["bcwind_available"]) == "Y":
        for on in ["bcwind_file"]:
            o[on] = get_run_option(on, od[on], mode, path=f"simulation.{on}")
        if not os.path.isfile(o["bcwind_file"]):
            raise FileNotFoundError(
                f"bcwind file {o['bcwind_file']!r} (simulation.bcwind_file) does not exist.")
    if yn(o["bcwind_available"]) == "Y" and "coupling" in args:
        # a coupled run keeps engage's window; the file's T=0 is checked later
        for on in ["start_date", "stop_date"]:
            o[on] = get_run_option(on, od[on], mode, path=f"simulation.{on}")
    elif yn(o["bcwind_available"]) == "Y":
        # the run dates are the wind file's own time range
        start_date, stop_date = fetch_bcwind_time_range(o["bcwind_file"])
        print(f"WARNING: bcwind_available=Y: the run dates are taken from the bcwind file "
              f"'{o['bcwind_file']}' [{start_date} .. {stop_date}], replacing any other dates "
              f"(bcwind_available=N keeps the given dates).")
        o["start_date"] = start_date
        o["stop_date"] = stop_date
    else:
        # also the window of the wind file cda2wind fetches from CDAWeb
        for on in ["start_date", "stop_date"]:
            o[on] = get_run_option(on, od[on], mode, path=f"simulation.{on}")

    # Duration includes engage's warm-up but not the spin-up (t < 0).
    t1 = parse_run_datetime(o["start_date"], "simulation.start_date")
    t2 = parse_run_datetime(o["stop_date"], "simulation.stop_date")
    o["start_date"] = start_date = t1.strftime(RUN_DATETIME_FORMAT)
    o["stop_date"] = stop_date = t2.strftime(RUN_DATETIME_FORMAT)
    if t2 <= t1:
        raise ValueError(
            f"stop_date ({stop_date}) must be after start_date ({start_date})."
        )
    if yn(o["bcwind_available"]) != "Y" and not args.get("grx"):
        err = cda2wind_window_error(start_date, stop_date)
        if err:
            raise ValueError(err)
    simulation_duration = (t2 - t1).total_seconds()
    standalone = "coupling" not in args and not args.get("grx")

    # the segment is asked after the walltime, below
    for on in ["gamera_grid_type", "gamera_grid_inner_radius",
               "gamera_grid_outer_radius", "hpc_system"]:
        o[on] = get_run_option(on, od[on], mode, path=f"simulation.{on}")
    node_type_od = node_type_description(option_descriptions, o["hpc_system"])
    if node_type_od is not None:
        o["node_type"] = get_run_option("node_type", node_type_od, mode, path="simulation.node_type")

    o = options["pbs"] = {}

    od = option_descriptions["pbs"]["_common"]

    hpc_for_account = options["simulation"]["hpc_system"]
    if od["account_name"].get("prompt") is not None:
        seed_account_prompt(od["account_name"], hpc_for_account)
        if not account_required(hpc_for_account):
            od["account_name"]["prompt"] = None
            od["account_name"]["default"] = None
    if "default" not in od.get("kaiju_install_directory", {}):
        od["kaiju_install_directory"]["default"] = os.environ["KAIJUHOME"]
    if "default" not in od.get("kaiju_build_directory", {}):
        od["kaiju_build_directory"]["default"] = kaiju_build_default(run_component(args))
    for on in od:
        if on == "num_segments":
            # derived below; the key keeps its place
            o[on] = None
            continue
        check = account_prompt_check(hpc_for_account) if on == "account_name" else None
        run_dir = (os.path.abspath(os.path.expanduser(str(o["run_directory"])))
                   if o.get("run_directory") not in (None, "") else None)
        o[on] = get_run_option(on, od[on], mode, check=check, path=f"pbs.{on}", run_dir=run_dir)

    hpc_platform = options["simulation"]["hpc_system"]
    gamera_grid_type = options["simulation"]["gamera_grid_type"]
    od = get_pbs_descriptions(option_descriptions, hpc_platform)
    od["num_helpers"]["default"] = _keyed_default(od["num_helpers"]["default"], gamera_grid_type)
    node_type = options["simulation"].get("node_type")
    for on in od:
        set_thread_default(od, on, o, hpc_platform, node_type)
        check = None
        if on == "walltime" and standalone:
            check = (lambda v: standalone_walltime_error(hpc_platform, o.get("queue"), v, gamera_grid_type,
                                                         announce=True))
        o[on] = get_run_option(on, od[on], mode, check=check, path=f"pbs.{on}")

    o = options["simulation"]
    od = option_descriptions["simulation"]
    wall, queue = options["pbs"].get("walltime"), options["pbs"].get("queue")
    if "default" not in od["segment_duration"]:
        od["segment_duration"]["default"] = (
            standalone_segment_default(hpc_platform, queue, wall, gamera_grid_type, simulation_duration)
            if standalone else str(simulation_duration))
    for on in ["use_segments"]:
        o[on] = get_run_option(on, od[on], mode, path=f"simulation.{on}")
    if yn(o["use_segments"]) == "Y":
        seg_check = None
        if standalone:
            def seg_check(v):
                try:
                    if float(v) <= 0:
                        raise ValueError
                except (TypeError, ValueError):
                    return None
                return standalone_walltime_error(hpc_platform, queue, wall, gamera_grid_type, float(v),
                                                 announce=True)
        for on in ["segment_duration"]:
            o[on] = get_run_option(on, od[on], mode, check=seg_check, path=f"simulation.{on}")
        try:
            if float(o["segment_duration"]) <= 0:
                raise ValueError
        except ValueError:
            raise ValueError(
                f"segment_duration must be a positive number, got: "
                f"{o['segment_duration']!r}"
            )
    else:
        o["segment_duration"] = str(simulation_duration) if standalone else od["segment_duration"]["default"]
        if standalone:
            err = standalone_walltime_error(hpc_platform, queue, wall, gamera_grid_type,
                                            simulation_duration, run_seconds=simulation_duration, announce=True)
            if err:
                raise ValueError(err)

    # simulation segments (t > 0), excluding the spin-up
    num_segments = 1
    if yn(o["use_segments"]) == "Y":
        num_segments = math.ceil(simulation_duration/float(o["segment_duration"]))
    options["pbs"]["num_segments"] = str(num_segments)
    o = options["pbs"]

    # second chunk: 1 voltron node + the helpers
    num_helpers = int(o["num_helpers"])
    select2 = 1 + num_helpers
    o["select2"] = str(select2)

    options["gamera"] = {}

    o = options["gamera"]["sim"] = {}
    od = option_descriptions["gamera"]["sim"]
    od["H5Grid"]["default"] = f"lfm{options['simulation']['gamera_grid_type']}.h5"
    od["runid"]["default"] = options["simulation"]["job_name"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.sim.{on}")

    o = options["gamera"]["floors"] = {}
    od = option_descriptions["gamera"]["floors"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.floors.{on}")

    o = options["gamera"]["timestep"] = {}
    od = option_descriptions["gamera"]["timestep"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.timestep.{on}")

    o = options["gamera"]["restart"] = {}
    od = option_descriptions["gamera"]["restart"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.restart.{on}")

    o = options["gamera"]["physics"] = {}
    od = option_descriptions["gamera"]["physics"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.physics.{on}")

    o = options["gamera"]["ring"] = {}
    od = option_descriptions["gamera"]["ring"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.ring.{on}")

    o = options["gamera"]["wind"] = {}
    od = option_descriptions["gamera"]["wind"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.wind.{on}")
    # a user wind file is linked in as bcwind.h5: kaiju cannot read a path here
    if yn(options["simulation"].get("bcwind_available", "N")) == "Y":
        o["tsfile"] = CDA2WIND_OUTPUT

    o = options["gamera"]["source"] = {}
    od = option_descriptions["gamera"]["source"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.source.{on}")

    o = options["gamera"]["iPdir"] = {}
    od = option_descriptions["gamera"]["iPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.iPdir.{on}")

    o = options["gamera"]["jPdir"] = {}
    od = option_descriptions["gamera"]["jPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.jPdir.{on}")

    o = options["gamera"]["kPdir"] = {}
    od = option_descriptions["gamera"]["kPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.kPdir.{on}")

    o = options["gamera"]["coupling"] = {}
    od = option_descriptions["gamera"]["coupling"]
    od["blockHalo"]["default"] = machine_keyed_default(od["blockHalo"]["default"], hpc_platform)
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"gamera.coupling.{on}")

    options["voltron"] = {}

    o = options["voltron"]["time"] = {}
    od = option_descriptions["voltron"]["time"]
    od["tFin"]["default"] = simulation_duration
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"voltron.time.{on}")

    o = options["voltron"]["spinup"] = {}
    od = option_descriptions["voltron"]["spinup"]
    component = run_component(args)
    for on in od:
        # a coupled GRT run needs the spin-up segment
        check = coupled_dospin_prompt_check if (on == "doSpin" and component == "grt") else None
        o[on] = get_run_option(on, od[on], mode, check=check, path=f"voltron.spinup.{on}")

    o = options["voltron"]["output"] = {}
    od = option_descriptions["voltron"]["output"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"voltron.output.{on}")

    o = options["voltron"]["coupling"] = {}
    od = option_descriptions["voltron"]["coupling"]
    od["doAsyncCoupling"]["default"] = machine_keyed_default(od["doAsyncCoupling"]["default"],
                                                             hpc_platform)
    # coupled runs need a whole-second dtCouple that divides 60
    coupled = "coupling" in args or bool(args.get("grx"))
    for on in od:
        check = dtcouple_prompt_check if (coupled and on == "dtCouple") else None
        o[on] = get_run_option(on, od[on], mode, check=check, path=f"voltron.coupling.{on}")

    o = options["voltron"]["restart"] = {}
    od = option_descriptions["voltron"]["restart"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"voltron.restart.{on}")

    o = options["voltron"]["imag"] = {}
    od = option_descriptions["voltron"]["imag"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"voltron.imag.{on}")

    o = options["voltron"]["helpers"] = {}
    od = option_descriptions["voltron"]["helpers"]
    od["numHelpers"]["default"] = num_helpers
    od["useHelpers"]["default"] = "T" if num_helpers else "F"
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"voltron.helpers.{on}")

    options["chimp"] = {}

    o = options["chimp"]["units"] = {}
    od = option_descriptions["chimp"]["units"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"chimp.units.{on}")

    o = options["chimp"]["fields"] = {}
    od = option_descriptions["chimp"]["fields"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"chimp.fields.{on}")

    o = options["chimp"]["domain"] = {}
    od = option_descriptions["chimp"]["domain"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"chimp.domain.{on}")

    o = options["chimp"]["tracer"] = {}
    od = option_descriptions["chimp"]["tracer"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"chimp.tracer.{on}")

    options["remix"] = {}

    o = options["remix"]["conductance"] = {}
    od = option_descriptions["remix"]["conductance"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"remix.conductance.{on}")

    o = options["remix"]["precipitation"] = {}
    od = option_descriptions["remix"]["precipitation"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"remix.precipitation.{on}")

    options["raiju"] = {}

    options["raiju"]["output"] = {}
    o = options["raiju"]["output"]
    od = option_descriptions["raiju"]["output"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.output.{on}")

    options["raiju"]["grid"] = {}
    o = options["raiju"]["grid"]
    od = option_descriptions["raiju"]["grid"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.grid.{on}")

    options["raiju"]["domain"] = {}
    o = options["raiju"]["domain"]
    od = option_descriptions["raiju"]["domain"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.domain.{on}")

    options["raiju"]["config"] = {}
    o = options["raiju"]["config"]
    od = option_descriptions["raiju"]["config"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.config.{on}")

    options["raiju"]["plasmasphere"] = {}
    o = options["raiju"]["plasmasphere"]
    od = option_descriptions["raiju"]["plasmasphere"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.plasmasphere.{on}")

    options["raiju"]["losses"] = {}
    o = options["raiju"]["losses"]
    od = option_descriptions["raiju"]["losses"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.losses.{on}")

    options["raiju"]["cpl"] = {}
    o = options["raiju"]["cpl"]
    od = option_descriptions["raiju"]["cpl"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.cpl.{on}")

    options["raiju"]["fluidIn1"] = {}
    o = options["raiju"]["fluidIn1"]
    od = option_descriptions["raiju"]["fluidIn1"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode, path=f"raiju.fluidIn1.{on}")

    return options


def validate_grid_file(path : str, grid_type : str):
    """True when the LFM grid file has the X/Y/Z shape of `grid_type`."""
    grid_shapes = {
        'D': (73, 57, 57),
        'Q': (137, 105, 105),
        'O': (265, 201, 201),
        'H': (521, 393, 393),
    }

    is_valid = False
    if grid_type not in grid_shapes:
        raise ValueError(
            f"Unknown grid type {grid_type!r}; expected one of {sorted(grid_shapes)}."
        )
    with h5py.File(path, 'r') as hf:
        if not all(d in hf for d in ("X", "Y", "Z")):
            raise ValueError(
                f"Grid file {path} is missing one of the X/Y/Z datasets."
            )
        if (hf['X'].shape == grid_shapes[grid_type] and
            hf['Y'].shape == grid_shapes[grid_type] and
            hf['Z'].shape == grid_shapes[grid_type]):
            is_valid = True

    return is_valid


def grid_file_radii(path: str):
    """(rin, rout) in Earth radii that an LFM grid file was built with (genLFM -Rin/-Rout)."""
    ng = LFM_GRID_GHOSTS
    with h5py.File(path, 'r') as hf:
        x, y, z = (hf[k][:, :, ng] for k in ("X", "Y", "Z"))
        rin = float(((x**2 + y**2 + z**2) ** 0.5).mean())
        rout = float(hf["X"][:, :, -ng - 1].max())
    return rin, rout


def requested_grid_radii(options: dict):
    """(rin, rout) the options ask for as floats; rout 0 means genLFM's default boundary."""
    sim = options["simulation"]
    keys = ("gamera_grid_inner_radius", "gamera_grid_outer_radius")
    try:
        rin, rout = (float(sim[k]) for k in keys)
    except (KeyError, TypeError, ValueError):
        raise ValueError(
            f"simulation.gamera_grid_inner_radius / gamera_grid_outer_radius must be numbers "
            f"(Earth radii), got {sim.get(keys[0])!r} / {sim.get(keys[1])!r}.")
    if rin <= 0 or rout < 0 or (rout > 0 and rout <= rin):
        raise ValueError(
            f"GAMERA grid radii must satisfy 0 < inner radius < outer radius (outer radius 0 = "
            f"genLFM's default), got inner {rin} / outer {rout}.")
    return rin, rout


def cda2wind_window_error(start_date, stop_date):
    """The error if the window is too short for cda2wind, else None."""
    t0 = parse_run_datetime(start_date, "simulation.start_date")
    t1 = parse_run_datetime(stop_date, "simulation.stop_date")
    if (t1 - t0).total_seconds() >= CDA2WIND_MIN_WINDOW_SECONDS:
        return None
    return (f"The run window {t0.strftime(RUN_DATETIME_FORMAT)} .. "
            f"{t1.strftime(RUN_DATETIME_FORMAT)} is shorter than "
            f"{CDA2WIND_MIN_WINDOW_SECONDS // 3600} hours: with bcwind_available=N the solar wind "
            f"file is fetched by cda2wind, which needs at least "
            f"{CDA2WIND_MIN_WINDOW_SECONDS // 3600} hours. Lengthen the run or supply a wind file "
            f"(bcwind_available=Y).")


def existing_grid_file_reusable(options: dict, clobber: bool = False) -> bool:
    """True if lfm<g>.h5 in the cwd matches the requested grid; other radii raise unless clobber."""
    grid_type = str(options["simulation"]["gamera_grid_type"])
    grid_file = f"lfm{grid_type}.h5"
    if not os.path.exists(grid_file):
        return False
    if not validate_grid_file(grid_file, grid_type):
        raise TypeError(f'Invalid grid file {grid_file} found, aborting.')
    rin, rout = requested_grid_radii(options)
    have_rin, have_rout = grid_file_radii(grid_file)
    same = math.isclose(have_rin, rin, rel_tol=GRID_RADIUS_RTOL) and (
        rout == 0 or math.isclose(have_rout, rout, rel_tol=GRID_RADIUS_RTOL))
    if same:
        return True
    msg = (f"Existing grid file {grid_file} was built with inner radius {have_rin:g} / sunward "
           f"outer radius {have_rout:g}, but the options ask for {rin:g} / {rout:g}")
    if clobber:
        print(f"{YELLOW}WARNING: {msg}; regenerating it (--clobber).{RESET}")
        return False
    raise ValueError(f"{msg}. Remove it, or rerun with --clobber to regenerate it.")


def check_preprocessing_inputs(options: dict, args: dict):
    """Refuse grid radii and solar wind inputs that preprocessing would fail on (run in the run dir)."""
    requested_grid_radii(options)
    sim = options["simulation"]
    if yn(sim.get("bcwind_available", "N")) == "N":
        err = cda2wind_window_error(sim["start_date"], sim["stop_date"])
        if err:
            raise ValueError(err)
        return
    if sim.get("bcwind_file"):
        err = (bcwind_t0_error(sim["bcwind_file"], sim["start_date"], args.get("coupling"))
               or bcwind_window_error(sim["bcwind_file"], sim["start_date"], sim["stop_date"]))
    else:
        here = os.path.join(os.getcwd(), CDA2WIND_OUTPUT)
        if not os.path.isfile(CDA2WIND_OUTPUT):
            raise ValueError(
                f"simulation.bcwind_available=Y but simulation.bcwind_file is not set and the run "
                f"directory {os.getcwd()} has no {CDA2WIND_OUTPUT}: every deck names "
                f"{CDA2WIND_OUTPUT} (gamera.wind.tsfile) and nothing would create it. Set "
                f"simulation.bcwind_file to the solar wind file (makeitso links it into the run "
                f"directory as {CDA2WIND_OUTPUT}), put the file at {here}, or use "
                f"bcwind_available=N to fetch it with cda2wind.")
        err = (bcwind_t0_error(CDA2WIND_OUTPUT, sim["start_date"], args.get("coupling"))
               or bcwind_window_error(CDA2WIND_OUTPUT, sim["start_date"], sim["stop_date"], name=here))
    if err:
        raise ValueError(err)


def check_bcwind_t0(options: dict, args: dict, start_date=None):
    """Raise bcwind_t0_error for simulation.bcwind_file before the run directory is created."""
    sim = options.get("simulation", {})
    p = sim.get("bcwind_file")
    if yn(sim.get("bcwind_available", "N")) != "Y" or not p or not os.path.isfile(p):
        return
    err = bcwind_t0_error(p, start_date or sim["start_date"], args.get("coupling"))
    if err:
        raise ValueError(err)


def bcwind_t0_error(bcwind_file, start_date, coupling=None):
    """The error if the wind file's T=0 is not the run's start, else None.

    gamera takes UT(T=0) from the wind file; in a coupled run the start is engage's start date
    minus gr_warm_up_time.
    """
    try:
        t0 = bcwind_file_info(bcwind_file)["t0"]
    except (OSError, ValueError, KeyError):
        return None
    start = parse_run_datetime(start_date, "simulation.start_date")
    if t0 is None or parse_run_datetime(t0, "the wind file's T=0") == start:
        return None
    name = os.path.basename(str(bcwind_file))
    start_s = start.strftime(RUN_DATETIME_FORMAT)
    warm = float((coupling or {}).get("gr_warm_up_time") or 0) if coupling else None
    if warm is None:
        return (f"the solar wind file {name} has T=0 at {t0}, not simulation.start_date {start_s} "
                f"(gamera takes UT from the file): start_date {t0} matches the file.")
    match = (parse_run_datetime(t0, "the wind file's T=0")
             + datetime.timedelta(seconds=warm)).strftime(RUN_DATETIME_FORMAT)
    return (f"the solar wind file {name} has T=0 at {t0}, not the GR start {start_s} (engage "
            f"simulation.start_date minus gr_warm_up_time {_secs(warm)} s; gamera takes UT from the "
            f"file): engage start_date {match} matches the file.")


def bcwind_window_error(bcwind_file, start_date, stop_date, name="simulation.bcwind_file"):
    """The error if the wind file is missing or does not cover the window, else None.

    The model holds the last sample for one spacing, so the run may stop at UT[-1] + spacing.
    """
    if not os.path.isfile(bcwind_file):
        return f"bcwind file {bcwind_file!r} ({name}) does not exist."
    info = bcwind_file_info(bcwind_file)
    w0, w1, dt = info["w0"], info["w1"], info["spacing"]
    t0 = parse_run_datetime(start_date, "simulation.start_date")
    t1 = parse_run_datetime(stop_date, "simulation.stop_date")
    last = w1 + datetime.timedelta(seconds=dt)
    if w0 <= t0 and t1 <= last:
        return None
    return (f"The run window {t0.strftime(RUN_DATETIME_FORMAT)} .. "
            f"{t1.strftime(RUN_DATETIME_FORMAT)} is not covered by the solar wind file "
            f"{bcwind_file!r} ({name}), which spans {w0.strftime(RUN_DATETIME_FORMAT)} .. "
            f"{w1.strftime(RUN_DATETIME_FORMAT)} at {_secs(dt)} s spacing (the model holds the last "
            f"sample for one spacing, so the run may stop at {last.strftime(RUN_DATETIME_FORMAT)} "
            f"at the latest). Set simulation.start_date / stop_date inside that range, or use "
            f"bcwind_available=N to fetch the window from CDAWeb.")


def restart_files_in_run_directory(options: dict, run_dir: str = os.curdir) -> list:
    """The run's voltron and gamera restart files in run_dir and its gr_spinup/."""
    runid = glob.escape(restart_runid(options))
    pats = [os.path.join(run_dir, d, f"{runid}{tail}") for d in ("", GR_SPINUP_DIR)
            for tail in (".volt.Res.*.h5", "*.gam.Res.*.h5")]
    return sorted({f for p in pats for f in glob.glob(p) if os.path.exists(f)})


def restart_files_error(options: dict, run_dir: str, clobber: bool = False):
    """The error if a segmented run's directory already holds its restart files (unless clobber)."""
    if clobber or not _is_segmented(options):
        return None
    found = restart_files_in_run_directory(options, run_dir)
    if not found:
        return None
    shown = ", ".join(os.path.relpath(f, run_dir) for f in found[:3]) + (", ..." if len(found) > 3 else "")
    return (f"The run directory {run_dir} already holds restart files of {restart_runid(options)} "
            f"({shown}): use a fresh run directory, or remove them, or rerun with --clobber to "
            f"regenerate the decks (the job chain then starts again from the spin-up).")


def _wind_tsfile(options: dict) -> str:
    return str(options.get("gamera", {}).get("wind", {}).get("tsfile") or CDA2WIND_OUTPUT)


def apply_bcwind_tsfile(options: dict):
    """With bcwind_available=Y, force gamera.wind.tsfile to the bare bcwind.h5."""
    sim = options.get("simulation", {})
    if yn(sim.get("bcwind_available", "N")) != "Y":
        return
    wind = options.setdefault("gamera", {}).setdefault("wind", {})
    if wind.get("tsfile") != CDA2WIND_OUTPUT:
        src = sim.get("bcwind_file")
        why = (f"makeitso links simulation.bcwind_file {src!r} into the run directory as "
               f"{CDA2WIND_OUTPUT}" if src else
               f"without simulation.bcwind_file the run directory's {CDA2WIND_OUTPUT} is used; set "
               f"simulation.bcwind_file to {wind.get('tsfile')!r} to have makeitso link that file in")
        _note(f"gamera.wind.tsfile {wind.get('tsfile')!r} -> {CDA2WIND_OUTPUT!r} "
              f"(bcwind_available=Y: {why}).")
        wind["tsfile"] = CDA2WIND_OUTPUT


def deck_string_error(options: dict):
    """The error for DECK_STRING_OPTIONS values kaiju cannot read (or empty ones), else None."""
    bad = []
    for path in DECK_STRING_OPTIONS:
        sect = options.get(path[0], {}).get(path[1], {})
        v = sect.get(path[2]) if isinstance(sect, dict) else None
        if v is None or (str(v).strip() and not _DECK_STRING_UNREADABLE.search(str(v))):
            continue
        name = ".".join(path)
        hint = ""
        if path == ("gamera", "wind", "tsfile"):
            hint = (" (for a solar wind file elsewhere use bcwind_available=Y with "
                    "simulation.bcwind_file: makeitso links it into the run directory as "
                    f"{CDA2WIND_OUTPUT})")
        elif path == ("gamera", "sim", "runid"):
            hint = " (it follows simulation.job_name)"
        bad.append(f"{name} = {str(v)!r}{hint}")
    if not bad:
        return None
    return (f"{'; '.join(bad)}: kaiju reads {'this deck attribute' if len(bad) == 1 else 'these deck attributes'} "
            f"with a Fortran list-directed read, which stops at '/', ',' or a blank (a value starting "
            f"with '/', or an empty value, keeps the model default), so the model would not see the "
            f"value. Use a non-empty bare name with no '/', ',' or blank.")


def check_deck_strings(options: dict):
    """Raise deck_string_error as ValueError."""
    err = deck_string_error(options)
    if err:
        raise ValueError(err)


def pin_bcwind_file(options: dict, run_dir: str):
    """Make a relative bcwind_file that exists here absolute, since it is re-read from run_dir."""
    sim = options.get("simulation", {})
    if yn(sim.get("bcwind_available", "N")) != "Y":
        return
    if os.path.realpath(run_dir) == os.path.realpath(os.getcwd()):
        return
    p = sim.get("bcwind_file")
    if p and not os.path.isabs(p) and os.path.isfile(p):
        sim["bcwind_file"] = os.path.abspath(p)


def stage_bcwind_file(options: dict, clobber: bool = False):
    """Link (or copy) simulation.bcwind_file into the cwd as bcwind.h5; a different one raises unless clobber."""
    sim = options["simulation"]
    if yn(sim.get("bcwind_available", "N")) != "Y" or not sim.get("bcwind_file"):
        return
    src = os.path.abspath(sim["bcwind_file"])
    dst = CDA2WIND_OUTPUT
    if os.path.lexists(dst):
        if os.path.exists(dst) and (os.path.samefile(src, dst) or filecmp.cmp(src, dst, shallow=False)):
            return
        if not clobber:
            raise FileExistsError(
                f"The run directory {os.getcwd()} already has a {dst} that is not "
                f"simulation.bcwind_file {sim['bcwind_file']!r} ({src}); every deck names {dst}. "
                f"Remove it, or rerun with --clobber to replace it with a link to {src}.")
        os.remove(dst)
    try:
        os.symlink(src, dst)
        how = "linked"
    except OSError:
        shutil.copy2(src, dst)
        how = "copied"
    print(f"makeitso: {how} the solar wind file {src} into the run directory as {dst}.")


def run_preprocessing_steps(options: dict, args: dict):
    """Run genLFM (if needed), cda2wind (bcwind_available=N) and genRAIJU in the run directory."""
    debug = args["debug"]
    sim = options["simulation"]

    grid_type = str(sim["gamera_grid_type"])
    grid_file = f'lfm{grid_type}.h5'
    if existing_grid_file_reusable(options, args.get("clobber", False)):
        print(f'Using existing grid file {grid_file}.')
    else:
        if os.path.exists(grid_file):
            os.remove(grid_file)
        # radii may be JSON numbers: argv needs str
        cmd = 'genLFM'
        cmd_args = [cmd, '-gid', grid_type,
                    '-Rin', str(sim['gamera_grid_inner_radius']),
                    '-Rout', str(sim['gamera_grid_outer_radius'])]
        if debug:
            print(f"cmd = {cmd_args}")
        subprocess.run(cmd_args, check=True)

    if yn(sim["bcwind_available"]) == "N":
        tsfile = _wind_tsfile(options)
        cmd = "cda2wind"
        cmd_args = [cmd, "-t0", str(sim["start_date"]), "-t1", str(sim["stop_date"]),
                    "-interp", "-bx", "-f107", "100", "-kp", "3"]
        if tsfile != CDA2WIND_OUTPUT:
            cmd_args += ["-o", tsfile]
        if debug:
            print(f"cmd = {cmd_args}")
        subprocess.run(cmd_args, check=True)

    cmd = "genRAIJU"
    cmd_args = [cmd]
    subprocess.run(cmd_args, check=True)


def _import_segment_solver():
    """Import waxrun's segment_solver from the kaiju checkout."""
    import sys
    kaijuhome = os.environ.get("KAIJUHOME") or os.path.join(SUPPORT_FILES_DIRECTORY, os.pardir, os.pardir)
    wax_dir = os.path.normpath(os.path.join(kaijuhome, "scripts", "waxrun"))
    if wax_dir not in sys.path:
        sys.path.append(wax_dir)
    import segment_solver
    return segment_solver


def _import_binident():
    """Load engage's binident.py by path, reusing the module engage already imported."""
    import sys
    import importlib.util
    kaijuhome = os.environ.get("KAIJUHOME") or os.path.join(SUPPORT_FILES_DIRECTORY, os.pardir, os.pardir)
    path = os.path.normpath(os.path.join(kaijuhome, "scripts", "engage", "binident.py"))
    loaded = sys.modules.get("binident")
    if loaded is not None and os.path.abspath(getattr(loaded, "__file__", "") or "") == path:
        return loaded
    spec = importlib.util.spec_from_file_location("makeitso_binident", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_voltron_exe(options: dict, args: dict = None):
    """Standalone GR: refuse a non-ELF voltron_mpi.x and warn about missing restart files."""
    if run_component(args) != "gr":
        return
    bi = _import_binident()
    exe = os.path.join(str(options.get("pbs", {}).get("kaiju_build_directory") or ""), "bin", "voltron_mpi.x")
    if os.path.isfile(exe) and not bi.is_elf(exe):
        raise ValueError(f"{exe} is not an ELF executable: the job scripts copy it in as voltron_mpi.x. "
                         f"Point pbs.kaiju_build_directory at a kaiju build whose bin/voltron_mpi.x is the "
                         f"compiled voltron.")
    gamera = options.get("gamera", {})
    if _is_segmented(options) or tf(gamera.get("restart", {}).get("doRes")) != "T":
        return
    runid = restart_runid(options)
    prefix = os.path.join(os.getcwd(), runid)
    found = sorted(f for f in glob.glob(glob.escape(prefix) + ".volt.Res.[0-9][0-9][0-9][0-9][0-9].h5")
                   if os.path.exists(f))
    if not found:
        print(f"{YELLOW}WARNING: the deck restarts from {runid}.volt.Res.*.h5, which the run "
              f"directory {os.getcwd()} does not hold yet: stage the restart files before "
              f"submitting.{RESET}")
        return
    try:
        index = int(gamera["restart"].get("nRes"))
    except (TypeError, ValueError, KeyError):
        index = -1
    if index < 0:
        index = int(found[-1].rsplit(".", 2)[-2])
    missing = [f for f in bi.RES_FAMILIES if f not in bi.RES_OPTIONAL and not bi.res_files(prefix, index, f)]
    if missing:
        names = [f"{runid}*.gam.Res.{index:05d}.h5" if f == "gam" else f"{runid}.{f}.Res.{index:05d}.h5"
                 for f in missing]
        print(f"{YELLOW}WARNING: the deck restarts from index {index:05d}, but the run directory "
              f"has no {', '.join(names)} (voltron stops with 'Restart file not found'): stage "
              f"them before submitting.{RESET}")
    for f in bi.RES_OPTIONAL:
        if not bi.res_files(prefix, index, f):
            print(f"{YELLOW}WARNING: no {runid}.{f}.Res.{index:05d}.h5: DragonKing cold-starts from the "
                  f"restart.{RESET}")


def dtcouple_prompt_check(value):
    """Prompt check of dtCouple in a coupled run (a whole number of seconds dividing 60)."""
    return _import_segment_solver().dtcouple_error(value)


def warn_dtcouple(options: dict):
    """Warn once per run about a dtCouple other than 5 s (gamera's inner boundary assumes 5 s)."""
    dtc = options.get("voltron", {}).get("coupling", {}).get("dtCouple")
    _import_segment_solver().warn_dtcouple_once(dtc, emit=lambda s: print(f"{YELLOW}{s}{RESET}"))


def run_component(args) -> str:
    """'grx', 'grt' (engage coupled runs) or 'gr' (standalone)."""
    args = args or {}
    if args.get("grx"):
        return "grx"
    return "grt" if "coupling" in args else "gr"


def _is_segmented(options: dict) -> bool:
    return yn(options.get("simulation", {}).get("use_segments", "N")) == "Y"


def grid_type_error(grid, component: str = "gr"):
    """The error if `grid` is not allowed for the run kind, else None."""
    grids = RUN_GRID_TYPES[component]
    if str(grid) in grids:
        return None
    if component == "grx":
        return (f"simulation.gamera_grid_type={grid!r} is refused for a GRX run: GRX is validated on "
                f"grid Q only (the only grid a GRX run has used). Set gamera_grid_type='Q'.")
    if component == "grt" and str(grid) == "H":
        return ("simulation.gamera_grid_type='H' is refused for a coupled GRT run: GRT is run and "
                "tested on grids D, Q and O only. H is a developer grid that voltron refuses without "
                "a developer override the tools do not write. Use D, Q or O (H stays available to "
                "a standalone GR run).")
    return (f"simulation.gamera_grid_type={grid!r} is not a GAMERA grid type for this run "
            f"({'|'.join(grids)}).")


NOSPIN_REMEDY = ("set voltron.spinup.doSpin='T' (recommended) or restart from an existing run "
                 "(gamera.restart.doRes='T')")


def coupled_nospin_error(options: dict, component: str):
    """The error for a GRT run without the spin-up (WARMUP-01 restarts from its index 0), else None."""
    if component != "grt":
        return None
    value = options.get("voltron", {}).get("spinup", {}).get("doSpin", "T")
    if tf(value) == "T":
        return None
    return (f"voltron.spinup.doSpin={value!r} is refused for a coupled GRT run: the GR warm-up needs "
            f"the spin-up restart (WARMUP-01 restarts from index 0, written only by the SPINUP "
            f"segment); set voltron.spinup.doSpin='T'.")


def coupled_dospin_prompt_check(value):
    return coupled_nospin_error({"voltron": {"spinup": {"doSpin": value}}}, "grt")


def _grid_rule(options: dict, component: str):
    grid = options.get("simulation", {}).get("gamera_grid_type")
    return grid_type_error(grid, component) if grid not in (None, "") else None


def standalone_dogcm_error(options: dict, component: str):
    """The error for doGCM='T' in a standalone run, else None."""
    if component != "gr":
        return None
    if tf(options.get("voltron", {}).get("coupling", {}).get("doGCM", "F")) != "T":
        return None
    return ("voltron.coupling.doGCM='T' is refused for a standalone run: makeitso launches no GCM. "
            "Use engage (--component grt or grx).")


# Each rule(options, component) returns an error message or None.
OPTION_RULES = (_grid_rule, coupled_nospin_error, standalone_dogcm_error)


def option_rules_error(options: dict, args: dict = None):
    """The first error from OPTION_RULES, else None."""
    component = run_component(args)
    for rule in OPTION_RULES:
        err = rule(options, component)
        if err:
            return err
    return None


# doGCM of the MAGE-only lead-in decks; coupled decks set their own.
COUPLED_OPTIONS_DOGCM = "F"


def apply_option_rules(options: dict, args: dict = None):
    """Raise option_rules_error as ValueError."""
    err = option_rules_error(options, args)
    if err:
        raise ValueError(err)
    return options


def grx_segment_plan(options: dict, stephours, delaycouplinghrs, t0, t1):
    """Per-segment GRX deck overrides [{tFin, nRes, doRes, doGCM[, doSpin]}] (all str).

    Raises ValueError for a cadence that breaks the coupled timing rules. Writes nothing.
    """
    seg = _import_segment_solver()
    step, delay = int(stephours), int(delaycouplinghrs)
    cpl = options["voltron"]["coupling"]
    coupling_start = seg._parse_dt(t0) + datetime.timedelta(hours=delay)
    err = (seg.cadence_error(cpl.get("dtCouple"), tSpin=options["voltron"]["spinup"].get("tSpin"),
                             gr_warm_up_time=delay * 3600, segment_duration=step * 3600,
                             coupled_start=coupling_start)
           or seg.coupling_start_error(coupling_start, cpl.get("dtCouple"),
                                       name="the GRX coupling start (simulation.start_date)"))
    if err:
        raise ValueError(err)
    dtCouple = float(cpl["dtCouple"])
    ic, iu, maxhours = seg.compute_segments(t0, t1, step, delay)
    dtRes = float(options["voltron"]["restart"]["dtRes"])
    plan = []
    for ii in range(0, iu):
        nRes = seg.nres(ii, ic, delay, step, dtRes=dtRes)
        overrides = {"tFin": str(seg.tfin(ii, ic, delay, step, dtCouple=dtCouple)), "nRes": str(nRes),
                     "doRes": "T" if nRes >= 0 else "F", "doGCM": "T" if ii > ic else "F"}
        if ii == 0:
            overrides["doSpin"] = "T"
        plan.append(overrides)
    return plan


def grx_is_lead_in(overrides: dict) -> bool:
    """True for a GRX lead-in deck (doGCM=F), which lives in gr_spinup/."""
    return str(overrides.get("doGCM", "T")).upper() == "F"


def gen_grx_voltron(options: dict, stephours, delaycouplinghrs, t0, t1, run_dir=".", args=None):
    """Write the per-segment GRX voltron decks and return the .xml paths.

    Lead-in decks go to <run_dir>/gr_spinup/, coupled decks to run_dir.
    """
    plan = grx_segment_plan(options, stephours, delaycouplinghrs, t0, t1)
    iu = len(plan)
    seg = _import_segment_solver()
    ic = seg.compute_segments(t0, t1, int(stephours), int(delaycouplinghrs))[0]
    job_name = options["simulation"]["job_name"]
    check_deck_strings(options)

    # StrictUndefined: a missing option must fail, not render as an empty attribute
    with open(INI_TEMPLATE, "r", encoding="utf-8") as f:
        template = Template(f.read(), undefined=StrictUndefined)

    # render all before writing any
    texts = []
    for ii, overrides in enumerate(plan):
        opt = copy.deepcopy(options)
        opt["voltron"]["time"]["tFin"] = overrides["tFin"]
        opt["gamera"]["restart"]["nRes"] = overrides["nRes"]
        opt["gamera"]["restart"]["doRes"] = overrides["doRes"]
        opt["voltron"]["coupling"]["doGCM"] = overrides["doGCM"]
        if "doSpin" in overrides:
            opt["voltron"]["spinup"]["doSpin"] = overrides["doSpin"]
        texts.append(template.render(opt, grx_deck=GRX_DECK_KEYS))

    os.makedirs(run_dir, exist_ok=True)
    ini_files = []
    for ii, text in enumerate(texts):
        where = lead_in_directory(run_dir) if grx_is_lead_in(plan[ii]) else run_dir
        ini_file = os.path.join(where, f"{seg.segment_id(job_name, ii, ic)}.ini")
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(text)
        ini_files.append(ini_file)
    volt_files = convert_ini_to_xml(ini_files, args if args is not None else {"debug": False})
    lead = sum(1 for p in plan if grx_is_lead_in(p))
    print(f"makeitso: wrote {iu} GRX voltron config files to {run_dir}"
          + (f" ({lead} GR lead-in decks in {GR_SPINUP_DIR}/)" if lead else ""))
    return volt_files


def prepare_run_directory(run_directory: str) -> str:
    """Create pbs.run_directory if needed and return its absolute path."""
    run_dir = os.path.abspath(os.path.expanduser(str(run_directory)))
    if os.path.exists(run_dir) and not os.path.isdir(run_dir):
        raise NotADirectoryError(
            f"pbs.run_directory {run_directory!r} ({run_dir}) exists but is not a directory.")
    try:
        os.makedirs(run_dir, exist_ok=True)
    except OSError as e:
        raise OSError(
            f"cannot create pbs.run_directory {run_directory!r} ({run_dir}): {e}") from e
    return run_dir


def run_directory_path(options: dict, filename: str) -> str:
    """Path of `filename` in the run directory, spelled from inside it (the cwd when writing)."""
    rd = str(options["pbs"]["run_directory"])
    return os.path.join(rd if os.path.isabs(rd) else os.curdir, filename)


# Decks carry bare file names, so gr_spinup/ gets relative symlinks to the run directory's inputs,
# and voltron reads restarts from its cwd, so the first main segment links them in.
_LEAD_IN_SEGMENT = re.compile(r"-(SPINUP|WARMUP-\d+)$")


def is_lead_in(segment_id) -> bool:
    """True for a GR lead-in deck (<job>-SPINUP or <job>-WARMUP-NN)."""
    return bool(_LEAD_IN_SEGMENT.search(str(segment_id)))


def deck_dir(segment_id) -> str:
    """The deck's directory relative to the run directory: GR_SPINUP_DIR for lead-in decks, else '.'."""
    return GR_SPINUP_DIR if is_lead_in(segment_id) else os.curdir


def deck_path(options: dict, segment_id, ext: str) -> str:
    """Path of a deck's <segment_id><ext> file in its deck_dir, spelled like run_directory_path."""
    sub = deck_dir(segment_id)
    name = f"{segment_id}{ext}"
    return run_directory_path(options, name if sub == os.curdir else os.path.join(sub, name))


def lead_in_directory(run_dir: str = os.curdir) -> str:
    """Create <run_dir>/gr_spinup and return its path."""
    path = os.path.join(run_dir, GR_SPINUP_DIR)
    os.makedirs(path, exist_ok=True)
    return path


def lead_in_input_names(options: dict) -> list:
    """The bare input file names a lead-in deck reads (grid, solar wind, RAIJU config)."""
    names = []
    for path in (("gamera", "sim", "H5Grid"), ("gamera", "wind", "tsfile"), ("raiju", "config", "fname")):
        v = options.get(path[0], {}).get(path[1], {}).get(path[2])
        if v not in (None, "") and str(v) not in names:
            names.append(str(v))
    return names


def lead_in_input_links(options: dict, clobber: bool = False, run_dir: str = os.curdir) -> list:
    """Link gr_spinup/<name> -> ../<name> for each lead-in input and return the link paths.

    A different real file already there raises FileExistsError unless clobber.
    """
    sub = lead_in_directory(run_dir)
    links = []
    for name in lead_in_input_names(options):
        link = os.path.join(sub, name)
        target = os.path.join(os.pardir, name)
        if os.path.islink(link):
            if os.readlink(link) == target:
                links.append(link)
                continue
            os.remove(link)
        elif os.path.exists(link):
            top = os.path.join(run_dir, name)
            if os.path.isfile(top) and filecmp.cmp(link, top, shallow=False):
                links.append(link)
                continue
            if not clobber:
                raise FileExistsError(
                    f"{os.path.abspath(link)} is a file, not the link {GR_SPINUP_DIR}/{name} -> "
                    f"../{name} that makeitso makes (the GR lead-in decks read the run directory's "
                    f"{name} through it). Remove it, or rerun with --clobber to replace it.")
            os.remove(link)
        os.symlink(target, link)
        links.append(link)
    return links


def handoff_lines(runid, index: int) -> list:
    """Bash lines that hard-link restart set `index` from gr_spinup/ into the job directory.

    A missing required restart exits 90 before the launch.
    """
    bi = _import_binident()
    raw, runid = str(runid), shlex.quote(str(runid))
    src = f"{GR_SPINUP_DIR}/"
    required = [f for f in bi.RES_FAMILIES if f not in bi.RES_OPTIONAL and f != "gam"]
    optional = [f for f in bi.RES_FAMILIES if f in bi.RES_OPTIONAL]
    return (
        [f"# ---- GR restart hand-off: link restart index n from {src}",
         f"n={int(index):05d}",
         f'pull() {{ s="{src}$1"; d="./$1"',
         '  if [ -e "$s" ]; then [ "$s" -ef "$d" ] || ln -f "$s" "$d" 2>/dev/null || ln -sf "$s" "$d"',
         '  elif [ -e "$d" ]; then :',
         '  else return 1; fi; }',
         'miss() { echo "GATE: $1"; exit 90; }',
         f'gam=$(for s in {src}{runid}*.gam.Res.$n.h5 ./{runid}*.gam.Res.$n.h5; do '
         'if [ -e "$s" ]; then echo "${s##*/}"; fi; done | sort -u)',
         f'[ -n "$gam" ] || miss "no {raw}*.gam.Res.$n.h5 in {src} nor here"',
         "for f in $gam " + " ".join(f"{runid}.{fam}.Res.$n.h5" for fam in required) + "; do",
         f'  pull "$f" || miss "restart $f is neither in {src} nor here"',
         "done"]
        + [f'pull {runid}.{fam}.Res.$n.h5 || echo "WARNING: {raw}.{fam}.Res.$n.h5 is neither in '
           f'{src} nor here: DragonKing cold-starts"' for fam in optional]
        + ["# ---- end GR restart hand-off"])


def json_relative_path(value, json_path: str, name: str = "pbs.run_directory"):
    """Resolve a relative path from an options JSON against the JSON's own directory."""
    if value in (None, "") or os.path.isabs(os.path.expanduser(str(value))):
        return value
    value = os.path.expanduser(str(value))
    jdir = os.path.dirname(os.path.abspath(json_path))
    by_json = os.path.normpath(os.path.join(jdir, value))
    by_cwd = os.path.abspath(value)
    if by_cwd == by_json:
        return value
    _note(f"{name} {value!r} is relative to the directory of {os.path.basename(json_path)}, not to the "
          f"current directory.")
    return by_json


def resolve_json_run_directory(options: dict, json_path: str) -> dict:
    pbs = options.get("pbs")
    if isinstance(pbs, dict) and pbs.get("run_directory") not in (None, ""):
        pbs["run_directory"] = json_relative_path(pbs["run_directory"], json_path)
    return options


def run_options_for_json(options: dict, json_dir: str = None, run_dir: str = None) -> dict:
    """A copy of options with pbs.run_directory spelled relative to json_dir when it is below it."""
    out = copy.deepcopy(options)
    pbs = out.get("pbs")
    if not isinstance(pbs, dict) or pbs.get("run_directory") in (None, ""):
        return out
    json_dir = os.path.abspath(json_dir or os.getcwd())
    rd = os.path.abspath(os.path.expanduser(str(run_dir if run_dir is not None else pbs["run_directory"])))
    pbs["run_directory"] = JobGen.run_relative(rd, json_dir, json_dir, dot=False)
    return out


def saved_json_spelling(run_dir: str, launch_dir: str) -> str:
    """PARAMETERS_JSON in run_dir, spelled relative to launch_dir when below it, else absolute."""
    path = os.path.join(os.path.abspath(run_dir), PARAMETERS_JSON)
    rel = os.path.relpath(path, os.path.abspath(launch_dir))
    return path if rel.split(os.sep)[0] == os.pardir else rel


def dump_run_options(path: str, options: dict, run_dir: str = None):
    """Write the options JSON with pbs.run_directory relative to the JSON's directory."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(run_options_for_json(options, os.path.dirname(os.path.abspath(path)),
                                       run_dir=run_dir or os.getcwd()), f, indent=JSON_INDENT)


def submit_entries(scripts: list, base: str = None) -> list:
    """(subdir, script) entries of a submit-all chain, subdir relative to base when below it."""
    base = os.path.abspath(base or os.getcwd())
    entries = []
    for p in scripts:
        ap = os.path.abspath(os.fspath(p))
        rel = os.path.relpath(ap, base)
        sub, name = os.path.split(rel)
        if rel.startswith(os.pardir + os.sep):
            sub = os.path.dirname(ap)
        entries.append((sub or os.curdir, name))
    return entries


def handoff_plan(options: dict, args: dict, decks: list = None) -> dict:
    """{segment_id: restart index} of the first deck after the GR lead-in, or {}."""
    if decks is None:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                decks = plan_segments(copy.deepcopy(options), args or {})
        except ValueError:
            return {}
    lead = [k for k, d in enumerate(decks) if is_lead_in(d["segment_id"])]
    if not lead or lead[-1] + 1 >= len(decks):
        return {}
    nxt = decks[lead[-1] + 1]
    if nxt.get("doRes") == "T" and nxt.get("nRes") is not None:
        return {nxt["segment_id"]: int(nxt["nRes"])}
    return {}


def _secs(value):
    """Format seconds for a message (7200 for 7200.0)."""
    return f"{value:.0f}" if float(value).is_integer() else f"{value}"


def _year_split(start_dt, t0, t1):
    """Offset (s from start_dt) of the Jan 1 at which window (t0, t1] must be cut, or None.

    Must match tiegcmrun's segment_time rule so both models split at the same instant.
    """
    dt0 = start_dt + datetime.timedelta(seconds=t0)
    dt1 = start_dt + datetime.timedelta(seconds=t1)
    if dt0.year == dt1.year:
        return None
    if (dt1.month, dt1.day, dt1.hour, dt1.minute, dt1.second) == (1, 1, 0, 0, 0):
        return None
    return (datetime.datetime(dt0.year + 1, 1, 1) - start_dt).total_seconds()


def _align_restart_cadence(options: dict, starts: list,
                           inputs: str = "segment_duration and gr_warm_up_time") -> float:
    """Return a dtRes (s) that divides every segment start, lowering it in options if needed.

    Raises ValueError if that needs more than a MAX_DTRES_REDUCTION-fold reduction.
    """
    dtRes = float(options["voltron"]["restart"]["dtRes"])
    if dtRes <= 0:
        raise ValueError(f"voltron.restart.dtRes must be a positive number of seconds, got "
                         f"{options['voltron']['restart']['dtRes']!r}.")
    off = [t for t in starts if t % dtRes]
    if not off:
        return dtRes
    values = [dtRes] + list(starts)
    cadence = 0
    if all(float(v).is_integer() for v in values):
        for v in values:
            cadence = math.gcd(cadence, int(v))

    shown = ", ".join(_secs(t) for t in off[:4]) + (", ..." if len(off) > 4 else "")
    if cadence * MAX_DTRES_REDUCTION < dtRes:
        raise ValueError(
            f"Every segment must start on a voltron restart dump, i.e. at a multiple of "
            f"voltron.restart.dtRes ({_secs(dtRes)} s), but segments would start at {shown} s "
            f"after the GR start. Choose a {inputs} that dtRes "
            f"divides (multiples of {_secs(dtRes)} s), or set voltron.restart.dtRes "
            f"(INTERMEDIATE) to a cadence that divides them.")
    print(f"{YELLOW}NOTE: voltron.restart.dtRes lowered from {_secs(dtRes)} s to {cadence} s so "
          f"that every segment starts on a restart dump (segments start at {shown} s after the "
          f"GR start, which {_secs(dtRes)} s does not divide).{RESET}")
    options["voltron"]["restart"]["dtRes"] = str(float(cadence))
    return float(cadence)


def restart_runid(options: dict) -> str:
    return str(options.get("gamera", {}).get("sim", {}).get("runid") or options["simulation"]["job_name"])


def plan_segments(options: dict, args: dict):
    """Plan the decks of a run in order: [{segment_id, tFin, doRes, nRes, doGCM, year_split}].

    Each segment restarts from the dump at its own start (nRes = start/dtRes). Also validates
    the run and sets derived values (num_segments, dtRes) in options.
    """
    runid = options["simulation"]["job_name"]
    coupled = "coupling" in args
    err = option_rules_error(options, args)
    if err:
        raise ValueError(err)
    do_spin = tf(options["voltron"]["spinup"]["doSpin"]) == "T"
    nospin = (
        "voltron.spinup.doSpin='F' cannot start: voltron stops with 'Spinup is required with deep "
        "coupling' unless the run starts from a restart.")

    if yn(options["simulation"]["use_segments"]) != "Y":
        if not coupled and not do_spin and tf(options["gamera"]["restart"].get("doRes")) != "T":
            raise ValueError(f"{nospin} To run it: {NOSPIN_REMEDY}.")
        options["pbs"]["num_segments"] = "1"
        return [{"segment_id": runid}]

    tFin = float(options["voltron"]["time"]["tFin"])
    segment_duration = float(options["simulation"]["segment_duration"])
    decks = []
    if do_spin:
        # tFin 1.0 so the restart dump at t=0 is written
        decks.append({"segment_id": f"{runid}-SPINUP", "tFin": 1.0, "doRes": "F"})
        if coupled:
            decks[-1]["doGCM"] = COUPLED_OPTIONS_DOGCM

    if not coupled:
        if not do_spin:
            raise ValueError(f"{nospin} Segment 01 of a segmented run without the SPINUP segment "
                             f"has no restart to start from: {NOSPIN_REMEDY} in an unsegmented "
                             f"deck (use_segments='N').")
        num_segments = math.ceil(tFin/segment_duration)
        given = options["pbs"].get("num_segments")
        if str(given) != str(num_segments):
            try:
                same = int(given) == num_segments
            except (TypeError, ValueError):
                same = False
            if not same and given is not None:
                print(f"{YELLOW}WARNING: pbs.num_segments={given!r} ignored: the run needs "
                      f"{num_segments} segment(s) of segment_duration={_secs(segment_duration)} "
                      f"s to reach tFin={_secs(tFin)} s.{RESET}")
            options["pbs"]["num_segments"] = str(num_segments)
        starts = [(job - 1)*segment_duration for job in range(1, num_segments + 1)]
        dtRes = _align_restart_cadence(options, starts, inputs="segment_duration")
        nres = [int(t/dtRes) for t in starts]
        for job in range(1, num_segments + 1):
            deck = {"segment_id": f"{runid}-{job:02d}",
                    # +1 s so the end dump is written; the last stops at the requested end
                    "tFin": job*segment_duration + 1.0 if job < num_segments else tFin,
                    "doRes": "T", "nRes": nres[job - 1]}
            decks.append(deck)
        return decks

    coupling = args["coupling"]
    gr_warm_up_time = float(coupling["gr_warm_up_time"])
    start_dt = datetime.datetime.fromisoformat(options["simulation"]["start_date"])

    volt = options["voltron"]
    err = _import_segment_solver().cadence_error(
        volt["coupling"].get("dtCouple"),
        tSpin=volt["spinup"].get("tSpin") if do_spin else None,
        gr_warm_up_time=gr_warm_up_time, segment_duration=segment_duration,
        coupled_start=start_dt + datetime.timedelta(seconds=gr_warm_up_time),
        dtOut=volt.get("output", {}).get("dtOut"), window=tFin - gr_warm_up_time)
    if err:
        raise ValueError(err)

    if gr_warm_up_time % segment_duration == 0:
        piece, num_warmup = segment_duration, int(gr_warm_up_time/segment_duration)
    elif gr_warm_up_time < segment_duration:
        piece, num_warmup = gr_warm_up_time, 1
    else:
        piece, num_warmup = segment_duration, math.ceil(gr_warm_up_time/segment_duration)
    warmup = [((job - 1)*piece, min(job*piece, gr_warm_up_time))
              for job in range(1, num_warmup + 1)]

    windows = []
    year_crossings = 0
    for job in range(1, math.ceil((tFin - gr_warm_up_time)/segment_duration) + 1):
        t0 = (job - 1)*segment_duration + gr_warm_up_time
        t1 = min(job*segment_duration + gr_warm_up_time, tFin)
        cut = _year_split(start_dt, t0, t1)
        if cut is None:
            windows.append((t0, t1, False))
        else:
            year_crossings += 1
            windows += [(t0, cut, True), (cut, t1, False)]
    # only one year crossing has been validated end to end
    if year_crossings >= 2:
        raise NotImplementedError(
            f"This coupled run crosses {year_crossings} calendar-year boundaries; the "
            f"segmenter only supports one. Split it into per-year pieces.")

    starts = [t0 for t0, _ in warmup] + [w[0] for w in windows]
    dtRes = _align_restart_cadence(options, starts)
    nres = [int(t/dtRes) for t in starts]

    for job, (t0, t1) in enumerate(warmup, start=1):
        decks.append({"segment_id": f"{runid}-WARMUP-{job:02d}", "tFin": t1 + 1.0,
                      "doRes": "T", "nRes": nres[job - 1], "doGCM": COUPLED_OPTIONS_DOGCM})

    # tfin_delta=T: all but the last coupled segment run dtCouple-1 s past their end, so the
    # exchange at the boundary happens once, inside the earlier segment.
    if tf(coupling["tfin_delta"]) == "T":
        tfin_coupling_padding = float(options["voltron"]["coupling"]["dtCouple"]) - 1
    else:
        tfin_coupling_padding = 0.0
    for job, (t0, t1, year_split) in enumerate(windows, start=1):
        tfin_padding = -1.0 if job == len(windows) else tfin_coupling_padding - 1.0
        decks.append({"segment_id": f"{runid}-{job:02d}", "tFin": (t1 + 1.0) + tfin_padding,
                      "doRes": "T", "nRes": nres[len(warmup) + job - 1],
                      "doGCM": "T", "year_split": year_split})
    return decks


def create_ini_files(options: dict, args: dict):
    """Render template.ini once per planned deck and return the .ini paths."""
    with open(INI_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    # StrictUndefined: a missing option must fail, not render as an empty attribute
    template = Template(template_content, undefined=StrictUndefined)

    ini_files = []
    for deck in plan_segments(options, args):
        opt = copy.deepcopy(options)
        opt["simulation"]["segment_id"] = deck["segment_id"]
        if "tFin" in deck:
            opt["voltron"]["time"]["tFin"] = str(deck["tFin"])
        if "doRes" in deck:
            opt["gamera"]["restart"]["doRes"] = deck["doRes"]
        if "nRes" in deck:
            opt["gamera"]["restart"]["nRes"] = str(deck["nRes"])
        if "doGCM" in deck:
            opt["voltron"]["coupling"]["doGCM"] = deck["doGCM"]
        if deck.get("year_split"):
            print(f"Segment {deck['segment_id']} crosses year boundary, adjusting segment end "
                  f"time to Jan 1 of next year.")
        try:
            ini_content = template.render(opt)
        except UndefinedError as e:
            raise ValueError(f"The run options lack a value template.ini needs ({e.message}); "
                             f"add it to the options JSON.") from e
        ini_file = deck_path(opt, deck["segment_id"], ".ini")
        if deck_dir(deck["segment_id"]) != os.curdir:
            lead_in_directory(os.path.dirname(os.path.dirname(ini_file)) or os.curdir)
        ini_files.append(ini_file)
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)

    return ini_files


def _wait_for_removal(path: str, timeout: float = 10.0):
    deadline = time.monotonic() + timeout
    while os.path.lexists(path) and time.monotonic() < deadline:
        time.sleep(0.05)


def convert_ini_to_xml(ini_files: list, args: dict):
    """Convert each .ini to .xml beside it, delete the .ini, and return the XML paths."""
    debug = args["debug"]

    xml_files = []
    for ini_file in ini_files:

        xml_file = os.path.splitext(ini_file)[0] + ".xml"

        # Each call gets its own scratch dir: XMLGenerator ends with an unwaited
        # `rm -rf .Settings/` that races the next call in a shared directory.
        cmd = "XMLGenerator"
        cmd_args = [cmd, os.path.abspath(ini_file), os.path.abspath(xml_file)]
        if debug:
            print(f"cmd = {cmd_args}")
        scratch = tempfile.mkdtemp(prefix=".xmlgen-",
                                   dir=os.path.dirname(os.path.abspath(xml_file)))
        try:
            subprocess.run(cmd_args, check=True, cwd=scratch)
            _wait_for_removal(os.path.join(scratch, ".Settings"))
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

        xml_files.append(xml_file)

        os.remove(ini_file)

    return xml_files


def render_pbs(options: dict, stack: str = "gr", handoff=None, job_dir: str = None) -> str:
    """Render a MAGE job script for options["simulation"]["segment_id"] via jobgen.

    stack selects the module set where machines.yaml has several; job_dir, if given, is cd'ed to
    by absolute path instead of the scheduler's submit directory.
    """
    sim = options["simulation"]
    pbs = options["pbs"]
    jg = JobGen(CONFIG_DIR, machine=sim["hpc_system"],
                node_type=sim.get("node_type"))

    # two task groups: GAMERA compute, then voltron + helpers
    compute = {"mpi_tasks": int(pbs["select"]) * int(pbs["mpiprocs"]),
               "tasks_per_node": int(pbs["mpiprocs"]),
               "threads": int(pbs["ompthreads"])}
    helper = {"mpi_tasks": int(pbs["select2"]) * int(pbs["helper_mpiprocs"]),
              "tasks_per_node": int(pbs["helper_mpiprocs"]),
              "threads": int(pbs["helper_ompthreads"])}
    task_groups = [compute, helper]

    segment_id = sim["segment_id"]
    jobid = jg.sched("jobid")
    jobname = jg.sched("jobname", job_name=segment_id)
    build_dir = pbs["kaiju_build_directory"]

    L = [jg.shebang(), "",
         f"# This script was generated to run on {sim['hpc_system']}.", ""]
    # SLURM: job-wide directives must precede the heterogeneous resource block
    pre = jg.directive_lines({
        "job_name": segment_id,
        "account": pbs.get("account_name"),
        "group_list": pbs.get("group_list"),
        "queue": pbs.get("queue"),
    })
    res = jg.resource_lines(task_groups, pbs["walltime"], queue=pbs.get("queue")) + constraint_lines(jg)
    post = jg.directive_lines({
        "priority": pbs.get("job_priority"),
        "mail": "abe",
        "join_out": True,
    })
    L += pre + post + res

    L += ["",
          f'echo "Job {jobid} started at `date` on `hostname`."', "",
          "# Specify the ID string for the run.",
          f"export RUNID={jobname}", "",
          "# Run in the job directory.",
          (f"cd {shlex.quote(os.path.abspath(job_dir))}" if job_dir else jg.sched("workdir_cd")), "",
          "# Load the required modules for kaiju."]
    L += stack_module_lines(jg, stack)
    L += ["module list"]
    L += [m for m in as_list(pbs.get("local_modules")) if m not in (None, "None")]

    L += ["",
          "# Define the kaiju installation location.",
          f"export KAIJU_INSTALL_DIR={shlex.quote(str(pbs['kaiju_install_directory']))}", "",
          "# Set KAIJUHOME and the other kaiju variables.",
          "source $KAIJU_INSTALL_DIR/scripts/setupEnvironment.sh", "",
          "# Add the kaiju build bin directory to PATH.",
          f"export PATH={shlex.quote(str(build_dir) + '/bin')}:$PATH", "",
          "# Set other environment variables.",
          "export MPI_TYPE_DEPTH=32",
          "export OMP_STACKSIZE=100M",
          f"export OMP_NUM_THREADS={pbs['helper_ompthreads']}", "",
          'echo "The active environment variables are:"',
          "printenv"]

    nodefile_src = jg.sched("nodefile")
    if nodefile_src:
        L += ["",
              "# Copy the node file for use by the placer command.",
              f"nodefile=nodefile.{jobid}",
              f"cp {nodefile_src} $nodefile"]

    L += ["",
          "# Copy the main executable.",
          f"[ -L ./voltron_mpi.x ] || cp {shlex.quote(str(build_dir) + '/bin/voltron_mpi.x')} ./voltron_mpi.x", ""]
    if handoff:
        L += list(handoff) + [""]
    L += ["# Run the model. Direct output from the program is saved in a text file.",
          "EXE=./voltron_mpi.x",
          f'OUTPUT_FILE="$EXE-{segment_id}.out"',
          'echo "Running $EXE on model $RUNID."',
          f'{launch_prefix(jg, task_groups)} $EXE "$RUNID.xml" &>> $OUTPUT_FILE',
          # exit with the launch's status: the next segment depends afterok on it
          "rc=$?", "",
          f'echo "Job {jobid} ended at `date` on `hostname`."',
          "exit $rc"]
    return "\n".join(L) + "\n"


def stack_module_lines(jg, stack=None):
    """Module purge/load lines of the machine's module set for a run type (stack)."""
    stacks = jg.m.get("stacks") or {}
    stack = stack if stack in stacks else None
    if stacks and stack is None:
        stack = "gr" if "gr" in stacks else next(iter(stacks))
    lines = jg.module_lines(stack=stack) if stacks else jg.module_lines()
    purge = jg.stack_option(stack, "module_purge") if stacks else None
    if purge and lines and lines[0] == jg.m.get("module_purge"):
        lines[0] = purge
    return lines


def launch_prefix(jg, task_groups):
    """The machine's launcher prefix with {total_tasks} filled in."""
    return jg.m["mpirun"].format(total_tasks=jg.total_tasks(task_groups))


def constraint_lines(jg):
    """The machine's node constraint directive (e.g. #SBATCH --constraint=cpu), or []."""
    constraint = jg.m.get("constraint")
    tmpl = (jg.b.get("resource") or {}).get("node_selector_directive")
    if not constraint or not tmpl:
        return []
    return [f"{jg.b['directive_prefix']} {tmpl.format(v=constraint)}"]


def create_pbs_scripts(xml_files: list, options: dict, args: dict):
    """Write one job script per deck and the submit-all script <job>_pbs.sh.

    Returns (pbs_scripts, submit_all_jobs_script, spinup_pbs_scripts, warmup_pbs_scripts).
    """
    debug = args["debug"]

    derive_select(options)

    handoff = handoff_plan(options, args)

    pbs_scripts = []
    for xml_file in xml_files:
        opt = copy.deepcopy(options)
        runid = opt["simulation"]["job_name"]
        filename = os.path.split(xml_file)[-1]
        if debug:
            print(f"filename = {filename}")
        segment_id = filename.replace(".xml", "")
        if debug:
            print(f"segment_id = {segment_id}")
        opt["simulation"]["segment_id"] = segment_id
        lines = None
        if segment_id in handoff:
            lines = handoff_lines(restart_runid(opt), handoff[segment_id])
        pbs_script = deck_path(opt, segment_id, ".pbs")
        job_dir = None
        if deck_dir(segment_id) != os.curdir:
            lead_in_directory(os.path.dirname(os.path.dirname(pbs_script)) or os.curdir)
            job_dir = os.path.dirname(os.path.abspath(pbs_script))
        pbs_content = render_pbs(opt, stack="grt" if "coupling" in args else "gr", handoff=lines,
                                 job_dir=job_dir)
        pbs_scripts.append(pbs_script)
        with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)

    # Submit-all script: each job is submitted from its own directory and waits on the previous one.
    submit_all_jobs_script = f"{options['simulation']['job_name']}_pbs.sh"
    if not pbs_scripts:
        raise ValueError(
            "No PBS scripts were generated (no .xml decks); cannot build the submit-all "
            "script. Check that segments/decks were produced.")
    jg = JobGen(CONFIG_DIR, machine=options["simulation"]["hpc_system"],
                node_type=options["simulation"].get("node_type"))
    ids = [os.path.basename(p)[:-len(".pbs")] for p in pbs_scripts]
    entries = submit_entries(pbs_scripts)
    lead = [e for e, i in zip(entries, ids) if is_lead_in(i)]
    rest = [e for e, i in zip(entries, ids) if not is_lead_in(i)]
    chains = [("GR lead-in", "job_id", lead, [])] if lead else []
    if rest:
        chains.append(("segments" if lead or len(rest) > 1 else None, "job_id", rest,
                       ["$job_id"] if lead else []))
    with open(submit_all_jobs_script, "w", encoding="utf-8") as f:
        f.write(jg.submit_chain_script(
            chains, header=f"Submit the {options['simulation']['job_name']} GR run. Run from anywhere.",
            jobids=f"{options['simulation']['job_name']}.jobids"))
    os.chmod(submit_all_jobs_script, 0o755)

    warmup_pbs_scripts = []
    spinup_pbs_scripts = []
    if "coupling" in args:
        # by name, not position: the SPINUP segment exists only with doSpin=T
        spinup_pbs_scripts = [p for p in pbs_scripts
                              if os.path.basename(p).endswith("-SPINUP.pbs")]
        warmup_pbs_scripts = [p for p in pbs_scripts
                              if "-WARMUP-" in os.path.basename(p)]
    return pbs_scripts, submit_all_jobs_script,spinup_pbs_scripts, warmup_pbs_scripts


def makeitso(args: dict = None):
    """Configure and write a MAGE run; the entry point for the command line and engage.

    Returns (options, spinup_pbs_scripts, warmup_pbs_scripts), (options, volt_files, []) for GRX,
    or None for --rederive list.
    """
    local_args = copy.deepcopy(args_default)
    if args is not None:
        local_args.update(args)
    args = local_args
    # programmatic callers bypass argparse's choices
    args["mode"] = normalize_mode(args.get("mode"))
    rederive = parse_rederive(args.get("rederive"))
    if rederive == ("list",):
        print(rederive_groups_text())
        return None
    require_env("KAIJUHOME")

    if args["debug"]:
        print(f"args = {args}")
    clobber = args["clobber"]
    debug = args["debug"]
    options_path = args["options_path"]
    if rederive and not options_path:
        raise ValueError("--rederive needs -o <options JSON>")
    verbose = args["verbose"]

    option_descriptions = load_option_descriptions()
    if debug:
        print(f"option_descriptions = {option_descriptions}")

    option_descriptions = update_option_descriptions(option_descriptions, args)
    if debug:
        print(f"option_descriptions = {option_descriptions}")

    if options_path:
        if verbose:
            print(f"Reading run options from {options_path}.")
        options = read_options_json(options_path)
        problem = section_problem(options)
        if problem:
            raise ValueError(f"{options_path}: {problem}")
        problem = old_format_problem(options)
        if problem:
            raise ValueError(f"{options_path}: {problem}")
        resolve_json_run_directory(options, options_path)
        rederive_pop(options, rederive)
        options = prepare_replay(options, option_descriptions, args)
    elif args.get("tui"):
        _tf = _import_own_tui_form()
        options = None
        try:
            options = _tf.run_tui_wizard(args)
        except _tf._FormCancelled:
            print(f"{YELLOW}makeitso wizard cancelled (Ctrl-C); makeitso wrote nothing.{RESET}")
            raise SystemExit(130) from None
        except _tf._FormUnavailable as e:
            print(f"Wizard unavailable ({e}); using the linear prompts.")
            if debug:
                import traceback
                traceback.print_exc()
        if options is None:
            options = prompt_user_for_run_options(option_descriptions, args)
        else:
            options = _tf.apply_makeitso_accept_actions(options, args)
    else:
        options = prompt_user_for_run_options(option_descriptions, args)
    if debug:
        print(f"options = {options}")

    # All checks run before anything is written.
    validate_options(options, option_descriptions, args)
    check_helper_layout(options)

    check_node_type(options)
    check_pe_layout(options)
    settle_account(options, "the run options")
    apply_bcwind_tsfile(options)
    check_deck_strings(options)
    apply_option_rules(options, args)

    # GRX: preprocess and write the GRX .xml decks only; engage writes the coupled job script.
    if args.get("grx"):
        grx = args["grx"]
        derive_select(options)
        options["pbs"].setdefault(
            "select2", str(1 + int(options["pbs"].get("num_helpers", 0) or 0)))
        # validates the cadence; gen_grx_voltron plans again
        grx_segment_plan(options, grx["stephours"], grx["delaycouplinghrs"], grx["t0"], grx["t1"])
        warn_dtcouple(options)
        if args.get("before_write"):
            args["before_write"](options)
        check_bcwind_t0(options, args, grx["t0"])
        # absolute: we chdir into it, so a relative path would double-nest later writes
        run_dir = prepare_run_directory(grx["run_dir"])
        options["pbs"]["run_directory"] = grx["run_dir"]
        options["simulation"]["start_date"] = grx["t0"]
        options["simulation"]["stop_date"] = grx["t1"]
        pin_bcwind_file(options, run_dir)
        os.chdir(run_dir)
        check_preprocessing_inputs(options, args)
        stage_bcwind_file(options, args.get("clobber", False))
        if args.get("save_parameters", True):
            dump_run_options(PARAMETERS_JSON, options)
        run_preprocessing_steps(options, args)
        lead_in_input_links(options, args.get("clobber", False))
        volt_files = gen_grx_voltron(options, grx["stephours"], grx["delaycouplinghrs"],
                                     grx["t0"], grx["t1"], run_dir, args)
        return options, volt_files, []

    # cda2wind accepts only the canonical date spelling
    sim = options.get("simulation", {})
    for k in ("start_date", "stop_date"):
        if sim.get(k) is not None:
            sim[k] = parse_run_datetime(sim[k], f"simulation.{k}").strftime(RUN_DATETIME_FORMAT)

    if "coupling" in args and yn(options["simulation"]["use_segments"]) != "Y":
        raise ValueError(
            "Coupled (engage GRT) runs require simulation.use_segments='Y' — the GR "
            "spin-up/warm-up phases are only generated for segmented runs. Set "
            "use_segments=Y (with a segment_duration).")

    # validates and settles num_segments / dtRes before the JSON is written
    decks = plan_segments(options, args)
    warn_dtcouple(options)
    derive_select(options)
    check_queue_limits(options)
    check_standalone_walltime(options, args)

    if args.get("before_write"):
        args["before_write"](options)
    check_bcwind_t0(options, args)
    # every write below is spelled from inside the run directory
    run_dir = prepare_run_directory(options["pbs"]["run_directory"])
    err = restart_files_error(options, run_dir, clobber)
    if err:
        raise FileExistsError(err)
    pin_bcwind_file(options, run_dir)
    launch_dir = os.getcwd()
    os.chdir(run_dir)

    check_preprocessing_inputs(options, args)
    check_voltron_exe(options, args)

    save = args.get("save_parameters", True)
    path = PARAMETERS_JSON
    if save and os.path.exists(path) and not clobber:
        raise FileExistsError(f"Options file {path} exists!")
    stage_bcwind_file(options, clobber)
    if save:
        dump_run_options(path, options)

    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(options, args)
    if any(is_lead_in(d["segment_id"]) for d in decks):
        lead_in_input_links(options, clobber)

    if verbose:
        print("Creating .ini file(s) for run.")
    ini_files = create_ini_files(options, args)
    if debug:
        print(f"ini_files = {ini_files}")

    if verbose:
        print("Converting .ini file(s) to .xml file(s).")
    xml_files = convert_ini_to_xml(ini_files, args)
    if debug:
        print(f"xml_files = {xml_files}")

    if verbose:
        print("Creating PBS job script(s) for run.")
    pbs_scripts, all_jobs_script, spinup_pbs_scripts, warmup_pbs_scripts = create_pbs_scripts(
        xml_files, options, args)
    if verbose:
        print(f"The PBS job scripts {pbs_scripts} are ready.")
    print(f"The PBS scripts {pbs_scripts} have been created, each with a "
          "corresponding XML file. To submit the jobs with the proper "
          "dependency (to ensure each segment runs in order), please run the "
          f"script {all_jobs_script} like this:\n"
          f"bash {all_jobs_script}")
    if save:
        saved = saved_json_spelling(run_dir, launch_dir)
        print(f"The options are saved in {saved}; replay them with makeitso.py -o {saved}.")

    # year-boundary splits may have added coupled segments
    if "coupling" in args:
        options["pbs"]["num_segments"] = str(len(pbs_scripts)-len(spinup_pbs_scripts)-len(warmup_pbs_scripts))
    return options, spinup_pbs_scripts, warmup_pbs_scripts


def main():
    """Command-line entry point."""
    parser = create_command_line_parser()

    args = parser.parse_args()
    args.tui = use_tui(args.no_tui)
    if args.debug:
        print(f"args = {args}")

    # a refusal prints one ERROR line (exit 1); --debug keeps the traceback
    try:
        makeitso(vars(args))
    except (ValueError, FileExistsError) as e:
        if args.debug:
            raise
        text = re.sub(r"\x1b\[[0-9;]*m", "", str(e))
        print(f"{RED}ERROR: {text}{RESET}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
