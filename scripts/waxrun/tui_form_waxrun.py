"""Arrow-navigable wizard for waxrun's options; asks what waxrun.prompt_options would ask.

Uses makeitso's tui_form engine, loaded by file path because tiegcmrun also ships a tui_form.py.
"""
import contextlib
import copy
import importlib.util
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ENGINE_FILE = os.path.join(_HERE, os.pardir, "makeitso", "tui_form.py")
_ENGINE_NAME = "waxrun_tui_engine"


def _engine():
    """makeitso's tui_form module, loaded once."""
    mod = sys.modules.get(_ENGINE_NAME)
    if mod is None:
        spec = importlib.util.spec_from_file_location(_ENGINE_NAME, os.path.normpath(_ENGINE_FILE))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[_ENGINE_NAME] = mod
        spec.loader.exec_module(mod)
    return mod


def _waxrun():
    import waxrun
    return waxrun


_FormUnavailable = _engine()._FormUnavailable
_FormCancelled = _engine()._FormCancelled

# asked last, defaults derived from the automatic restart pick
START_FILE = "ic.start_file_date"
SPINUP_START = "ic.spinup_start"
IC_DATES = (START_FILE, SPINUP_START)
BOOL_VALIDS = ["Y", "N"]


def _shown(value, spec):
    if isinstance(spec.get("default"), bool):
        return "" if value is None else ("Y" if value else "N")
    return "" if value is None else str(value)


def _prompted_leaves(desc, opt, mode, hidden):
    """(path, key, spec, section) of every leaf waxrun.prompt_options would ask, in its order."""
    wx = _waxrun()
    out, last = [], []

    def walk(dsec, prefix, section):
        for key, spec in dsec.items():
            if key.startswith("_") or not isinstance(spec, dict):
                continue
            path = f"{prefix}.{key}"
            if "LEVEL" in spec:
                if not spec.get("prompt") or path in hidden or not wx._at_or_below(spec["LEVEL"], mode):
                    continue
                if path in IC_DATES:
                    last.append((path, key, spec, section))
                else:
                    out.append((path, key, spec, section))
            else:
                walk(spec, path, section)

    for section, dsec in desc.items():
        if section.startswith("_") or section == "version" or not isinstance(dsec, dict):
            continue
        walk(dsec, section, section)
    return out + sorted(last, key=lambda leaf: IC_DATES.index(leaf[0]))


def _get(opt, path):
    cur = opt
    for part in path.split("."):
        if not isinstance(cur, dict):
            return None
        cur = cur.get(part)
    return cur


def _put(opt, path, value):
    parts = path.split(".")
    cur = opt
    for part in parts[:-1]:
        cur = cur.setdefault(part, {})
    cur[parts[-1]] = value


def srcroot_error(value, run_dir=None):
    """None, or why value cannot be simulation.srcroot; a relative value is read against run_dir."""
    s = os.path.expanduser(str(value or "").strip())
    if s and run_dir and not os.path.isabs(s):
        s = os.path.normpath(os.path.join(run_dir, s))
    if not s:
        return ("simulation.srcroot is required: the CAM tree you cloned; both WACCM-X cases are "
                "created from it (export WAXHOME=<tree> to make it the default)")
    waxhome = _waxrun().tree_solver.waxhome()
    origin = " (the $WAXHOME default)" if os.path.normpath(os.path.abspath(s)) == waxhome else ""
    if not os.path.isdir(s):
        return f"{s}{origin} is not a directory"
    if not os.path.isfile(os.path.join(s, "cime", "scripts", "create_newcase")):
        return f"{s}{origin} is not a CAM clone (no cime/scripts/create_newcase)"
    return None


class WaxrunFormState(_engine().FormState):
    """FormState over waxrun leaves; yes/no leaves are shown as Y/N."""

    def __init__(self, opt, desc, mode, hidden=()):
        eng = _engine()
        self._opt, self._desc = opt, desc
        self._specs = {}
        leaves = _prompted_leaves(desc, opt, mode, set(hidden))
        names = {}                   # over all modes, so a name never depends on the mode
        for _, key, _, _ in _prompted_leaves(desc, opt, "EXPERT", set()):
            names[key] = names.get(key, 0) + 1
        fields = []
        self._pick_cache, self._choices_cache = {}, {}
        for path, key, spec, section in leaves:
            self._specs[path] = spec
            # the full path when two leaves share a name; the simulation leaf keeps the bare name
            name = key if names[key] == 1 or path == f"simulation.{key}" else path
            is_bool = isinstance(spec.get("default"), bool)
            shown = _shown(_get(opt, path), spec)
            if path == START_FILE and shown:
                shown = _waxrun().resolve_start_file_date(opt, shown) or shown
            if path == "simulation.project" and shown in (None, ""):
                shown = _waxrun().project_prompt_default(self._hpc())
            fields.append(eng.Field(name, path, spec["LEVEL"], shown,
                                    BOOL_VALIDS if is_bool else spec.get("valids"),
                                    spec["prompt"], section, spec.get("description"),
                                    spec.get("warning")))
        valid = {path: self._validator(path) for path in self._specs}
        fns = {START_FILE: self._start_file_default, SPINUP_START: self._spinup_start_default}
        derive = {p: (fn, ()) for p, fn in fns.items() if p in self._specs and _get(opt, p) is None}
        super().__init__(fields, derive, valid, mode, {}, None,
                         {START_FILE: self._start_file_note, SPINUP_START: self._spinup_start_note})
        # offered, not enforced
        self.choices = {"simulation.project": lambda v: _waxrun().project_choices(self._hpc())}
        for path, spec in self._specs.items():
            if spec.get("type") == "sourcemods":
                self.choices[path] = lambda v: _waxrun().sourcemods_choices()

    def _hpc(self):
        return (self._opt.get("simulation") or {}).get("hpc_system", "derecho")

    def visible(self):
        """The fields of the mode that waxrun.asked keeps given the answers so far."""
        wx = _waxrun()
        return [f for f in self.fields if wx._at_or_below(f.level, self.mode)
                and wx.asked(f.path, self._section_now(f.path))]

    def _section_now(self, path):
        prefix = path.rsplit(".", 1)[0]
        sec = dict(_get(self._opt, prefix) or {})
        for p, v in self.values.items():
            if p.rsplit(".", 1)[0] == prefix:
                sec[p.rsplit(".", 1)[1]] = v
        return sec

    def set(self, path, value):
        spec = self._specs.get(path, {})
        if isinstance(spec.get("default"), bool):
            b = _waxrun()._bool_of(value)
            if b is not None:
                value = "Y" if b else "N"
        if path == START_FILE and str(value).strip():
            value = _waxrun().resolve_start_file_date(self.current_options(), str(value).strip()) or value
        super().set(path, value)

    def rederive(self):
        f = self.by_path.get(START_FILE)
        if f is not None:
            f.valids = self._start_file_choices() or None
        super().rederive()

    def _start_file_choices(self):
        """The restart times of ic.start_file_dir as answered so far."""
        opt = self.current_options()
        key = repr(_get(opt, "ic.start_file_dir"))
        if key not in self._choices_cache:
            self._choices_cache[key] = _waxrun().start_file_choices(opt)
        return self._choices_cache[key]

    def _convert(self, path, text):
        """The shown text of one leaf -> its value; blank or none/null is None."""
        wx = _waxrun()
        spec = self._specs[path]
        s = str(text).strip()
        numeric = isinstance(spec.get("default"), (int, float)) and not isinstance(spec.get("default"), bool)
        if s != "" and spec.get("type") == "sourcemods":
            return wx._leaf_value(path, spec, s)
        if s == "" or (not numeric and s.lower() in ("none", "null")):
            if spec.get("default") is not None:
                raise ValueError(f"{path} needs a value (Ctrl-R restores the default)")
            return None
        return wx._leaf_value(path, spec, s, run_dir=self._run_dir())

    def _run_dir(self):
        """The run folder as answered so far."""
        text = self.values.get("simulation.run_directory")
        if text is None or not str(text).strip():
            return _waxrun()._run_dir_of(self._opt)
        return os.path.abspath(os.path.expanduser(str(text).strip()))

    def _validator(self, path):
        def check(_values, text):
            if path == "simulation.srcroot":
                return srcroot_error(text, self._run_dir())
            if path == "simulation.project" and str(text).strip() in ("", "None", "none", "null"):
                hpc = (self._opt.get("simulation") or {}).get("hpc_system", "derecho")
                err = _waxrun().project_error(None, hpc)
                if err:
                    return err
            try:
                value = self._convert(path, text)
            except ValueError as e:
                return str(e)
            if path in IC_DATES:
                return self._ic_date_error(path, value)
            return None
        return check

    def answers(self):
        """The nested answers of every shown leaf; an IC date left at its derived default is None."""
        out = {}
        for f in self.visible():
            if f.path in self.derive and f.path not in self.overridden:
                _put(out, f.path, None)
                continue
            _put(out, f.path, self._convert(f.path, self.values[f.path]))
        return out

    def current_options(self):
        """The options dict with the answers so far applied."""
        opt = copy.deepcopy({k: v for k, v in self._opt.items() if k != "_engage"})
        opt["_engage"] = self._opt.get("_engage")
        for f in self.visible():
            if f.path in IC_DATES:
                continue
            try:
                _put(opt, f.path, self._convert(f.path, self.values[f.path]))
            except ValueError:
                pass
        return opt

    def _start_file(self):
        """The start file date answered so far; None while it is the automatic pick."""
        if START_FILE not in self._specs or (START_FILE in self.derive and START_FILE not in self.overridden):
            return _get(self._opt, START_FILE)
        try:
            return self._convert(START_FILE, self.values[START_FILE])
        except ValueError:
            return None

    def _auto_pick(self, start_file_date=None):
        """(day, spinup_start, note) of waxrun.auto_pick with the answers so far, output silenced."""
        opt = self.current_options()
        key = repr([_get(opt, p) for p in ("ic.start_file_dir", "spinup.target_seconds", "spinup.min_hours",
                                            "spinup.max_hours", "spinup.dtime", "simulation.res",
                                            "simulation.coupling_start", "simulation.stop_time",
                                            "coupled.stephours", "mage.delaycouplinghrs")]
                   + [start_file_date])
        if key not in self._pick_cache:
            sink = io.StringIO()
            try:
                with contextlib.redirect_stderr(sink), contextlib.redirect_stdout(sink):
                    self._pick_cache[key] = _waxrun().auto_pick(opt, start_file_date)
            except (SystemExit, Exception) as e:
                self._pick_cache[key] = (None, None, _waxrun().auto_pick_failure(e))
        return self._pick_cache[key]

    def _start_file_default(self, _values):
        day = self._auto_pick()[0]
        if day is None:
            raise _engine()._NotApplicable()
        return day

    def _spinup_start_default(self, _values):
        start = self._auto_pick(self._start_file())[1]
        if start is None:
            raise _engine()._NotApplicable()
        return start

    def _start_file_note(self, _values):
        sfd = self._start_file()
        return self._auto_pick(sfd)[2] if sfd is not None else self._auto_pick()[2]

    def _spinup_start_note(self, _values):
        """The spinup length the answer gives, or why the default was moved."""
        if SPINUP_START in self.derive and SPINUP_START not in self.overridden:
            note = self._auto_pick(self._start_file())[2]
            if note and "default moved" in note:
                return note
        try:
            value = self._convert(SPINUP_START, self.values.get(SPINUP_START, ""))
            cs = _waxrun()._run_dt(self.current_options()["simulation"]["coupling_start"], "coupling_start")
        except (ValueError, SystemExit, KeyError):
            return None
        if value is None:
            value = self._auto_pick(self._start_file())[1]
        if value is None:
            return None
        hours = (cs - _waxrun()._run_dt(value, "ic.spinup_start")).total_seconds() / 3600
        return f"{hours:g} h of spinup"

    def _ic_date_error(self, path, value):
        opt = self.current_options()
        default = self._auto_pick(self._start_file())[1] if path == SPINUP_START else None
        sink = io.StringIO()
        try:
            with contextlib.redirect_stderr(sink), contextlib.redirect_stdout(sink):
                return _waxrun().ic_date_error(opt, path, value, default)
        except (SystemExit, Exception):
            return None


def build_state(opt, mode="BASIC", desc=None, hidden=()):
    """The waxrun FormState for opt; engage-filled and hidden leaves are not asked."""
    wx = _waxrun()
    desc = desc or wx.load_descriptions()
    hide = set(opt.get("_engage_filled", ())) | set(hidden)
    return WaxrunFormState(opt, desc, mode, hide)


def run_waxrun_wizard(opt, mode="BASIC", desc=None, hidden=()):
    """Run the waxrun wizard and return the nested answers.

    Raises _FormUnavailable without a TTY or prompt_toolkit, _FormCancelled on Ctrl-C."""
    eng = _engine()
    eng.require_tui()
    state = build_state(opt, mode, desc, hidden)
    if not state.visible():
        return {}
    return eng.run_form(state, "waxrun", collect=lambda st: st.answers())
