"""Arrow-navigable wizard for engage's own options (simulation, pbs, coupling).

Returns the same options dict as engage's linear prompts (`--no-tui`).
"""
import copy
import datetime
import importlib.util
import os
import sys

LEVELS = {"BASIC": 0, "INTERMEDIATE": 1, "EXPERT": 2}

# makeitso's tui_form engine, loaded by file path because tiegcmrun ships a tui_form.py too.
_ENGINE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "makeitso",
                            "tui_form.py")
_ENGINE_NAME = "engage_tui_engine"


def _engine():
    mod = sys.modules.get(_ENGINE_NAME)
    if mod is None:
        spec = importlib.util.spec_from_file_location(_ENGINE_NAME, os.path.normpath(_ENGINE_FILE))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[_ENGINE_NAME] = mod
        spec.loader.exec_module(mod)
    return mod


_FormUnavailable = _engine()._FormUnavailable
_FormCancelled = _engine()._FormCancelled


def _arg(args, key, default=None):
    if isinstance(args, dict):
        return args.get(key, default)
    return getattr(args, key, default)


def _is_set(x):
    s = str(x).strip().lower()
    return s not in ("", "none", "null", "[null]", "[none]", "[]")


def _try_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


_LIST_FIELDS = ("local_modules",)


def _parse_list_field(name, text):
    """An edited list field's text -> list, with the linear prompts' parser."""
    import makeitso
    return makeitso.as_list(text, tokens=name in makeitso.TOKEN_LIST_FIELDS)


class Field:
    __slots__ = ("name", "path", "level", "default", "valids", "prompt", "section",
                 "description", "warning")

    def __init__(self, name, path, level, default, valids, prompt, section,
                 description=None, warning=None):
        self.name, self.path, self.level = name, path, level
        self.default, self.valids, self.prompt, self.section = default, valids, prompt, section
        self.description, self.warning = description, warning


class FormState:
    def __init__(self, fields, derive, valid, mode="BASIC", context=None, cond=None, warn=None):
        self.fields = fields
        self.by_path = {f.path: f for f in fields}
        self.derive = derive
        self.valid = valid
        self.mode = mode
        self.cond = cond
        self.warn = warn or {}
        self.context = context or {}
        self.values = {f.path: ("" if f.default is None else f.default) for f in fields}
        self.overridden = set()
        self.derive_errors = {}
        # path -> fn(values) -> suggestions that left/right step through
        self.choices = {}
        self.rederive()

    def _all(self):
        d = dict(self.values)
        d.update(self.context)
        return d

    def visible(self):
        allv = self._all()
        out = []
        for f in self.fields:
            # not asked, but still emitted by options_dict
            if f.prompt is None:
                continue
            if LEVELS.get(f.level, 99) > LEVELS.get(self.mode, 99):
                continue
            if self.cond is not None and not self.cond(f, allv):
                continue
            out.append(f)
        return out

    def set(self, path, value):
        self.values[path] = self._normalize(path, value)
        self.overridden.add(path)
        self.rederive()

    def _normalize(self, path, value):
        """Strip whitespace from every field but the list field, as the linear prompts do."""
        f = self.by_path.get(path)
        if isinstance(value, str) and (f is None or f.name not in _LIST_FIELDS):
            return value.strip()
        return value

    def reset(self, path):
        """Ctrl-R: back to the derived or static default."""
        self.overridden.discard(path)
        f = self.by_path.get(path)
        if f is not None:
            self.values[path] = _static_default(f)
        self.rederive()

    def rederive(self):
        # two passes: a derive may read a field listed after it
        for _ in range(2):
            self.derive_errors = {}
            for f in self.fields:
                if f.path in self.overridden:
                    continue
                spec = self.derive.get(f.path)
                if not spec:
                    self.values[f.path] = _static_default(f)
                    continue
                fn, _deps = spec
                try:
                    self.values[f.path] = fn(self._all())
                except Exception as e:
                    self.derive_errors[f.path] = str(e)

    def errors(self):
        out = {}
        for f in self.visible():
            vfn = self.valid.get(f.path)
            if vfn:
                e = vfn(self._all(), self.values[f.path])
                if e:
                    out[f.path] = e
                    continue
            if f.valids:                               # static valids membership (case-sensitive)
                v = str(self.values[f.path]).strip()
                if v != "":
                    vv = [str(x).strip() for x in f.valids]
                    ok = v in vv
                    if not ok:                         # numeric tolerance ("2.0" matches 2)
                        fv = _try_float(v)
                        ok = fv is not None and fv in [g for g in (_try_float(x) for x in vv)
                                                       if g is not None]
                    if not ok:
                        out[f.path] = f"must be one of {f.valids}"
        return out

    def choices_for(self, path):
        fn = self.choices.get(path)
        if not fn:
            return []
        try:
            return [str(x) for x in (fn(self._all()) or [])]
        except Exception:
            return []

    def warnings_for(self, path):
        fn = self.warn.get(path)
        if not fn:
            return None
        try:
            return fn(self._all()) or None
        except Exception:
            return None

    def options_dict(self):
        """The nested options dict; values stay raw strings (T/F/Y/N are not coerced)."""
        out = {}
        for f in self.fields:
            cur = out
            parts = f.path.split(".")
            for p in parts[:-1]:
                cur = cur.setdefault(p, {})
            val = self.values[f.path]
            if isinstance(val, str) and f.name in _LIST_FIELDS:
                # a blank box on a null-default field is the untouched default (None)
                if val.strip() == "" and f.default is None:
                    val = None
                else:
                    val = _parse_list_field(f.name, val)
            elif isinstance(val, str):
                s = val.strip()
                if s == "":
                    val = None if f.default is None else val
                elif s.lower() in ("[null]", "[none]"):
                    val = [None]
            cur[parts[-1]] = val
        return out


def _static_default(f):
    return "" if f.default is None else f.default


class MakeitsoFormState(FormState):
    """FormState whose PBS machine block (pbs.*) is rebuilt when hpc_system changes."""
    def __init__(self, option_descriptions, derive, valid, mode, context, cond, warn):
        self._od = option_descriptions
        self._hpc = context.get("__hpc0")
        self._grid = context.get("__grid0")
        self._grx = bool(context.get("__grx"))
        self._machine_error = None
        fields = load_makeitso_fields(option_descriptions, self._hpc, self._grid, self._grx)
        super().__init__(fields, derive, valid, mode, context, cond, warn)

    def set(self, path, value):
        value = self._normalize(path, value)
        self.values[path] = value
        self.overridden.add(path)
        if path == "simulation.hpc_system" and value != self._hpc:
            self._rebuild_machine(value)
        self.rederive()

    def reset(self, path):
        f = self.by_path.get(path)
        if path == "simulation.hpc_system" and f is not None:
            self.overridden.discard(path)
            default = _static_default(f)
            self.values[path] = default
            if default != self._hpc:
                self._rebuild_machine(default)
            self.rederive()
            return
        super().reset(path)

    def errors(self):
        out = super().errors()
        if self._machine_error and "simulation.hpc_system" not in out:
            out["simulation.hpc_system"] = self._machine_error
        return out

    def _rebuild_machine(self, new_hpc):
        self._hpc = new_hpc
        old_vals, old_over = dict(self.values), set(self.overridden)
        grid = self.values.get("simulation.gamera_grid_type", self._grid)
        try:
            new_fields = load_makeitso_fields(self._od, new_hpc, grid, self._grx, strict=True)
        except Exception as e:
            self._machine_error = f"no PBS settings for hpc_system {new_hpc!r}: {e}"
            return
        self._machine_error = None
        self.fields = new_fields
        self.by_path = {f.path: f for f in new_fields}
        nv = {}
        for f in new_fields:
            if f.path in old_over:
                nv[f.path] = old_vals.get(f.path)
            else:
                nv[f.path] = "" if f.default is None else f.default
        self.values = nv
        self.overridden = {p for p in old_over if p in self.by_path}


def _iso(s):
    import makeitso
    return makeitso.parse_run_datetime(s)


def _field_from(meta, name, path, section):
    return Field(name, path, meta.get("LEVEL", "EXPERT"), meta.get("default"),
                 meta.get("valids"), meta.get("prompt", name), section,
                 meta.get("description"), meta.get("warning"))


# The simulation fields engage's linear prompts collect; the rest belong to makeitso.
_ENGAGE_SIM_FIELDS = ("job_name", "start_date", "stop_date", "use_segments",
                      "segment_duration", "gamera_grid_type", "hpc_system", "node_type")
# Asked after the PBS block: their check and default depend on the walltime.
_SEGMENT_FIELDS = ("use_segments", "segment_duration")


def load_makeitso_fields(option_descriptions, hpc_system, gamera_grid_type, grx=False, strict=False):
    """The ordered Field list for the chosen machine and grid."""
    import makeitso
    fields = []

    def is_field(meta):
        return isinstance(meta, dict) and "LEVEL" in meta and ("prompt" in meta or "default" in meta)

    for name, meta in option_descriptions.get("simulation", {}).items():
        if name == "node_type":
            # only on a machine with several node types
            try:
                meta = makeitso.node_type_description(option_descriptions, hpc_system)
            except Exception:
                meta = None
            if meta is None:
                continue
        if name in _ENGAGE_SIM_FIELDS and name not in _SEGMENT_FIELDS and is_field(meta):
            fields.append(_field_from(meta, name, f"simulation.{name}", "simulation"))
    for name, meta in option_descriptions.get("pbs", {}).get("_common", {}).items():
        if is_field(meta):
            fields.append(_field_from(meta, name, f"pbs.{name}", "pbs"))

    try:
        mblock = makeitso.synthesize_pbs_descriptions(hpc_system, option_descriptions)
    except Exception:
        if strict:
            raise
        mblock = option_descriptions.get("pbs", {}).get(hpc_system, {})
    # engage's coupled-run local_modules field replaces makeitso's where both exist
    try:
        import grx_options
        _here = os.path.dirname(os.path.abspath(__file__))
        _eng_pbs = grx_options.load_engage_descriptions(
            os.path.join(_here, "option_engage_descriptions.json")).get("pbs", {}).get(hpc_system, {})
    except Exception:
        _eng_pbs = {}
    for name, meta in (mblock or {}).items():
        if not is_field(meta):
            continue
        if name == "local_modules" and name in _eng_pbs:
            meta = _eng_pbs[name]
        d = dict(meta)
        # grid-keyed default, e.g. num_helpers {D:..,Q:..}
        if isinstance(d.get("default"), dict):
            d["default"] = d["default"].get(gamera_grid_type, next(iter(d["default"].values())))
        if grx and name == "walltime":
            import grx_options
            d["default"] = grx_options.grx_defaults(hpc_system)["pbs.walltime"]
        fields.append(_field_from(d, name, f"pbs.{name}", "pbs"))

    for name in _SEGMENT_FIELDS:
        meta = option_descriptions.get("simulation", {}).get(name)
        if is_field(meta):
            fields.append(_field_from(meta, name, f"simulation.{name}", "simulation"))

    # component blocks (gamera, voltron, ...) are left to the makeitso wizard
    for name, meta in option_descriptions.get("coupling", {}).items():
        if is_field(meta):
            fields.append(_field_from(meta, name, f"coupling.{name}", "coupling"))
    return fields


_ENGAGE_CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")


def _account_required(hpc):
    """engage.account_required without importing engage (this module runs inside its __main__)."""
    import makeitso
    return makeitso.account_required(hpc, _ENGAGE_CONFIG)


def _walltime_check(component):
    """Walltime validator; an unsegmented GRT run's whole window must fit too."""
    def check(v, val):
        import grx_options
        seg = None
        if component == "grt" and str(v.get("simulation.use_segments", "Y")).upper() != "Y":
            try:
                seg = (_iso(v["simulation.stop_date"]) - _iso(v["simulation.start_date"])).total_seconds()
            except (KeyError, TypeError, ValueError):
                seg = None
        return grx_options.walltime_error(seg, val, str(v.get("simulation.hpc_system", "")),
                                          v.get("pbs.queue"), component,
                                          v.get("simulation.gamera_grid_type"))
    return check


def _segment_choices(v, component):
    import grx_options
    return grx_options.segment_choices(
        component, str(v.get("simulation.hpc_system", "")), v.get("pbs.queue"), v.get("pbs.walltime"),
        v.get("simulation.gamera_grid_type"), v.get("simulation.start_date"),
        v.get("simulation.stop_date"), v.get("coupling.gr_warm_up_time") if component == "grx" else None)


def _segment_derive(component, fallback):
    """segment_duration default: the largest valid length in seconds, else `fallback(v)`."""
    def fn(v):
        import grx_options
        try:
            valid, info, _ = _segment_choices(v, component)
        except Exception:
            valid, info = [], {}
        return grx_options.suggested_seconds(valid, info) if valid else fallback(v)
    return fn


_SEGMENT_DEPS = ["simulation.start_date", "simulation.stop_date", "simulation.gamera_grid_type",
                 "simulation.hpc_system", "pbs.queue", "pbs.walltime"]


def _account_check(v, val):
    hpc = str(v.get("simulation.hpc_system", ""))
    if not _account_required(hpc):
        return None
    import makeitso
    return makeitso.account_prompt_check(hpc)(val)


_HELPERS = {"D": "0", "Q": "2", "O": "4", "H": "8"}


def _derive_registry():
    import makeitso

    def _grid(v):
        return str(v.get("simulation.gamera_grid_type", "Q"))

    def _hpc(v):
        return str(v.get("simulation.hpc_system", "aitken"))

    def _eng(v, key):
        cfg = v.get("__engage_cfg")
        return (cfg.get(key) if isinstance(cfg, dict) else None)

    def _threads(v, ranks_path):
        import makeitso
        nt = str(v.get("simulation.node_type") or "").strip() or None
        t = makeitso.thread_default(_hpc(v), nt, v.get(ranks_path))
        if t is None:
            raise ValueError(f"{ranks_path} is not a positive integer")
        return t

    def _sim_seconds(v):
        return (_iso(v["simulation.stop_date"]) - _iso(v["simulation.start_date"])).total_seconds()

    def _num_segments(v):
        if str(v.get("simulation.use_segments", "N")).upper() != "Y":
            return "1"
        import math
        return str(math.ceil(_sim_seconds(v) / float(v["simulation.segment_duration"])))

    return {
        # fallback: the whole run, as a float string (e.g. 7200.0)
        "simulation.segment_duration": (_segment_derive("grt", lambda v: str(_sim_seconds(v))),
                                        list(_SEGMENT_DEPS)),
        "pbs.account_name": (lambda v: (_eng(v, "account_name")
                                        or (makeitso.account_prompt_default(_hpc(v), _ENGAGE_CONFIG)
                                            if _account_required(_hpc(v)) else None)),
                             ["simulation.hpc_system"]),
        "pbs.kaiju_install_directory": (lambda v: _eng(v, "kaiju_install_directory")
                                        or os.environ.get("KAIJUHOME", ""), []),
        "pbs.kaiju_build_directory": (lambda v: _eng(v, "kaiju_build_directory")
                                      or makeitso.kaiju_build_default("grx" if v.get("__grx") else "grt"), []),
        "coupling.conda_env": (lambda v: os.environ.get("CONDA_DEFAULT_ENV", ""), []),
        "pbs.num_segments": (_num_segments, ["simulation.use_segments",
                             "simulation.segment_duration", "simulation.start_date", "simulation.stop_date"]),
        "pbs.num_helpers": (lambda v: _HELPERS.get(_grid(v), "2"), ["simulation.gamera_grid_type"]),
        # cores per node // ranks per node
        "pbs.ompthreads": (lambda v: _threads(v, "pbs.mpiprocs"),
                           ["simulation.hpc_system", "simulation.node_type", "pbs.mpiprocs"]),
        "pbs.helper_ompthreads": (lambda v: _threads(v, "pbs.helper_mpiprocs"),
                                  ["simulation.hpc_system", "simulation.node_type",
                                   "pbs.helper_mpiprocs"]),
    }


def _validate_registry():
    def iso(v, val):
        try:
            _iso(val)
            return None
        except (ValueError, TypeError):
            return "Invalid date (expected yyyy-mm-ddThh:mm:ss)"

    def stop_after(v, val):
        err = iso(v, val)
        if err:
            return err
        try:
            if _iso(val) <= _iso(v["simulation.start_date"]):
                return "stop_date must be after start_date"
        except (ValueError, TypeError, KeyError):
            return None
        if v.get("__grx"):
            return None
        import grx_options
        return grx_options.grt_window_error(v.get("simulation.start_date"), val)

    def seg_dur(v, val):
        if str(v.get("simulation.use_segments", "N")).upper() != "Y" and not v.get("__grx"):
            return None
        f = _try_float(val)
        if f is None or f <= 0:
            return "segment_duration must be a positive number"
        if v.get("__grx"):
            return None                               # the GRX registry's own check follows
        import grx_options
        return grx_options.segment_hours_error(
            "grt", str(v.get("simulation.hpc_system", "")), v.get("pbs.queue"), v.get("pbs.walltime"),
            v.get("simulation.gamera_grid_type"), v.get("simulation.start_date"),
            v.get("simulation.stop_date"), val)

    return {
        "simulation.start_date": iso,
        "simulation.stop_date": stop_after,
        "simulation.segment_duration": seg_dur,
        "pbs.account_name": _account_check,
        "pbs.walltime": _walltime_check("grt"),
    }


def _grx_validate_registry():
    """GRX validators. BASIC hides gr_warm_up_time, so there the segment check also covers it."""
    import grx_options
    base = _validate_registry()

    def machine(v, val):
        return grx_options.machine_error(str(val).strip())

    def segment(v, val):
        err = base["simulation.segment_duration"](v, val)
        if err:
            return err
        hpc = str(v.get("simulation.hpc_system", ""))
        err = (grx_options.segment_hours_error(
                   "grx", hpc, v.get("pbs.queue"), v.get("pbs.walltime"),
                   v.get("simulation.gamera_grid_type"), v.get("simulation.start_date"),
                   v.get("simulation.stop_date"), val, v.get("coupling.gr_warm_up_time"))
               or grx_options.segment_error(v.get("simulation.start_date"), v.get("simulation.stop_date"),
                                            val, v.get("coupling.gr_warm_up_time"), v.get("pbs.walltime"),
                                            hpc, v.get("pbs.queue")))
        if err or LEVELS.get("INTERMEDIATE") <= LEVELS.get(v.get("__mode"), 0):
            return err
        return grx_options.warmup_error(v.get("coupling.gr_warm_up_time"), val)

    def warmup(v, val):
        return grx_options.warmup_error(val, v.get("simulation.segment_duration"),
                                        v.get("simulation.start_date"), v.get("simulation.stop_date"),
                                        v.get("pbs.walltime"), str(v.get("simulation.hpc_system", "")),
                                        v.get("pbs.queue"))

    def start(v, val):
        return base["simulation.start_date"](v, val) or grx_options.coupling_start_error(val)

    return dict(base, **{"simulation.hpc_system": machine, "simulation.segment_duration": segment,
                         "coupling.gr_warm_up_time": warmup, "simulation.start_date": start,
                         "pbs.walltime": _walltime_check("grx")})


def _cond_registry(grx=False):
    def cond(f, v):
        p = f.path
        if p == "simulation.use_segments" and grx:
            return False                               # a GRX run is always segmented
        if p == "simulation.segment_duration":
            return grx or str(v.get("simulation.use_segments", "N")).upper() == "Y"
        if p == "coupling.conda_env" and grx:
            return False                               # GRX: CONDA_DEFAULT_ENV, never asked
        if p == "pbs.account_name":
            return _account_required(str(v.get("simulation.hpc_system", "")))
        return True
    return cond


def _warn_registry():
    return {}


def _grx_warn_registry():
    return dict(_warn_registry())


def build_state(args, mode="BASIC"):
    """The engage FormState (simulation, pbs and coupling fields), with GRX defaults for --component grx."""
    import json
    import makeitso
    here = os.path.dirname(os.path.abspath(__file__))
    option_descriptions = makeitso.load_option_descriptions()
    import grx_options
    engage_od = grx_options.load_engage_descriptions(os.path.join(here, "option_engage_descriptions.json"))
    option_descriptions["coupling"] = engage_od.get("coupling", {})   # engage-only section
    if "use_segments" in option_descriptions.get("simulation", {}):
        option_descriptions["simulation"]["use_segments"]["default"] = "Y"
    grx = _arg(args, "component") == "grx"
    # only machines with an engage PBS block
    hod = option_descriptions["simulation"]["hpc_system"]
    hod["valids"] = [m for m in makeitso.machine_names() if m in engage_od.get("pbs", {})]
    if hod.get("default") not in hod["valids"]:
        hod["default"] = hod["valids"][0] if hod["valids"] else None
    god = option_descriptions["simulation"]["gamera_grid_type"]
    god["valids"] = list(makeitso.RUN_GRID_TYPES["grx" if grx else "grt"])
    if god.get("default") not in god["valids"]:
        god["default"] = god["valids"][0]
    derive, valid, cond, warn = (_derive_registry(), _validate_registry(), _cond_registry(),
                                 _warn_registry())
    if grx:
        import grx_options
        gd = grx_options.grx_defaults()
        sod = option_descriptions["simulation"]
        sod["start_date"]["default"] = gd["simulation.start_date"]
        sod["stop_date"]["default"] = gd["simulation.stop_date"]
        sod["hpc_system"]["default"] = gd["simulation.hpc_system"]
        sod["hpc_system"]["valids"] = grx_options.grx_machines()
        option_descriptions["coupling"]["gcm_spin_up_time"]["default"] = gd["coupling.gcm_spin_up_time"]
        derive["simulation.segment_duration"] = (
            _segment_derive("grx", lambda v: gd["simulation.segment_duration"]),
            list(_SEGMENT_DEPS) + ["coupling.gr_warm_up_time"])
        valid, cond, warn = _grx_validate_registry(), _cond_registry(grx=True), _grx_warn_registry()
    hpc0 = option_descriptions["simulation"]["hpc_system"].get("default", "aitken")
    grid0 = option_descriptions["simulation"]["gamera_grid_type"].get("default", "Q")
    ctx = {"__engage_cfg": None, "__engage": False, "__coupling": True,
           "__hpc0": hpc0, "__grid0": grid0, "__grx": grx, "__mode": mode}
    state = MakeitsoFormState(option_descriptions, derive, valid, mode, ctx, cond, warn)
    state.choices = {"pbs.account_name": lambda v: makeitso.account_choices(
        str(v.get("simulation.hpc_system", "")), _ENGAGE_CONFIG)}
    return state


def _seed_engage(state, eng):
    """Pin the values a replay supplied so derives do not overwrite them."""
    if not isinstance(eng, dict):
        return
    for sec in ("simulation", "pbs"):
        block = eng.get(sec)
        if not isinstance(block, dict):
            continue
        for k, val in block.items():
            if val is None:
                continue
            path = f"{sec}.{k}"
            if path in state.by_path:
                state.set(path, val)


def apply_makeitso_accept_actions(options, args):
    """Drop an unused bcwind_file and add pbs.select2 (1 + num_helpers), as the linear prompts do."""
    sim = options.get("simulation", {})
    pbs = options.get("pbs", {})
    if str(sim.get("bcwind_available", "N")).upper() != "Y":
        sim.pop("bcwind_file", None)
    try:
        pbs["select2"] = str(1 + int(pbs.get("num_helpers", 0)))
    except (TypeError, ValueError):
        pass
    return options


def run_tui_wizard(args):
    """Run the wizard and return the options dict.

    Raises _FormUnavailable without prompt_toolkit or a TTY, _FormCancelled (a subclass) on Ctrl-C.
    """
    eng = _engine()
    eng.require_tui()
    mode = _arg(args, "mode") or "BASIC"
    state = build_state(args, mode=mode)
    stack = "GRX" if state.context.get("__grx") else "GRT"
    return eng.run_form(state, f"engage ({stack})")


def apply_engage_accept_actions(options, args):
    """Drop a stray bcwind_file; no pbs.select2 here (the coupled PBS reads makeitso's)."""
    options.get("simulation", {}).pop("bcwind_file", None)
    return options


def run_engage_wizard(args):
    """Run the engage wizard and return the engage_options dict."""
    return apply_engage_accept_actions(run_tui_wizard(args), args)
