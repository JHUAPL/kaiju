"""Arrow-navigable option-entry wizard for makeitso, and the screen (run_form) engage and waxrun reuse.

Returns the same options dict as the linear prompts (`--no-tui`).
"""
import copy
import datetime
import os

LEVELS = {"BASIC": 0, "INTERMEDIATE": 1, "EXPERT": 2}


class _FormUnavailable(Exception):
    """Raised when the TUI can't run (no prompt_toolkit / no TTY)."""


class _FormCancelled(_FormUnavailable):
    """Raised when the user cancels the wizard (Ctrl-C)."""


class _NotApplicable(Exception):
    """Raised by a derive that does not apply to the current answers: the field keeps its default."""


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
                 "description", "warning", "meta")

    def __init__(self, name, path, level, default, valids, prompt, section,
                 description=None, warning=None, meta=None):
        self.name, self.path, self.level = name, path, level
        self.default, self.valids, self.prompt, self.section = default, valids, prompt, section
        self.description, self.warning = description, warning
        self.meta = meta


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
        # fn(field, value) -> None or a type error; checked after the field's own validator
        self.coerce_error = None
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
            if f.path in self.context.get("__engage_provided", ()):
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
        self.derive_errors = {}
        if self.cond is not None:
            # a derived field hidden by another answer drops its typed value
            allv = self._all()
            for p in list(self.overridden):
                f = self.by_path.get(p)
                if f is not None and p in self.derive and not self.cond(f, allv):
                    self.overridden.discard(p)
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
                except _NotApplicable:
                    self.values[f.path] = _static_default(f)
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
            if self.coerce_error is not None:
                e = self.coerce_error(f, self.values[f.path])
                if e:
                    out[f.path] = e
                    continue
            if f.valids:
                v = str(self.values[f.path]).strip()
                if v != "":
                    vv = [str(x).strip() for x in f.valids]
                    # exact spelling only: the value is emitted raw
                    ok = v in vv
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
        vis = set(id(f) for f in self.visible())
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
            elif id(f) in vis and val is not None and not isinstance(val, list):
                # the linear prompts str() an asked field; a hidden one keeps its raw type
                val = str(val)
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
        self._machine_error = None
        fields = load_makeitso_fields(option_descriptions, self._hpc, self._grid)
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
            new_fields = load_makeitso_fields(self._od, new_hpc, grid, strict=True)
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
                 meta.get("description"), meta.get("warning"), meta)


# Asked after the PBS block: their check and default depend on the walltime.
_SEGMENT_FIELDS = ("use_segments", "segment_duration")


def load_makeitso_fields(option_descriptions, hpc_system, gamera_grid_type, strict=False):
    """The ordered Field list for the chosen machine and grid; strict raises when no PBS block can be built."""
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
        if name not in _SEGMENT_FIELDS and is_field(meta):
            fields.append(_field_from(meta, name, f"simulation.{name}", "simulation"))
    for name, meta in option_descriptions.get("pbs", {}).get("_common", {}).items():
        if is_field(meta):
            fields.append(_field_from(meta, name, f"pbs.{name}", "pbs"))

    # get_pbs_descriptions returns the block engage may already have cached, else synthesizes it
    try:
        mblock = makeitso.get_pbs_descriptions(option_descriptions, hpc_system)
    except Exception:
        if strict:
            raise
        mblock = option_descriptions.get("pbs", {}).get(hpc_system, {})
    for name, meta in (mblock or {}).items():
        if is_field(meta):
            d = dict(meta)
            # grid-keyed default, e.g. num_helpers {D:..,Q:..}
            if isinstance(d.get("default"), dict):
                d["default"] = d["default"].get(gamera_grid_type, next(iter(d["default"].values())))
            fields.append(_field_from(d, name, f"pbs.{name}", "pbs"))

    for name in _SEGMENT_FIELDS:
        meta = option_descriptions.get("simulation", {}).get(name)
        if is_field(meta):
            fields.append(_field_from(meta, name, f"simulation.{name}", "simulation"))

    for sec in ("gamera", "voltron", "chimp", "dragonking", "raiju"):
        node = option_descriptions.get(sec, {})
        for sub, subnode in node.items():
            if not isinstance(subnode, dict):
                continue
            for name, meta in subnode.items():
                if is_field(meta):
                    fields.append(_field_from(meta, name, f"{sec}.{sub}.{name}", sec))
    return fields


_BCWIND_RANGES = {}


def _bcwind_range(path):
    """(start, stop) dates of a bcwind file, memoized per (path, mtime, size); read errors are memoized too."""
    st = os.stat(path)
    key = (os.path.abspath(path), st.st_mtime_ns, st.st_size)
    if key not in _BCWIND_RANGES:
        import warnings
        import makeitso
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _BCWIND_RANGES[key] = makeitso.fetch_bcwind_time_range(path)
        except Exception as e:
            _BCWIND_RANGES[key] = e
    r = _BCWIND_RANGES[key]
    if isinstance(r, Exception):
        raise r
    return r


def _segmented(v):
    return str(v.get("simulation.use_segments", "N")).upper() == "Y"


def _component(v):
    return "grx" if v.get("__grx") else ("grt" if v.get("__coupling") else "gr")


_HELPERS = {"D": "0", "Q": "2", "O": "4", "H": "8"}
_IJDIR = {"D": "2", "Q": "4", "O": "8", "H": "16"}
_KDIR = {"D": "1", "Q": "1", "O": "1", "H": "1"}


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

    def _bcwind_dates(idx):
        def fn(v):
            if str(v.get("simulation.bcwind_available", "N")).upper() != "Y":
                raise _NotApplicable()
            if v.get("__coupling"):
                raise _NotApplicable()                # coupled: engage's warm-up-shifted dates
            bw = v.get("simulation.bcwind_file", "")
            if not bw or not os.path.isfile(bw):
                raise FileNotFoundError(bw)
            return _bcwind_range(bw)[idx]
        return fn

    def _segment_duration(v):
        sim = _eng(v, "simulation")
        if isinstance(sim, dict) and sim.get("segment_duration") is not None:
            return sim["segment_duration"]
        if _component(v) == "gr" and _segmented(v):
            return makeitso.standalone_segment_default(_hpc(v), v.get("pbs.queue"), v.get("pbs.walltime"),
                                                       _grid(v), _sim_seconds(v))
        return str(_sim_seconds(v))

    def _num_segments(v):
        if str(v.get("simulation.use_segments", "N")).upper() != "Y":
            return "1"
        import math
        return str(math.ceil(_sim_seconds(v) / float(v["simulation.segment_duration"])))

    return {
        "simulation.start_date": (_bcwind_dates(0), ["simulation.bcwind_available", "simulation.bcwind_file"]),
        "simulation.stop_date":  (_bcwind_dates(1), ["simulation.bcwind_available", "simulation.bcwind_file"]),
        "simulation.segment_duration": (_segment_duration,
                                        ["simulation.start_date", "simulation.stop_date",
                                         "simulation.gamera_grid_type", "simulation.hpc_system",
                                         "pbs.queue", "pbs.walltime", "simulation.use_segments"]),
        "pbs.account_name": (lambda v: (_eng(v, "account_name")
                                        or (makeitso.account_prompt_default(_hpc(v)) if makeitso.account_required(_hpc(v)) else None)),
                             ["simulation.hpc_system"]),
        "pbs.kaiju_install_directory": (lambda v: _eng(v, "kaiju_install_directory")
                                        or os.environ.get("KAIJUHOME", ""), []),
        "pbs.kaiju_build_directory": (lambda v: _eng(v, "kaiju_build_directory")
                                      or makeitso.kaiju_build_default(_component(v)), []),
        "pbs.num_segments": (_num_segments, ["simulation.use_segments",
                             "simulation.segment_duration", "simulation.start_date", "simulation.stop_date"]),
        "pbs.num_helpers": (lambda v: _HELPERS.get(_grid(v), "2"), ["simulation.gamera_grid_type"]),
        # cores per node // ranks per node
        "pbs.ompthreads": (lambda v: _threads(v, "pbs.mpiprocs"),
                           ["simulation.hpc_system", "simulation.node_type", "pbs.mpiprocs"]),
        "pbs.helper_ompthreads": (lambda v: _threads(v, "pbs.helper_mpiprocs"),
                                  ["simulation.hpc_system", "simulation.node_type",
                                   "pbs.helper_mpiprocs"]),
        "gamera.sim.H5Grid":  (lambda v: f"lfm{_grid(v)}.h5", ["simulation.gamera_grid_type"]),
        "gamera.sim.runid":   (lambda v: v["simulation.job_name"], ["simulation.job_name"]),
        # a user wind file is linked into the run directory under this name; kaiju cannot read a path here
        "gamera.wind.tsfile": (lambda v: makeitso.CDA2WIND_OUTPUT,
                               ["simulation.bcwind_available", "simulation.bcwind_file"]),
        "gamera.iPdir.N":     (lambda v: _IJDIR.get(_grid(v), "4"), ["simulation.gamera_grid_type"]),
        "gamera.jPdir.N":     (lambda v: _IJDIR.get(_grid(v), "4"), ["simulation.gamera_grid_type"]),
        "gamera.kPdir.N":     (lambda v: _KDIR.get(_grid(v), "1"), ["simulation.gamera_grid_type"]),
        "gamera.coupling.blockHalo": (lambda v: makeitso.machine_science_default(
                                          "gamera.coupling.blockHalo", _hpc(v)), ["simulation.hpc_system"]),
        # tFin/numHelpers are native float/int; options_dict str()s them only when visible
        "voltron.time.tFin": (lambda v: float(_sim_seconds(v)), ["simulation.start_date", "simulation.stop_date"]),
        "voltron.coupling.doAsyncCoupling": (lambda v: makeitso.machine_science_default(
                                                 "voltron.coupling.doAsyncCoupling", _hpc(v)),
                                             ["simulation.hpc_system"]),
        "voltron.helpers.numHelpers": (lambda v: int(v["pbs.num_helpers"]), ["pbs.num_helpers"]),
        "voltron.helpers.useHelpers": (lambda v: "F" if int(v["pbs.num_helpers"]) == 0 else "T",
                                       ["pbs.num_helpers"]),
    }


def _validate_registry():
    def bcwind_file(v, val):
        if str(v.get("simulation.bcwind_available", "N")).upper() != "Y":
            return None
        s = str(val).strip()
        if not (s and os.path.isfile(s)):
            return f"bcwind file {val!r} does not exist"
        try:
            _bcwind_range(s)
        except Exception as e:
            return f"bcwind file {val!r} is not a usable bcwind file: {e}"
        return None

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
            # cda2wind refuses a window under 2 hours
            if str(v.get("simulation.bcwind_available", "N")).upper() != "Y":
                import makeitso
                return makeitso.cda2wind_window_error(v["simulation.start_date"], val)
            return None
        except (ValueError, TypeError, KeyError):
            return None

    def seg_dur(v, val):
        if str(v.get("simulation.use_segments", "N")).upper() != "Y":
            return None
        f = _try_float(val)
        if f is None or f <= 0:
            return "segment_duration must be a positive number"
        if _component(v) != "gr":
            return None
        import makeitso
        return makeitso.standalone_walltime_error(str(v.get("simulation.hpc_system", "")), v.get("pbs.queue"),
                                                  v.get("pbs.walltime"), v.get("simulation.gamera_grid_type"), f)

    def walltime(v, val):
        if _component(v) != "gr":
            return None
        import makeitso
        return makeitso.standalone_walltime_error(str(v.get("simulation.hpc_system", "")), v.get("pbs.queue"),
                                                  val, v.get("simulation.gamera_grid_type"))

    def dtcouple(v, val):
        if not v.get("__coupling"):
            return None
        import makeitso
        return makeitso.dtcouple_prompt_check(val)

    def dospin(v, val):
        if _component(v) != "grt":
            return None
        import makeitso
        return makeitso.coupled_dospin_prompt_check(str(val).strip())

    def account(v, val):
        import makeitso
        return makeitso.account_prompt_check(str(v.get("simulation.hpc_system", "")))(val)

    return {
        "pbs.account_name": account,
        "simulation.bcwind_file": bcwind_file,
        "simulation.start_date": iso,
        "simulation.stop_date": stop_after,
        "simulation.segment_duration": seg_dur,
        "pbs.walltime": walltime,
        "voltron.coupling.dtCouple": dtcouple,
        "voltron.spinup.doSpin": dospin,
    }


def _cond_registry():
    def cond(f, v):
        p = f.path
        bc = str(v.get("simulation.bcwind_available", "N")).upper() == "Y"
        if p == "simulation.bcwind_file":
            return bc
        if p in ("simulation.start_date", "simulation.stop_date"):
            return not bc
        # always the bare 'bcwind.h5', set by its derive
        if p == "gamera.wind.tsfile":
            return False
        if p == "simulation.segment_duration":
            return str(v.get("simulation.use_segments", "N")).upper() == "Y"
        if p == "voltron.coupling.doGCM":
            return False                              # never asked (engage sets it per deck)
        # under engage the account is engage's
        if p == "pbs.account_name" and not v.get("__engage"):
            import makeitso
            return makeitso.account_required(str(v.get("simulation.hpc_system", "")))
        return True
    return cond


def _warn_registry():
    def bcwind_warn(v):
        if str(v.get("simulation.bcwind_available", "N")).upper() != "Y" or v.get("__coupling"):
            return None                       # a coupled run keeps engage's dates
        bw = v.get("simulation.bcwind_file", "")
        rng = ""
        try:
            if bw and os.path.isfile(bw):
                a, b = _bcwind_range(bw)
                rng = f" [{a} .. {b}]"
        except Exception:
            pass
        return (f"WARNING: bcwind_available=Y -> run start/stop dates are taken from the bcwind "
                f"file '{bw}'{rng}, OVERRIDING any otherwise-set dates. Use bcwind_available=N "
                f"to keep caller-provided dates.")

    def dtcouple_warn(v):
        import makeitso
        return makeitso._import_segment_solver().dtcouple_warning(v.get("voltron.coupling.dtCouple"))
    return {"simulation.bcwind_file": bcwind_warn, "voltron.coupling.dtCouple": dtcouple_warn}


def build_state(args, mode="BASIC"):
    import makeitso
    option_descriptions = makeitso.load_option_descriptions()
    try:
        makeitso.update_option_descriptions(option_descriptions, args if isinstance(args, dict) else vars(args))
    except Exception:
        pass
    eng = args if isinstance(args, dict) else None
    hpc0 = ((eng or {}).get("simulation", {}).get("hpc_system")
            or option_descriptions["simulation"]["hpc_system"].get("default", "aitken"))
    grid0 = ((eng or {}).get("simulation", {}).get("gamera_grid_type")
             or option_descriptions["simulation"]["gamera_grid_type"].get("default", "Q"))
    # every simulation/pbs field engage supplied is hidden (its value is still emitted)
    provided = set()
    if eng and ("coupling" in eng or "pbs" in eng):
        for sec in ("simulation", "pbs"):
            for k in (eng.get(sec) or {}):
                provided.add(f"{sec}.{k}")
    ctx = {"__engage_cfg": eng, "__engage": bool(eng and ("coupling" in eng or "pbs" in eng)),
           "__coupling": bool(eng and "coupling" in eng), "__grx": bool(eng and eng.get("grx")),
           "__hpc0": hpc0, "__grid0": grid0,
           "__engage_provided": provided}
    state = MakeitsoFormState(option_descriptions, _derive_registry(), _validate_registry(),
                              mode, ctx, _cond_registry(), _warn_registry())
    state.choices = {"pbs.account_name": lambda v: makeitso.account_choices(str(v.get("simulation.hpc_system", "")))}

    def coerce_error(f, value):
        # a relative kaiju directory is read against the run directory typed so far
        if f.meta is None or str(value).strip() == "":
            return None
        rd = state.values.get("pbs.run_directory")
        run_dir = os.path.abspath(os.path.expanduser(str(rd))) if str(rd or "").strip() else None
        return makeitso.option_value_error(f.path, f.meta, value, run_dir)
    state.coerce_error = coerce_error
    _seed_engage(state, eng)
    return state


# Derived here, not taken from engage: num_segments covers the warm-up-extended span.
_ENGAGE_NOT_TRANSFERRED = ("pbs.num_segments",)


def _seed_engage(state, eng):
    """Pin the values engage supplied; a coupled start_date pins the warm-up-shifted default instead."""
    if not isinstance(eng, dict):
        return
    coupling = "coupling" in eng
    for sec in ("simulation", "pbs"):
        block = eng.get(sec)
        if not isinstance(block, dict):
            continue
        for k, val in block.items():
            path = f"{sec}.{k}"
            if path not in state.by_path or path in _ENGAGE_NOT_TRANSFERRED:
                continue
            if coupling and path == "simulation.start_date":
                val = state.by_path[path].default
            if val is None:
                continue
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
    """Run the makeitso wizard and return the options dict."""
    require_tui()
    mode = _arg(args, "mode") or "BASIC"
    return run_form(build_state(args, mode=mode), "makeitso")


def require_tui():
    """Raise _FormUnavailable without a TTY on stdin/stdout or an importable prompt_toolkit."""
    import sys
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        raise _FormUnavailable("no TTY")
    try:
        import prompt_toolkit  # noqa: F401
    except ImportError as e:
        raise _FormUnavailable(str(e))


def run_form(state, title, collect=None):
    """The wizard screen over any FormState, one field at a time.

    Returns collect(state) (default state.options_dict()) once every visible field is valid;
    raises _FormUnavailable when it cannot run and _FormCancelled on Ctrl-C.
    """
    require_tui()
    try:
        from prompt_toolkit import Application
        from prompt_toolkit.application import get_app
        from prompt_toolkit.key_binding import KeyBindings
        from prompt_toolkit.filters import Condition
        from prompt_toolkit.layout import Layout
        from prompt_toolkit.layout.containers import ConditionalContainer, HSplit, VSplit, Window
        from prompt_toolkit.layout.controls import FormattedTextControl
        from prompt_toolkit.layout.dimension import D
        from prompt_toolkit.widgets import TextArea
        from prompt_toolkit.styles import Style
        try:
            from prompt_toolkit.layout import ScrollablePane
        except ImportError:
            ScrollablePane = None
    except ImportError as e:
        raise _FormUnavailable(str(e))

    vis0 = state.visible()
    if not vis0:
        raise _FormUnavailable("no visible fields")
    cur = {"path": vis0[0].path}
    blocked = {"on": False}
    result = {"opts": None}

    def vis_paths():
        return [f.path for f in state.visible()]

    input_area = TextArea(text=str(state.values[cur["path"]]), multiline=False, height=1)
    input_area.buffer.cursor_position = len(input_area.text)

    def load(path):
        cur["path"] = path
        input_area.text = str(state.values[path])
        input_area.buffer.cursor_position = len(input_area.text)

    def commit():
        p = cur["path"]
        typed = input_area.text
        if typed != str(state.values[p]):
            state.set(p, typed)

    def picks(p):
        valids = [str(x) for x in (state.by_path[p].valids or [])]
        if valids:
            return valids
        fn = getattr(state, "choices_for", None)
        return fn(p) if fn else []

    def cycle_valid(delta):
        # sets the value at once so downstream defaults update
        p = cur["path"]
        valids = picks(p)
        if not valids:
            return
        cur_text = input_area.text.strip()
        try:
            i = valids.index(cur_text)
        except ValueError:
            i = -1 if delta > 0 else 0
        i = (i + delta) % len(valids)
        input_area.text = valids[i]
        input_area.buffer.cursor_position = len(input_area.text)
        state.set(p, valids[i])

    def finish():
        commit()
        errs = state.errors()
        if not errs:
            result["opts"] = collect(state) if collect else state.options_dict()
            app.exit()
        else:
            for p in vis_paths():
                if p in errs:
                    load(p)
                    return

    def move(delta):
        commit()
        if delta > 0 and cur["path"] in state.errors():
            blocked["on"] = True
            return
        blocked["on"] = False
        vp = vis_paths()
        p = cur["path"]
        if p not in vp:
            load(vp[0])
            return
        j = vp.index(p) + delta
        if j < 0:
            load(vp[0])
        elif j >= len(vp):
            finish()
        else:
            load(vp[j])

    def tag(p):
        if p in state.overridden:
            return "your value"
        return "derived" if p in state.derive else "default"

    def header():
        vp = vis_paths()
        p = cur["path"]
        i = (vp.index(p) + 1) if p in vp else 0
        f = state.by_path[p]
        return [("class:title", f" {title} wizard — {state.mode}    step {i}/{len(vp)}    "
                                f"[{f.section}] ")]

    def _shown(raw):
        return "—" if (raw is None or str(raw).strip() in ("", "[None]", "[none]", "None",
                                                           "[null]")) else str(raw)

    def overview():
        vp = vis_paths()
        p = cur["path"]
        ci = vp.index(p) if p in vp else 0
        lo, hi = max(0, ci - 9), min(len(vp), ci + 3)
        out = []
        if lo > 0:
            out.append(("class:more", f"     … {lo} earlier\n"))
        for i in range(lo, hi):
            f = state.by_path[vp[i]]
            val = _shown(state.values[vp[i]])
            if len(val) > 48:
                val = val[:47] + "…"
            if i == ci:
                out.append(("class:cur", f"  ▶ {f.name:<22} {val}"))
                out.append(("class:cur", "\n"))
            else:
                mark = "✓" if vp[i] in state.overridden else " "
                out.append(("class:donemark", f"  {mark} "))
                out.append(("class:donename", f"{f.name:<22} "))
                out.append(("class:doneval", f"{val}\n"))
        if hi < len(vp):
            out.append(("class:more", f"     … {len(vp) - hi} more\n"))
        return out

    def body():
        p = cur["path"]
        f = state.by_path[p]
        out = [("class:name", f"  {f.name}"), ("class:leveltag", f"   [{f.level}]\n")]
        for ln in str(f.prompt).split("\n"):
            out.append(("class:prompttext", f"  {ln}\n"))
        if f.description:
            out.append(("class:help", f"  {f.description}\n"))
        shown = picks(p)
        if shown:
            curv = str(state.values[p]).strip()
            chips = "  ".join((f"[{v}]" if str(v) == curv else f" {v} ") for v in shown)
            how = "(←/→ to choose)" if f.valids else "(←/→ to choose, or type another)"
            out.append(("class:valids", f"  {chips}   {how}\n"))
        out.append(("class:dim", f"  default ({tag(p)}):  "))
        out.append(("class:default", f"{_shown(state.values[p])}\n"))
        if f.warning:
            out.append(("class:warn", f"  ⚠ {f.warning}\n"))
        dw = state.warnings_for(p)
        if dw and dw != f.warning:
            for ln in str(dw).split("\n"):
                out.append(("class:warn", f"  ⚠ {ln}\n"))
        if p in state.derive_errors:
            out.append(("class:warn", f"  (could not derive — check upstream: "
                                      f"{state.derive_errors[p]})\n"))
        return out

    def errbox():
        # its own window above the input row, so a long help text cannot clip the error
        p = cur["path"]
        out = []
        err = state.errors().get(p)
        if err:
            out.append(("class:err", f"  ✗ {err}\n"))
            if blocked["on"]:
                out.append(("class:err", "  ↳ fix this to continue (or ↑ to go back)\n"))
        if out:
            out[-1] = (out[-1][0], out[-1][1].rstrip("\n"))
        return out

    def footer():
        # status first: a narrow terminal truncates the end of the line
        n = len(state.errors())
        try:
            cols = get_app().output.get_size().columns
        except Exception:
            cols = 120
        wide = cols >= 100
        status = (("class:err", f"  ({n} field(s) need fixing)") if n
                  else ("class:ok", "  (all valid — Ctrl-S to finish)" if wide
                        else "  (all valid)"))
        legend = ("   ↑ back   ↓/Enter next   Ctrl-S finish   Ctrl-R reset   Ctrl-C cancel"
                  if wide else "  ↑ back  ↓/Enter next  ^S finish  ^R reset  ^C cancel")
        return [status, ("class:foot", legend)]

    left = HSplit([
        Window(FormattedTextControl(body), dont_extend_height=True, wrap_lines=True),
        ConditionalContainer(
            Window(FormattedTextControl(errbox), dont_extend_height=True, wrap_lines=True),
            filter=Condition(lambda: bool(errbox()))),
        VSplit([Window(FormattedTextControl([("class:caret", "  › ")]), width=4), input_area]),
        Window(),
    ], width=D(weight=3))
    overview_win = Window(FormattedTextControl(overview), dont_extend_height=True)
    right = (ScrollablePane(overview_win, width=D(weight=2)) if ScrollablePane
             else HSplit([overview_win], width=D(weight=2)))
    root = HSplit([
        Window(FormattedTextControl(header), height=1, style="class:title"),
        VSplit([left, Window(width=1, char="│", style="class:rule"), right]),
        Window(FormattedTextControl(footer), height=1, style="class:title"),
    ])

    kb = KeyBindings()

    @kb.add("up")
    def _(e):
        move(-1)

    @kb.add("down")
    @kb.add("enter")
    def _(e):
        move(1)

    # only on fields with picks; free-text fields keep left/right cursor movement
    _on_valids = Condition(lambda: bool(picks(cur["path"])))

    @kb.add("left", filter=_on_valids)
    def _(e):
        cycle_valid(-1)

    @kb.add("right", filter=_on_valids)
    def _(e):
        cycle_valid(1)

    @kb.add("c-s")
    def _(e):
        finish()

    @kb.add("c-r")
    def _(e):
        p = cur["path"]
        state.reset(p)
        blocked["on"] = False
        load(p)

    @kb.add("c-c")
    def _(e):
        e.app.exit()

    style = Style.from_dict({
        "title": "reverse", "cur": "reverse bold", "donemark": "ansigreen",
        "donename": "ansibrightblack", "doneval": "ansicyan", "more": "ansibrightblack italic",
        "rule": "ansibrightblack", "name": "bold ansicyan", "leveltag": "ansibrightblack",
        "prompttext": "", "help": "ansibrightblack", "valids": "ansimagenta",
        "dim": "ansibrightblack", "default": "ansigreen bold", "warn": "ansiyellow",
        "err": "ansired bold", "caret": "ansigreen bold", "foot": "ansibrightblack", "ok": "ansigreen",
    })
    app = Application(layout=Layout(root), key_bindings=kb, style=style, full_screen=True)
    app.run()
    if result["opts"] is None:
        raise _FormCancelled("cancelled")
    return result["opts"]
