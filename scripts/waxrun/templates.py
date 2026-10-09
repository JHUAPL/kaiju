"""The waxrun user_nl_cam templates.

Rendering a template with the reference values reproduces that case's user_nl_cam byte-for-byte.
"""

from pathlib import Path

from jinja2 import Environment, StrictUndefined, UndefinedError, meta

HERE = Path(__file__).resolve().parent
CONFIG_DIR = HERE / "config"

TEMPLATES = {
    "spinup": CONFIG_DIR / "user_nl_cam.spinup.j2",
    "coupled": CONFIG_DIR / "user_nl_cam.coupled.j2",
}

# Also used by the coupled case: the coupled template's own files end in 2024.
REFERENCE_SOLAR = {
    "spinup": {
        "solar_irrad_data_file": "/glade/campaign/cesm/development/wawg/joemci/inputdata/atm/waccm/solar/SolarForcingNNLSSI1-ssi_v03r00-preliminary_daily_s18740509_e20260331_260716.nc",
        "solar_parms_data_file": "/glade/campaign/cesm/development/wawg/joemci/inputdata/atm/waccm/solar/wax_solar_3hr_1960-Feb2026_c260504.nc",
        "solar_wind_data_file": "/glade/campaign/cesm/development/wawg/joemci/inputdata/solar/solar_wind_imf_OMNI_WACCMX_2000001-2026212_c260814.nc",
    },
}

_NUDGE_WINDOW = ["nudge_beg_year", "nudge_beg_month", "nudge_beg_day",
                 "nudge_end_year", "nudge_end_month", "nudge_end_day"]
_COMMON = ["ncdata", "inithist", "nudge_model", "nudge_path", "model_times_per_day", *_NUDGE_WINDOW,
           "nudge_vwin_lindex", "nudge_vwin_hindex", "solar_irrad_data_file", "solar_parms_data_file"]

# render() refuses anything else
VARIABLES = {
    "spinup": ["ncdata", "nudge_t", "inithist", "solar_irrad_data_file", "solar_parms_data_file",
               "nudge_model", "nudge_path", "model_times_per_day", *_NUDGE_WINDOW,
               "nudge_vwin_lindex", "nudge_vwin_hindex"],
    "coupled": ["ncdata", "solar_wind_data_file", "inithist", "dtime", "solar_irrad_data_file",
                "solar_parms_data_file", "nudge_model", "nudge_path", "model_times_per_day", *_NUDGE_WINDOW,
                "nudge_vwin_lindex", "nudge_vwin_hindex"],
}
# taken from the spec, not spec['nl_vars']
SPEC_VARIABLES = {"spinup": [], "coupled": ["dtime"]}

_INT = ["model_times_per_day", *_NUDGE_WINDOW, "nudge_vwin_lindex", "nudge_vwin_hindex", "dtime"]


def _kind(kind):
    if kind not in TEMPLATES:
        raise ValueError(f"unknown template kind {kind!r}: {sorted(TEMPLATES)}")
    return kind


def environment():
    return Environment(undefined=StrictUndefined, keep_trailing_newline=True,
                       trim_blocks=False, lstrip_blocks=False, autoescape=False)


def template_text(kind):
    return TEMPLATES[_kind(kind)].read_text()


def variables(kind):
    """The variable names the template file references."""
    return meta.find_undeclared_variables(environment().parse(template_text(kind)))


def _check(kind, vars):
    unknown = sorted(set(vars) - set(VARIABLES[kind]))
    if unknown:
        raise ValueError(f"{kind} template: unknown variable(s) {unknown}; it takes {VARIABLES[kind]}")
    for k in _INT:
        if k in vars and (isinstance(vars[k], bool) or not isinstance(vars[k], int)):
            raise ValueError(f"{kind} template: {k} must be an int, got {vars[k]!r}")
    if "nudge_model" in vars and vars["nudge_model"] not in (".true.", ".false."):
        raise ValueError(f"{kind} template: nudge_model must be '.true.' or '.false.', got {vars['nudge_model']!r}")
    if "nudge_t" in vars and not isinstance(vars["nudge_t"], bool):
        raise ValueError(f"{kind} template: nudge_t must be a bool, got {vars['nudge_t']!r}")
    if "nudge_path" in vars and not str(vars["nudge_path"]).endswith("/"):
        raise ValueError(f"{kind} template: nudge_path must end with '/' (CAM concatenates "
                         f"Nudge_Path//Nudge_File), got {vars['nudge_path']!r}")
    for k, v in vars.items():
        if k not in _INT and k != "nudge_t" and not isinstance(v, str):
            raise ValueError(f"{kind} template: {k} must be a str, got {v!r}")
    if "dtime" in vars:
        dt = vars["dtime"]
        if dt <= 0 or 60 % dt:
            raise ValueError(f"{kind} template: dtime={dt} s must divide 60 (h1 nhtfrq = 60 // dtime)")
        if "model_times_per_day" in vars and vars["model_times_per_day"] != 86400 // dt:
            raise ValueError(f"{kind} template: model_times_per_day={vars['model_times_per_day']} "
                             f"contradicts dtime={dt} s (86400 // dtime = {86400 // dt})")


def render(kind, vars):
    """The user_nl_cam text for kind ('spinup' | 'coupled') from exactly its VARIABLES.

    Paths are passed unquoted; solar_wind_data_file '' drops that line."""
    _check(_kind(kind), vars)
    try:
        return environment().from_string(template_text(kind)).render(**vars)
    except UndefinedError as e:
        raise UndefinedError(f"{kind} template: {e} (it takes {VARIABLES[kind]})") from None
