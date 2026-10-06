"""Check that the solar, LBC and SST forcing files cover a waxrun window.

CAM does not extrapolate these files, so a window past their last date fails hours into the run.
"""

import os
import re
import sys
import xml.etree.ElementTree as ET
from datetime import date, datetime, timedelta

import netCDF4

FORCING_KEYS = ("solar_irrad_data_file", "solar_parms_data_file", "solar_wind_data_file",
                "flbc_file")
OPTIONAL_KEYS = frozenset(["solar_wind_data_file"])
SOLAR_KEYS = ("solar_irrad_data_file", "solar_parms_data_file", "solar_wind_data_file")
SHORT = {"solar_irrad_data_file": "solar_irrad", "solar_parms_data_file": "solar_parms",
         "solar_wind_data_file": "solar_wind", "flbc_file": "flbc"}

_NL_LINE = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=\s*)(['\"])([^'\"]*)\4(.*)$")


def _warn(s):
    print(s, file=sys.stderr)


def _dt(v):
    if isinstance(v, datetime):
        return v
    if isinstance(v, date):
        return datetime(v.year, v.month, v.day)
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%m-%d-%Y %H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            continue
    raise ValueError(f"unrecognized datetime: {v!r}")


def namelist_files(text, keys=FORCING_KEYS):
    """{key: path or None} from a namelist text; the last assignment wins."""
    out = {k: None for k in keys}
    for line in text.splitlines():
        if line.lstrip().startswith("!"):
            continue
        m = _NL_LINE.match(line)
        if m and m.group(2) in out:
            out[m.group(2)] = m.group(5)
    return out


def _ymd(d, sec=0):
    d = int(d)
    return datetime(d // 10000, (d // 100) % 100, d % 100) + timedelta(seconds=int(sec))


def file_span(path):
    """(first, last, n) of a forcing file's time axis (date/datesec, else time)."""
    with netCDF4.Dataset(str(path)) as ds:
        ds.set_auto_mask(False)
        if "date" in ds.variables:
            d = ds.variables["date"]
            n = d.shape[0]
            if "datesec" in ds.variables:
                s = ds.variables["datesec"]
                return _ymd(d[0], s[0]), _ymd(d[-1], s[-1]), n
            return _ymd(d[0]), _ymd(d[-1]), n
        if "time" in ds.variables:
            t = ds.variables["time"]
            cal = getattr(t, "calendar", "standard")
            first, last = netCDF4.num2date([t[0], t[-1]], t.units, cal,
                                           only_use_cftime_datetimes=False)
            return _dt(first.strftime("%Y-%m-%d %H:%M:%S")), \
                _dt(last.strftime("%Y-%m-%d %H:%M:%S")), t.shape[0]
    raise ValueError(f"{path}: no date or time axis")


def _xml_entries(path, ids):
    root = ET.parse(str(path)).getroot()
    out = {}
    for e in root.iter("entry"):
        if e.get("id") in ids:
            out[e.get("id")] = e.get("value")
    return out


def sst_settings(env_run, docn_streams=None):
    """The SST file, years and taxmode from env_run.xml and docn.streams.xml (taxmode None if absent)."""
    ent = _xml_entries(env_run, {"SSTICE_DATA_FILENAME", "SSTICE_YEAR_START", "SSTICE_YEAR_END",
                                 "SSTICE_YEAR_ALIGN", "DOCN_MODE"})
    sst = {"file": ent.get("SSTICE_DATA_FILENAME"),
           "year_first": int(ent["SSTICE_YEAR_START"]) if ent.get("SSTICE_YEAR_START") else None,
           "year_last": int(ent["SSTICE_YEAR_END"]) if ent.get("SSTICE_YEAR_END") else None,
           "year_align": int(ent["SSTICE_YEAR_ALIGN"]) if ent.get("SSTICE_YEAR_ALIGN") else None,
           "taxmode": None, "docn_mode": ent.get("DOCN_MODE")}
    if docn_streams and os.path.isfile(str(docn_streams)):
        root = ET.parse(str(docn_streams)).getroot()
        for si in root.iter("stream_info"):
            if si.get("name") != "prescribed":
                continue
            for tag, key, cast in (("taxmode", "taxmode", str), ("year_first", "year_first", int),
                                   ("year_last", "year_last", int),
                                   ("year_align", "year_align", int)):
                el = si.find(tag)
                if el is not None and el.text and el.text.strip():
                    sst[key] = cast(el.text.strip())
            f = si.find("datafiles/file")
            if f is not None and f.text:
                sst["file"] = f.text.strip()
    return sst


def cycled_year(year, year_first, year_last, year_align=None):
    """The data year CDEPS taxmode=cycle actually reads for model `year`."""
    align = year_first if year_align is None else year_align
    return year_first + (year - align) % (year_last - year_first + 1)


def coverage(caseroot, window, atm_in=None, env_run=None, docn_streams=None, strict=True):
    """Report whether every forcing file covers window = (run_startdate, stop + 1 d).

    An uncovered file is an error (SystemExit when strict); cycled SST years are a warning."""
    if caseroot is not None:
        caseroot = str(caseroot)
        atm_in = atm_in or os.path.join(caseroot, "CaseDocs", "atm_in")
        env_run = env_run or os.path.join(caseroot, "env_run.xml")
        docn_streams = docn_streams or os.path.join(caseroot, "CaseDocs", "docn.streams.xml")
    if not atm_in or not env_run:
        raise ValueError("coverage needs a caseroot or explicit atm_in and env_run paths")
    start, end = (_dt(window[0]), _dt(window[1]))
    if end <= start:
        raise ValueError(f"forcing window end {end} is not after start {start}")

    with open(atm_in) as f:
        files = namelist_files(f.read())
    report = {"ok": True, "window": [start.strftime("%Y-%m-%d %H:%M:%S"),
                                     end.strftime("%Y-%m-%d %H:%M:%S")],
              "findings": []}
    for key in FORCING_KEYS:
        short, path = SHORT[key], files[key]
        if path is None:
            level = "info" if key in OPTIONAL_KEYS else "error"
            report[short] = None
            report["findings"].append({"level": level, "key": key, "file": None,
                                       "message": f"{key} is not set in {atm_in}"})
            continue
        if not os.path.isfile(path):
            report[short] = {"file": path, "first": None, "last": None, "n": 0}
            report["findings"].append({"level": "error", "key": key, "file": path,
                                       "message": f"{key} = {path} does not exist"})
            continue
        first, last, n = file_span(path)
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

    sst = sst_settings(env_run, docn_streams)
    report["sst"] = sst
    if sst["year_first"] is not None and sst["year_last"] is not None:
        outside = sorted({y for y in range(start.year, end.year + 1)
                          if not sst["year_first"] <= y <= sst["year_last"]})
        if outside and sst["taxmode"] in (None, "cycle"):
            used = ", ".join(f"{y} -> {cycled_year(y, sst['year_first'], sst['year_last'], sst['year_align'])}"
                             for y in outside)
            how = (f"taxmode={sst['taxmode']}" if sst["taxmode"]
                   else "docn.streams.xml not found; CDEPS default taxmode=cycle assumed")
            msg = (f"SST {os.path.basename(sst['file'] or '?')} covers {sst['year_first']}-"
                   f"{sst['year_last']} ({how}); the run uses cycled SST years: {used}")
            report["findings"].append({"level": "warning", "key": "sst", "file": sst["file"],
                                       "message": msg})
            _warn(f"[waxrun] WARNING: {msg}")
        elif outside:
            report["findings"].append({"level": "warning", "key": "sst", "file": sst["file"],
                                       "message": f"SST years {outside} outside "
                                                  f"{sst['year_first']}-{sst['year_last']} with "
                                                  f"taxmode={sst['taxmode']}"})
    else:
        report["findings"].append({"level": "info", "key": "sst", "file": None,
                                   "message": f"no SSTICE_YEAR_START/END in {env_run}"})

    errors = [f for f in report["findings"] if f["level"] == "error"]
    report["ok"] = not errors
    if errors and strict:
        raise SystemExit("forcing coverage gate FAILED for window "
                         f"{start:%Y-%m-%d %H:%M} .. {end:%Y-%m-%d %H:%M}:\n  "
                         + "\n  ".join(e["message"] for e in errors)
                         + "\nCAM does not extrapolate these files: use ones that cover the window.")
    return report


def harmonize(coupled_nl_text, spinup_clone_nl_text):
    """coupled_nl_text with its solar_* files set to those of the spinup clone."""
    src = namelist_files(spinup_clone_nl_text, SOLAR_KEYS)
    want = {k: v for k, v in src.items() if v}
    if not want:
        return coupled_nl_text
    lines = coupled_nl_text.splitlines()
    done, last_solar = set(), None
    for i, line in enumerate(lines):
        if line.lstrip().startswith("!"):
            continue
        m = _NL_LINE.match(line)
        if not m or m.group(2) not in want:
            continue
        lead, key, eq, q, _, tail = m.groups()
        lines[i] = f"{lead}{key}{eq}{q}{want[key]}{q}{tail}"
        done.add(key)
        last_solar = i
    missing = [k for k in SOLAR_KEYS if k in want and k not in done]
    if missing:
        template = lines[last_solar] if last_solar is not None else None
        m = _NL_LINE.match(template) if template else None
        lead, width, eq, q = (m.group(1), len(m.group(2)), m.group(3), m.group(4)) if m else (" ", 0, " = ", "'")
        new = [f"{lead}{k}{' ' * max(0, width - len(k))}{eq}{q}{want[k]}{q}" for k in missing]
        at = last_solar + 1 if last_solar is not None else len(lines)
        lines[at:at] = new
    out = "\n".join(lines)
    if coupled_nl_text.endswith("\n"):
        out += "\n"
    return out
