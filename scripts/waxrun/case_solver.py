"""CIME case plumbing for waxrun: create_newcase, xmlchange, namelist patches, modules and build.

Runners take dry_run to record commands instead of executing them; nothing here submits a job.
"""

import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG_DIR = HERE / "config"
MODULE_SYSTEM_XML = CONFIG_DIR / "module_system.derecho_2309.xml"
SOURCEMODS_DIR = CONFIG_DIR / "sourcemods"

SPLICE_GUARD = ".waxrun_spliced"
BUILD_LOG = "build_2309.log"
CAM_CONFIG_CACHE = os.path.join("Buildconf", "camconf", "config_cache.xml")

DEFAULT_COMPSET = "FXHIST"
DEFAULT_RES = "ne30pg3_ne30pg3_mg17"
DEFAULT_MACHINE = "derecho"

# -nlev 145 is required by the L145 IC and nudging data
CAM_CONFIG_OPTS_APPEND = "-nlev 145 -pcols 2"
CPPDEF_OF_OPT = {"nlev": "PLEV", "pcols": "PCOLS"}

# integer literals: xmlchange rejects '$ATM_NCPL' for an int entry. WAV is a stub.
NCPL_COMPCLASSES = ("ATM", "LND", "ICE", "OCN", "GLC", "ROF")
NTASKS_COMPCLASSES = ("ATM", "LND", "ICE", "OCN", "ROF", "CPL")
# the tree's config_pes.xml default: 16 nodes x 128 tasks
TREE_PE_NODES = 16
TREE_NTASKS_ATM = 2048

# lookup order: the first file wins for a duplicated id
ENV_FILES = ["env_run.xml", "env_build.xml", "env_case.xml", "env_mach_pes.xml",
             "env_workflow.xml", "env_batch.xml", "env_archive.xml", "env_postprocessing.xml"]

# case values matched against the <modules> conditions of env_mach_specific.xml
DEFAULT_CASE_SETTINGS = {"COMPILER": "intel", "MPILIB": "mpich", "DEBUG": "FALSE",
                         "GPU_TYPE": "none", "COMP_INTERFACE": "nuopc"}


def _warn(s):
    print(s, file=sys.stderr)


def _run(cmd, **kw):
    print("+ " + (cmd if isinstance(cmd, str) else " ".join(map(str, cmd))))
    return subprocess.run(cmd, check=True, **kw)


def finding(level, check, detail):
    """One finding; level is 'OK', 'WARN' or 'FAIL'."""
    return {"level": level, "check": check, "detail": detail}


def has_fail(findings):
    return any(f["level"] == "FAIL" for f in findings)


def create_newcase_cmd(spec):
    """argv for create_newcase from the spec.

    --run-unsupported is needed because FXHIST on ne30pg3 is not a supported pair."""
    cmd = [os.path.join(spec["srcroot"], "cime", "scripts", "create_newcase")]
    if spec.get("run_unsupported", True):
        cmd.append("--run-unsupported")
    cmd += ["--compset", str(spec.get("compset") or DEFAULT_COMPSET),
            "--res", str(spec.get("res") or DEFAULT_RES),
            "--mach", str(spec.get("machine") or spec.get("cime_machine") or DEFAULT_MACHINE),
            "--compiler", str(spec.get("compiler") or DEFAULT_CASE_SETTINGS["COMPILER"]),
            "--mpilib", str(spec.get("mpilib") or DEFAULT_CASE_SETTINGS["MPILIB"]),
            "--project", str(spec["project"]),
            "--case", str(spec["caseroot"]),
            "--output-root", str(spec["cime_output_root"]),
            "--handle-preexisting-dirs", "u"]
    return cmd


def xmlchange_command(item):
    """argv for one xmlchange item: KEY=VALUE strings, optionally led by --subgroup <g> or --append."""
    item = list(item)
    argv = ["./xmlchange"]
    if item and item[0] == "--subgroup":
        argv += item[:2]
        item = item[2:]
    elif item and item[0] == "--append":
        argv.append(item[0])
        item = item[1:]
        if len(item) != 1:
            raise ValueError(f"--append takes exactly one KEY=VALUE, got {item!r}")
    if not item:
        raise ValueError("empty xmlchange item")
    for kv in item:
        if "=" not in kv:
            raise ValueError(f"xmlchange item is not KEY=VALUE: {kv!r}")
    argv.append(",".join(item))
    return argv


_KV_TOKEN = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$", re.S)


def _quote_token(tok):
    """shlex.quote, but KEY=VALUE needing quotes becomes KEY="VALUE"."""
    if shlex.quote(tok) == tok:
        return tok
    m = _KV_TOKEN.match(tok)
    if m and not set('"$`\\\n') & set(m.group(2)):
        return f'{m.group(1)}="{m.group(2)}"'
    return shlex.quote(tok)


def shell_line(argv):
    return " ".join(_quote_token(str(a)) for a in argv)


def run_xmlchanges(caseroot, items, dry_run=False, log=None):
    """Run the xmlchange items in caseroot and return their shell lines (appended to log if given)."""
    lines = []
    for item in items:
        argv = xmlchange_command(item)
        lines.append(shell_line(argv))
        if not dry_run:
            _run(argv, cwd=str(caseroot))
    if log:
        with open(log, "a") as fh:
            fh.write("".join(l + "\n" for l in lines))
    return lines


def ncpl_from_dtime(dtime):
    """Coupling intervals per hour for a physics timestep in seconds."""
    dtime = int(dtime)
    if dtime <= 0 or 3600 % dtime or 86400 % dtime:
        raise ValueError(f"dtime={dtime} s must divide 3600 s and 86400 s")
    return 3600 // dtime


def config_xmlchanges(kind, spec):
    """The ordered pre-setup xmlchange items of a new case.

    CAM_CONFIG_OPTS must come first: a preview_namelists before it leaves stale CAM_CPPDEFS."""
    if kind not in ("spinup", "coupled"):
        raise ValueError(f"kind must be 'spinup' or 'coupled', not {kind!r}")
    cfg = spec.get("config_xml") or {}
    items = [["--append", "CAM_CONFIG_OPTS="
              + (cfg.get("cam_config_opts_append") or CAM_CONFIG_OPTS_APPEND)]]
    if kind == "spinup" and cfg.get("pio_netcdf_format"):
        items.append([f"PIO_NETCDF_FORMAT={cfg['pio_netcdf_format']}"])
    if not cfg.get("atm_ncpl") and "dtime" not in spec:
        raise KeyError(f"the {kind} spec needs 'dtime' (or config_xml['atm_ncpl']) for the NCPL "
                       f"xmlchanges")
    ncpl = int(cfg["atm_ncpl"]) if cfg.get("atm_ncpl") else ncpl_from_dtime(spec["dtime"])
    if cfg.get("rof_ncpl") and int(cfg["rof_ncpl"]) != ncpl:
        raise ValueError(f"config_xml rof_ncpl={cfg['rof_ncpl']} != atm_ncpl={ncpl}: MOSART joins "
                         f"the fast loop at the ATM cadence in both phases")
    items.append([f"NCPL_BASE_PERIOD={cfg.get('ncpl_base_period') or 'hour'}"])
    items.append([f"{c}_NCPL={ncpl}" for c in NCPL_COMPCLASSES])
    # explicit None checks so pe_nodes 0 is refused, not defaulted
    nodes = cfg.get("pe_nodes")
    if nodes is None:
        nodes = (spec.get("pe") or {}).get("nodes")
    nodes = int(TREE_PE_NODES if nodes is None else nodes)
    if nodes < 1:
        raise ValueError(f"{kind}.pe_nodes={nodes} must be >= 1 node")
    if nodes != TREE_PE_NODES:
        _warn(f"[waxrun] WARNING: {kind}.pe_nodes={nodes} != the tree default {TREE_PE_NODES}: "
              f"NTASKS_* are overridden and the measured wall_s_per_simday no longer applies")
        items.append([f"NTASKS_{c}=-{nodes}" for c in NTASKS_COMPCLASSES])
    return items


def spinup_xmlchanges(spec):
    """The spinup xmlchange items as {'pre_setup': [...], 'post_setup': [...]}."""
    account = spec["project"]
    charge = spec.get("charge_account") or account
    # a chunked spinup configures the first chunk; the POSTRUN_SCRIPT sets the last one's STOP_N.
    # CMEPS writes a restart at the end of every run regardless of REST_N.
    chunks = spec.get("chunks") or []
    stop_n = chunks[0]["stop_n"] if len(chunks) > 1 else spec["stop_n"]
    pre = config_xmlchanges("spinup", spec) + [
        ["RUN_TYPE=startup"],
        [f"RUN_STARTDATE={spec['run_startdate']}", f"START_TOD={spec['start_tod']}"],
        [f"STOP_OPTION={spec['stop_option']}", f"STOP_N={stop_n}", "RESUBMIT=0"],
        ["REST_OPTION=ndays", "REST_N=5"],
        ["DOUT_S=FALSE"],
        [f"PROJECT={account}", f"CHARGE_ACCOUNT={charge}"],
    ]
    if len(chunks) > 1:
        pre.append([f"POSTRUN_SCRIPT={os.path.join(spec['caseroot'], 'waxrun_postrun.sh')}"])
    if (spec.get("nudging") or {}).get("enabled"):
        pre.append([f"PRERUN_SCRIPT={os.path.join(spec['caseroot'], 'waxrun_prerun.sh')}"])
    post = [["--subgroup", "case.run", f"JOB_QUEUE={spec['queue']}",
             f"JOB_WALLCLOCK_TIME={spec['walltime_est']}"]]
    if spec.get("job_priority"):
        # case-global: CIME aborts on JOB_PRIORITY inside a --subgroup item
        post.append([f"JOB_PRIORITY={spec['job_priority']}"])
    post.append(["--subgroup", "case.st_archive", "JOB_WALLCLOCK_TIME=00:20:00"])
    return {"pre_setup": pre, "post_setup": post}


def coupled_xmlchanges(spec):
    """The coupled xmlchange items; no batch settings, since the coupled case is never case.submit-ed."""
    account = spec["project"]
    charge = spec.get("charge_account") or account
    pre = config_xmlchanges("coupled", spec) + [
        ["RUN_TYPE=startup"],
        [f"RUN_STARTDATE={spec['run_startdate']}", f"START_TOD={spec['start_tod']}"],
        ["STOP_OPTION=nhours", f"STOP_N={spec['stephours']}", "RESUBMIT=0", "CONTINUE_RUN=FALSE"],
        # literal value: xmlchange rejects '$STOP_N' for an int entry
        ["REST_OPTION=nhours", f"REST_N={spec['stephours']}"],
        ["DOUT_S=FALSE"],
        [f"PROJECT={account}", f"CHARGE_ACCOUNT={charge}"],
        [f"RUNDIR={spec['rundir']}"],
    ]
    return {"pre_setup": pre, "post_setup": []}


_REF = re.compile(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?")


# negative values mean nodes, so the raw text reads back wrong; use pe_layout for these
_PE_ENTRIES = ("NTASKS", "NTHRDS", "ROOTPE", "PSTRID")


def _compclass_values(entry):
    """{'<id>_<COMPCLASS>': value} for a per-component entry."""
    eid = entry.get("id")
    if eid in _PE_ENTRIES:
        return {}
    out = {}
    for value in entry.iter("value"):
        cc = value.get("compclass")
        if cc:
            out[f"{eid}_{cc.upper()}"] = (value.text or "").strip()
    return out


def _env_values(caseroot):
    """({id: value}, {group id: {id: value}}) over the caseroot's env_*.xml."""
    top, groups = {}, {}
    for name in ENV_FILES:
        path = os.path.join(caseroot, name)
        if not os.path.exists(path):
            continue
        root = ET.parse(path).getroot()
        for group in root.iter("group"):
            gid = group.get("id")
            gmap = groups.setdefault(gid, {})
            for entry in group.iter("entry"):
                eid, val = entry.get("id"), entry.get("value")
                gmap.setdefault(eid, val)
                top.setdefault(eid, val)
                for cid, cval in _compclass_values(entry).items():
                    gmap.setdefault(cid, cval)
                    top.setdefault(cid, cval)
        for entry in root.iter("entry"):
            top.setdefault(entry.get("id"), entry.get("value"))
            for cid, cval in _compclass_values(entry).items():
                top.setdefault(cid, cval)
    return top, groups


def _resolve(value, table, depth=0):
    if value is None or "$" not in value or depth > 10:
        return value

    def sub(m):
        ref = table.get(m.group(1))
        return m.group(0) if ref is None else ref
    new = _REF.sub(sub, value)
    return new if new == value else _resolve(new, table, depth + 1)


def xmlquery(caseroot, names, subgroup=None, resolve=True):
    """{name: value or None} read from the caseroot's env_*.xml without importing CIME."""
    top, groups = _env_values(caseroot)
    table = dict(top)
    if subgroup:
        table.update(groups.get(subgroup, {}))
    out = {}
    for n in ([names] if isinstance(names, str) else names):
        v = table.get(n)
        out[n] = _resolve(v, table) if resolve else v
    return out


def dtime_from_case(caseroot, expect=None):
    """CAM physics timestep in seconds (3600 / ATM_NCPL); raises if it differs from expect."""
    q = xmlquery(caseroot, ["ATM_NCPL", "NCPL_BASE_PERIOD"])
    if q["NCPL_BASE_PERIOD"] != "hour":
        raise ValueError(f"{caseroot}: NCPL_BASE_PERIOD={q['NCPL_BASE_PERIOD']!r}, expected 'hour'")
    ncpl = int(q["ATM_NCPL"])
    if ncpl <= 0 or 3600 % ncpl:
        raise ValueError(f"{caseroot}: ATM_NCPL={ncpl} does not divide 3600 s")
    dtime = 3600 // ncpl
    if 86400 % dtime:
        raise ValueError(f"{caseroot}: dtime={dtime} s does not divide 86400 s")
    if expect is not None and dtime != int(expect):
        raise ValueError(f"{caseroot}: ATM_NCPL={ncpl} -> dtime={dtime} s, the spec says "
                         f"{int(expect)} s (the NCPL xmlchange did not take)")
    return dtime


def pe_components(caseroot):
    """{COMPCLASS: {NTASKS, ROOTPE, PSTRID}} from env_mach_pes.xml; negative values mean nodes."""
    path = os.path.join(caseroot, "env_mach_pes.xml")
    if not os.path.exists(path):
        return {}
    raw = {}
    for e in ET.parse(path).getroot().iter("entry"):
        eid = e.get("id")
        if eid not in ("NTASKS", "ROOTPE", "PSTRID"):
            continue
        per = {v.get("compclass").upper(): (v.text or "").strip()
               for v in e.iter("value") if v.get("compclass")}
        raw[eid] = (e.get("value"), per)
    top_n, per_n = raw.get("NTASKS", (None, {}))
    comps = sorted(set(per_n)) or (["ATM"] if top_n is not None else [])
    out = {}
    for c in comps:
        row = {}
        for eid, dflt in (("NTASKS", None), ("ROOTPE", 0), ("PSTRID", 1)):
            top, per = raw.get(eid, (None, {}))
            v = per.get(c, top)
            row[eid] = int(v) if v not in (None, "") else dflt
        out[c] = row
    return out


def _world_check(caseroot, comps, per_node, ntasks_atm):
    """Raise unless every CESM rank runs ATM (ROOTPE_ATM 0, PSTRID_ATM 1, no component beyond it).

    mage_init joins MPI_COMM_WORLD collectives from the ATM ranks only, so any other rank hangs the job."""
    def tasks(n):
        return -n * per_node if n < 0 else n
    atm = comps.get("ATM")
    if not atm:
        return
    root_atm = tasks(atm["ROOTPE"])
    bad = []
    if root_atm != 0:
        bad.append(f"ROOTPE_ATM={atm['ROOTPE']} (must be 0)")
    if atm["PSTRID"] != 1:
        bad.append(f"PSTRID_ATM={atm['PSTRID']} (must be 1)")
    ends = {c: tasks(r["ROOTPE"]) + (tasks(r["NTASKS"]) - 1) * r["PSTRID"] + 1
            for c, r in comps.items() if r["NTASKS"] is not None}
    world = max(ends.values()) if ends else ntasks_atm
    if world != ntasks_atm:
        beyond = sorted(c for c, e in ends.items() if e > ntasks_atm)
        bad.append(f"the components span {world} ranks (max ROOTPE + NTASKS*PSTRID"
                   + (f", reached by {', '.join(beyond)}" if beyond else "") + f"), NTASKS_ATM={ntasks_atm}")
    if bad:
        raise RuntimeError(
            f"{caseroot}: PE layout puts CESM ranks outside ATM: {'; '.join(bad)}. Keep every "
            f"component inside the ATM ranks (ROOTPE=0, NTASKS <= NTASKS_ATM), or the coupled job hangs.")


def pe_layout(caseroot, expect=None):
    """{ntasks_atm, nodes, tasks_per_node} of a case; raises if it differs from expect."""
    per_node = int(xmlquery(caseroot, "MAX_MPITASKS_PER_NODE")["MAX_MPITASKS_PER_NODE"] or 128)
    comps = pe_components(caseroot)
    ntasks = (comps.get("ATM") or {}).get("NTASKS")
    if ntasks is None:
        raise ValueError(f"{caseroot}: NTASKS for ATM not found in env_mach_pes.xml")
    if ntasks < 0:
        nodes, ntasks = -ntasks, -ntasks * per_node
    else:
        nodes = -(-ntasks // per_node)
    pe = {"ntasks_atm": ntasks, "nodes": nodes, "tasks_per_node": per_node}
    if expect is not None:
        want = expect if isinstance(expect, int) else int(expect["ntasks_atm"])
        if ntasks != want:
            why = ("the tree's config_pes.xml default changed; set pe_nodes explicitly and "
                   "re-measure wall_s_per_simday" if want == TREE_NTASKS_ATM else
                   f"the pe_nodes override's NTASKS_* xmlchange did not take (the case still has "
                   f"{'the tree default' if ntasks == TREE_NTASKS_ATM else 'another layout'})")
            raise RuntimeError(f"{caseroot}: NTASKS_ATM={ntasks} ({nodes} nodes x {per_node}), "
                               f"expected {want}: {why}")
    _world_check(caseroot, comps, per_node, ntasks)
    return pe


def assert_calendar(caseroot, want="GREGORIAN"):
    """Raise unless the case calendar is Gregorian, which waxrun's date arithmetic assumes."""
    cal = (xmlquery(caseroot, "CALENDAR")["CALENDAR"] or "").strip()
    if cal.upper() != want:
        raise RuntimeError(f"{caseroot}: CALENDAR={cal or '(unset)'}, waxrun needs {want}; use a HIST "
                           f"compset (FXHIST) or ./xmlchange CALENDAR={want}")
    return cal


def _append_pairs(opts_append):
    """'-nlev 145 -pcols 2' -> [('nlev', '145'), ('pcols', '2')]."""
    toks = str(opts_append).split()
    if len(toks) % 2 or not all(t.startswith("-") for t in toks[::2]):
        raise ValueError(f"CAM_CONFIG_OPTS append is not '-opt value' pairs: {opts_append!r}")
    return [(o.lstrip("-"), v) for o, v in zip(toks[::2], toks[1::2])]


def clear_config_cache(caseroot, dry_run=False):
    """Remove a stale CAM config_cache.xml; returns its path, or None.

    Neither case.setup --reset nor case.build --clean-all refreshes it after CAM_CONFIG_OPTS changes."""
    path = os.path.join(caseroot, CAM_CONFIG_CACHE)
    if not os.path.exists(path):
        return None
    print(f"[waxrun] {'(dry-run) ' if dry_run else ''}rm {path} (stale CAM config cache)")
    if not dry_run:
        os.remove(path)
    return path


def assert_cam_config(caseroot, opts_append=CAM_CONFIG_OPTS_APPEND):
    """Findings checking CAM_CONFIG_OPTS, CAM_CPPDEFS and the config cache carry opts_append once."""
    pairs = _append_pairs(opts_append)
    q = xmlquery(caseroot, ["CAM_CONFIG_OPTS", "CAM_CPPDEFS"])
    opts = q["CAM_CONFIG_OPTS"] or ""
    findings = []
    for opt, value in pairs:
        n = len(re.findall(rf"-{re.escape(opt)}\s+{re.escape(value)}(?![\w.])", opts))
        if n == 1:
            findings.append(finding("OK", "cam_config_opts", f"CAM_CONFIG_OPTS has -{opt} {value}"))
        else:
            findings.append(finding("FAIL", "cam_config_opts",
                                    f"{caseroot}: CAM_CONFIG_OPTS has '-{opt} {value}' {n} times, "
                                    f"expected once: {opts!r}"))
    cppdefs = q["CAM_CPPDEFS"] or ""
    for opt, value in pairs:
        macro = CPPDEF_OF_OPT.get(opt)
        if macro is None:
            continue
        if re.search(rf"-D{macro}={re.escape(value)}(?![\w.])", cppdefs):
            findings.append(finding("OK", "cam_cppdefs", f"CAM_CPPDEFS has -D{macro}={value}"))
        else:
            findings.append(finding("FAIL", "cam_cppdefs",
                                    f"{caseroot}: CAM_CPPDEFS lacks -D{macro}={value} — a "
                                    f"preview_namelists that ran before the CAM_CONFIG_OPTS append "
                                    f"leaves the stale value; remove {CAM_CONFIG_CACHE} "
                                    f"and preview again"))
    cache = os.path.join(caseroot, CAM_CONFIG_CACHE)
    if not os.path.exists(cache):
        findings.append(finding("WARN", "config_cache", f"{cache} is absent: preview_namelists "
                                                        f"has not run in this case"))
    else:
        try:
            cached = {e.get("id"): e.get("value") for e in ET.parse(cache).getroot().iter("entry")}
        except ET.ParseError as e:
            cached = None
            findings.append(finding("WARN", "config_cache", f"{cache}: {e}"))
        if cached is not None:
            for opt, value in pairs:
                if opt not in cached:
                    continue
                if cached[opt] == value:
                    findings.append(finding("OK", "config_cache", f"{opt}={value} in {cache}"))
                else:
                    findings.append(finding("FAIL", "config_cache",
                                            f"{cache}: {opt}={cached[opt]!r}, expected {value!r} "
                                            f"(stale cache)"))
    for f in findings:
        print(f"[waxrun] CAM {f['level']:4s} {f['check']}: {f['detail'].splitlines()[0]}")
    return findings


def case_settings(caseroot):
    """The case values the <modules> conditions are matched against, with defaults."""
    q = xmlquery(caseroot, list(DEFAULT_CASE_SETTINGS))
    return {k: (q[k] if q[k] is not None else v) for k, v in DEFAULT_CASE_SETTINGS.items()}


_STMT = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=\s*)(.*)$")
_COMMENT_ONLY = re.compile(r"^\s*!")


def _split_comment(s):
    """(value, sep, comment) of a namelist value line; '!' inside quotes is not a comment."""
    q = None
    for i, ch in enumerate(s):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "!":
            value = s[:i].rstrip()
            return value, s[len(value):i], s[i:]
    return s.rstrip(), s[len(s.rstrip()):], ""


def _statements(lines):
    """[(name, first line, last line)] of every statement, continuation lines included."""
    out = []
    i = 0
    while i < len(lines):
        m = _STMT.match(lines[i])
        if m:
            j = i + 1
            while j < len(lines):
                ln = lines[j]
                if not ln.strip() or _COMMENT_ONLY.match(ln) or _STMT.match(ln):
                    break
                j += 1
            out.append((m.group(2), i, j - 1))
            i = j
        else:
            i += 1
    return out


def _find(stmts, key):
    for name, i, j in stmts:
        if name.lower() == key.lower():
            return name, i, j
    return None


def _render_value(value):
    """(value, comment or None) of a raw RHS value or a [value, comment] pair."""
    if isinstance(value, (list, tuple)):
        v, c = value
        return str(v), f"  ! {c}"
    return str(value), None


def patch_user_nl(text, patches, fincl_edits=None, insert_after=None):
    """Replace or insert key = value in a user_nl_* text, changing nothing else.

    patches values are raw RHS text or [value, comment]; fincl_edits is {key: {"remove", "add"}}."""
    lines = text.split("\n")
    insert_after = insert_after or {}
    for key, value in (patches or {}).items():
        new_value, new_comment = _render_value(value)
        hit = _find(_statements(lines), key)
        if hit:
            name, i, j = hit
            m = _STMT.match(lines[i])
            prefix = m.group(1) + m.group(2) + m.group(3)
            if i == j:
                _, sep, comment = _split_comment(m.group(4))
            else:
                sep, comment = "", ""
            if new_comment is not None:
                sep, comment = new_comment[:2], new_comment[2:]
            lines[i:j + 1] = [prefix + new_value + (sep + comment if comment else "")]
        else:
            anchor = insert_after.get(key)
            hit = _find(_statements(lines), anchor) if anchor else None
            new_line = f"{key} = {new_value}" + (new_comment or "")
            if hit:
                _, i, j = hit
                indent = _STMT.match(lines[i]).group(1)
                lines.insert(j + 1, indent + new_line)
            else:
                if lines and lines[-1] != "":
                    lines.append("")
                lines[-1:] = [new_line, ""]
    for key, edit in (fincl_edits or {}).items():
        hit = _find(_statements(lines), key)
        if not hit:
            raise KeyError(f"{key} not in namelist text; cannot edit its list")
        name, i, j = hit
        m = _STMT.match(lines[i])
        prefix = m.group(1) + m.group(2) + m.group(3)
        region = "\n".join([m.group(4)] + lines[i + 1:j + 1])
        for item in edit.get("remove", []):
            region = _list_remove(region, item)
        for item in edit.get("add", []):
            region = _list_add(region, item)
        lines[i:j + 1] = (prefix + region).split("\n")
    return "\n".join(lines)


def _list_remove(region, item):
    pat = re.compile(r"(['\"])" + re.escape(item) + r"\1")
    m = pat.search(region)
    if not m:
        return region
    start, end = m.start(), m.end()
    after = re.match(r"[ \t]*,[ \t]*", region[end:])
    if after:
        end += after.end()
    else:
        before = re.search(r",\s*$", region[:start])
        if before:
            start = before.start()
    region = region[:start] + region[end:]
    return "\n".join(ln for k, ln in enumerate(region.split("\n")) if k == 0 or ln.strip())


def _list_add(region, item):
    value, sep, comment = _split_comment(region)
    if re.fullmatch(r"\s*(?:(['\"])\s*\1)?\s*", value):
        value = f"'{item}'"
    else:
        value = value + f", '{item}'"
    return value + sep + comment


# a SourceMods set: src.<component>/ folders of files, an optional tree.md5 (md5sum output)
_SET_NAME = re.compile(r"^[A-Za-z0-9_-][A-Za-z0-9._-]*$")
_COMP_DIR = re.compile(r"^src\.[A-Za-z0-9_]+$")
_MD5_LINE = re.compile(r"^([0-9a-fA-F]{32}) [ *](?:\./)?(\S.*)$")
TREE_MD5 = "tree.md5"


def _md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def sourcemods_sets(config_dir=SOURCEMODS_DIR):
    """The set names in config_dir (its subfolders)."""
    d = Path(config_dir)
    return sorted(p.name for p in d.iterdir() if p.is_dir() and not p.name.startswith(".")) if d.is_dir() else []


def _set_dir(where, value, config_dir):
    text = os.path.expanduser(str(value))
    if os.path.isabs(text):
        if not os.path.isdir(text):
            raise ValueError(f"{where}: not a directory")
        return Path(os.path.normpath(text))
    if not _SET_NAME.match(text):
        raise ValueError(f"{where}: give a set name or an absolute path")
    root = Path(config_dir, text)
    if not root.is_dir():
        sets = ", ".join(["none"] + sourcemods_sets(config_dir))
        raise ValueError(f"{where}: no SourceMods set of that name (sets: {sets}); add "
                         f"config/sourcemods/{text}/ or give an absolute path")
    return root


def _tree_md5(where, path):
    """{tree path: md5} from a tree.md5 file."""
    bases = {}
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        m = _MD5_LINE.match(line.rstrip())
        if not m or os.path.isabs(m.group(2)):
            raise ValueError(f"{where}: tree.md5 line {n} is not '<md5>  <path in the CAM tree>'")
        bases[m.group(2)] = m.group(1).lower()
    return bases


def sourcemods(kind, value, config_dir=SOURCEMODS_DIR):
    """The SourceMods set `value` names (none, a set name in config_dir or an absolute path).

    none: empty files/bases/replaces. A set adds dir and md5; bases is {tree path: md5} from tree.md5
    and replaces maps a file to the tree.md5 path of the same name. Raises ValueError."""
    if kind not in ("spinup", "coupled"):
        raise ValueError(f"kind must be 'spinup' or 'coupled', not {kind!r}")
    if value in ("none", None):
        return {"files": {}, "bases": {}, "replaces": {}}
    where = f"{kind}.sourcemods={value!r}"
    root = _set_dir(where, value, config_dir)
    files, md5 = {}, {}
    for top in sorted(root.iterdir()):
        if top.name.startswith(".") or top.name.startswith("README") or top.name == TREE_MD5:
            continue
        if not top.is_dir():
            raise ValueError(f"{where}: {top.name!r} is not inside a src.<component>/ folder (e.g. src.cam/)")
        if not _COMP_DIR.match(top.name):
            raise ValueError(f"{where}: {top.name!r} is not a src.<component>/ folder (e.g. src.cam/)")
        for f in sorted(top.iterdir()):
            if f.name.startswith("."):
                continue
            rel = f"{top.name}/{f.name}"
            # CIME compiles only the files directly in SourceMods/src.<component>/
            if f.is_dir():
                raise ValueError(f"{where}: {rel!r} is a folder; SourceMods are files directly in "
                                 f"src.<component>/")
            if not f.is_file():
                raise ValueError(f"{where}: {rel!r} is not a regular file")
            files[rel], md5[rel] = str(f), _md5(f)
    if not files:
        raise ValueError(f"{where}: {root} holds no src.<component>/ files")
    bases = _tree_md5(where, root / TREE_MD5) if (root / TREE_MD5).is_file() else {}
    by_name = {Path(t).name: t for t in bases}
    replaces = {rel: by_name[Path(rel).name] for rel in files if Path(rel).name in by_name}
    return {"dir": str(root), "files": files, "md5": md5, "bases": bases, "replaces": replaces}


def copy_sourcemods(caseroot, kind, policy, dry_run=False, config_dir=SOURCEMODS_DIR, expect=None):
    """Copy a SourceMods set into <caseroot>/SourceMods/ and return the relative paths copied.

    expect ({rel: md5} from the plan) refuses a set whose files changed since."""
    sm = sourcemods(kind, policy, config_dir)
    where = f"{kind}.sourcemods={policy!r}"
    if expect is not None:
        got = sm.get("md5", {})
        changed = sorted(r for r in set(got) | set(expect) if got.get(r) != expect.get(r))
        if changed:
            raise RuntimeError(f"{where}: {', '.join(changed)} changed since waxrun planned the case; "
                               f"run waxrun again")
    if not dry_run:
        have = sorted(p.name for p in Path(caseroot, "SourceMods").glob("src.*") if p.is_dir())
        missing = sorted({rel.split("/")[0] for rel in sm["files"]} - set(have))
        if missing:
            raise FileNotFoundError(f"{where}: the case has no SourceMods/{missing[0]}/ (it has "
                                    f"{', '.join(have) or 'none'})")
    copied = []
    for rel, src in sorted(sm["files"].items()):
        if not dry_run:
            shutil.copy2(src, Path(caseroot, "SourceMods", rel))
        print(f"[waxrun] {'(dry-run) ' if dry_run else ''}SourceMods/{rel} <- {src}")
        copied.append(rel)
    return copied


_MODSYS = re.compile(r"^([ \t]*)<module_system\b.*?</module_system>", re.S | re.M)


def _module_system_block(block_xml):
    text = Path(block_xml).read_text()
    root = ET.fromstring(text)
    if root.tag != "module_system":
        raise ValueError(f"{block_xml}: root element is <{root.tag}>, expected <module_system>")
    m = _MODSYS.search(text)
    return text[m.start() + len(m.group(1)):m.end()]


def splice_module_system(caseroot, block_xml=MODULE_SYSTEM_XML, path=None):
    """Replace only the <module_system> element of env_mach_specific.xml with block_xml's.

    Writes the .waxrun_spliced guard: case.setup --reset would regenerate the file from the tree."""
    target = Path(path or os.path.join(caseroot, "env_mach_specific.xml"))
    text = target.read_text()
    root = ET.fromstring(text)
    found = root.findall("module_system")
    if len(found) != 1:
        raise ValueError(f"{target}: expected exactly one <module_system>, found {len(found)}")
    block = _module_system_block(block_xml)
    m = _MODSYS.search(text)
    new = text[:m.start() + len(m.group(1))] + block + text[m.end():]
    ET.fromstring(new)
    bak = target.with_name(target.name + ".bak_tree")
    if not bak.exists():
        shutil.copy2(target, bak)
    target.write_text(new)
    # no time stamp on a dry-run copy, so a replay reproduces it byte for byte
    stamp = "" if path else f"{datetime.now().isoformat(timespec='seconds')}\n"
    Path(caseroot, SPLICE_GUARD).write_text(f"{os.path.abspath(block_xml)}\n{stamp}")


def _match(value, xml_value):
    if xml_value.startswith("!"):
        return re.match(xml_value[1:] + "$", str(value)) is None
    return re.match(xml_value + "$", str(value)) is not None


def _group_applies(node, settings):
    for attrib, xml_value in node.attrib.items():
        if attrib in ("name", "source", "queue", "unit_testing"):
            continue
        if not _match(settings.get(attrib.upper(), ""), xml_value):
            return False
    return True


def module_commands(xml_path, settings=None):
    """[(name, argument)] module commands in xml_path whose conditions match settings, as CIME selects them."""
    settings = {**DEFAULT_CASE_SETTINGS, **(settings or {})}
    root = ET.parse(xml_path).getroot()
    ms = root if root.tag == "module_system" else root.find("module_system")
    if ms is None:
        raise ValueError(f"{xml_path}: no <module_system>")
    out = []
    for group in ms.findall("modules"):
        if not _group_applies(group, settings):
            continue
        for cmd in group.findall("command"):
            if _group_applies(cmd, settings):
                out.append((cmd.get("name"), (cmd.text or "").strip()))
    return out


def runtime_modules(xml_path, settings=None):
    """The de-duplicated module load list implied by module_commands.

    Loads before a purge are kept: they are sticky modules that survive it under Lmod."""
    mods = []
    for name, arg in module_commands(xml_path, settings):
        if name == "purge":
            continue
        if name == "load":
            if arg not in mods:
                mods.append(arg)
        elif name in ("unload", "rm"):
            mods = [m for m in mods if m != arg]
        elif name in ("switch", "swap"):
            parts = arg.split()
            mods = [m for m in mods if m != parts[0]] + [parts[-1]]
        else:
            raise ValueError(f"{xml_path}: unsupported module command {name!r}")
    return mods


def module_script(commands):
    """The module commands as one shell line, in order, consecutive loads grouped."""
    out, cur = [], None
    for name, arg in commands:
        if name == "purge":
            out.append("module purge")
            cur = None
        elif cur and cur[0] == name:
            cur[1].append(arg)
        else:
            cur = [name, [arg]]
            out.append(cur)
    return "; ".join(o if isinstance(o, str) else f"module {o[0]} " + " ".join(o[1]) for o in out)


def module_load_line(modules, strict=False):
    """A shell line that purges and loads modules; a str passes through.

    strict loads one module per command, since Lmod exits 0 when one of several loads fails."""
    if isinstance(modules, str):
        return modules
    if strict:
        return " && ".join(["module purge"] + [f"module load {m}" for m in modules])
    return "module purge; module load " + " ".join(modules)


def f7_preflight(caseroot, srcroot, settings=None):
    """Findings before a build: the spliced modules load, and the tree's SFC is ifort."""
    findings = []
    ems = os.path.join(caseroot, "env_mach_specific.xml")
    settings = {**case_settings(caseroot), **(settings or {})}
    mods = runtime_modules(ems, settings)
    line = module_load_line(mods, strict=True)
    r = subprocess.run(["bash", "-lc", line], capture_output=True, text=True)
    if r.returncode == 0:
        findings.append(finding("OK", "module_dry_load", line))
    else:
        findings.append(finding("FAIL", "module_dry_load",
                                f"rc={r.returncode}: {line}\n{r.stderr.strip()}"))
    cmake = os.path.join(srcroot, "ccs_config", "machines", "cmake_macros", "intel.cmake")
    try:
        sfc = re.findall(r'^\s*set\(SFC\s+"?([^")\s]+)"?\)', Path(cmake).read_text(), re.M)
    except OSError as e:
        sfc, cmake = [], f"{cmake} ({e.strerror})"
    if sfc == ["ifort"]:
        findings.append(finding("OK", "sfc_ifort", f"{cmake}: set(SFC \"ifort\")"))
    else:
        findings.append(finding("FAIL", "sfc_ifort", f"{cmake}: set(SFC ...) -> {sfc or 'absent'}, "
                                f"expected ['ifort'] (an ifx flip breaks the classic-compiler build)"))
    for f in findings:
        print(f"[waxrun] preflight {f['level']:4s} {f['check']}: {f['detail'].splitlines()[0]}")
    return findings


_LDD = re.compile(r"^\s*(\S+)\s+=>\s+(.*?)(?:\s+\(0x[0-9a-f]+\))?\s*$")
# release tags of ncarenv and Cray PE install paths
_STACK = re.compile(r"/glade/u/apps/(?:derecho|common)/(\d+\.\d+)/|/opt/cray/pe/([a-z-]+)/(\d[\d.]*)/")


def _release_tags(s):
    return {("ncar", m[0]) if m[0] else (m[1], m[2]) for m in _STACK.findall(s)}


def _ldd(exe, module_list):
    line = module_load_line(module_list) if module_list else "module purge"
    r = subprocess.run(["bash", "-lc", f'{line}; echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"; '
                                       f"ldd {shlex.quote(str(exe))}"],
                       capture_output=True, text=True)
    if r.returncode != 0 and "=>" not in r.stdout:
        raise RuntimeError(f"ldd {exe} failed rc={r.returncode}: {r.stderr.strip()}")
    out, ldpath = {}, ""
    for ln in r.stdout.splitlines():
        if ln.startswith("LD_LIBRARY_PATH="):
            ldpath = ln[len("LD_LIBRARY_PATH="):]
            continue
        m = _LDD.match(ln)
        if m:
            out[m.group(1)] = None if m.group(2).strip() == "not found" else m.group(2).strip()
    return out, ldpath


def ldd_map(exe, module_list):
    """{soname: resolved path or None} from ldd exe after loading module_list in a login shell."""
    return _ldd(exe, module_list)[0]


def ldd_not_found(exe, module_list):
    return sum(1 for p in ldd_map(exe, module_list).values() if p is None)


def ldd_foreign_stack(exe, module_list):
    """Resolved libraries from a software-stack release the loaded modules do not put on LD_LIBRARY_PATH.

    RPATH lets a binary from another stack resolve every library, so ldd_not_found alone misses it."""
    resolved, ldpath = _ldd(exe, module_list)
    allowed = _release_tags(ldpath)
    families = {fam for fam, _ in allowed}
    if not allowed:
        return []
    return sorted(p for p in resolved.values()
                  if p and any(fam in families and (fam, ver) not in allowed
                               for fam, ver in _release_tags(p)))


def _newer_sources(srcroot, exe):
    t = os.path.getmtime(exe)
    hits = []
    for dirpath, dirnames, filenames in os.walk(srcroot):
        dirnames[:] = [d for d in dirnames if d != ".git"]
        for fn in filenames:
            if fn.endswith((".F90", ".f90")):
                p = os.path.join(dirpath, fn)
                try:
                    if os.path.getmtime(p) > t:
                        hits.append(p)
                except OSError:
                    pass
    return hits


def build_plan(spec):
    """Shell lines of the build / verify_only policy for commands.sh."""
    b = spec.get("build") or {}
    mode = b.get("mode", "build")
    caseroot = spec["caseroot"]
    exe = spec.get("cesm_exe") or os.path.join(spec["exeroot"], "cesm.exe")
    mods = spec.get("runtime_modules") or []
    lines = [f"cd {shlex.quote(caseroot)}"]
    if mode == "build":
        if b.get("clean_all", True):
            lines.append("./case.build --clean-all")
        lines.append(f"./case.build 2>&1 | tee {BUILD_LOG}")
        exe = os.path.join(spec["exeroot"], "cesm.exe")
    if mode in ("build", "verify_only"):
        lines.append(f"# assert case_solver.ldd_not_found({exe}, runtime_modules) == 0")
        if spec.get("voltron_exe"):
            lines.append(f"# assert case_solver.ldd_not_found({spec['voltron_exe']}, runtime_modules) == 0")
        if mods:
            lines.append(f"bash -lc {shlex.quote(module_load_line(mods) + '; ldd ' + exe + ' | grep -c \"not found\"')}")
    if mode == "verify_only" and b.get("no_newer_sources"):
        same_mods = ", SourceMods equal" if (spec.get("sourcemods") or {}).get("files") else ""
        lines.append(f"# build.share_exe: {exe} is the {b.get('shared_from', 'other')} case's build: same SRCROOT, "
                     f"no *.F90 under SRCROOT newer than it{same_mods}")
    return lines


def build(spec, dry_run=False, build_window=3600):
    """Run the build or verify_only policy and check the result; RuntimeError on a failed check.

    dry_run returns build_plan(spec) without running anything."""
    plan = build_plan(spec)
    if dry_run:
        return plan
    b = spec.get("build") or {}
    mode = b.get("mode", "build")
    caseroot = spec["caseroot"]
    mods = spec.get("runtime_modules") or []
    if mode not in ("build", "verify_only"):
        raise ValueError(f"unknown build mode {mode!r}")
    if mode == "build":
        if b.get("clean_all", True):
            _run(["./case.build", "--clean-all"], cwd=caseroot)
        _run(["bash", "-c", f"./case.build 2>&1 | tee {BUILD_LOG}; exit ${{PIPESTATUS[0]}}"],
             cwd=caseroot)
        exe = os.path.join(spec["exeroot"], "cesm.exe")
    else:
        exe = spec["cesm_exe"]
    for e in [exe] + ([spec["voltron_exe"]] if spec.get("voltron_exe") else []):
        if e != exe and not os.path.exists(e):
            # voltron may not be built yet; the coupled job's preflight re-checks it
            print(f"[waxrun] WARNING: {e} does not exist yet; skipping its library check. "
                  f"Build voltron_mpi.x before submitting the coupled deck.", file=sys.stderr)
            continue
        n = ldd_not_found(e, mods)
        if n:
            raise RuntimeError(f"{e}: {n} libraries not found under {module_load_line(mods)}")
    if mode == "build":
        atm_in = os.path.join(caseroot, "CaseDocs", "atm_in")
        if os.path.exists(atm_in) and os.path.getmtime(exe) <= os.path.getmtime(atm_in) - build_window:
            raise RuntimeError(f"{exe} is older than {atm_in} by more than {build_window} s: "
                               f"namelists and binary are not from one build")
    else:
        _assert_same_srcroot(b.get("source_caseroot"), caseroot)
        if b.get("no_newer_sources"):
            srcroot = xmlquery(caseroot, "SRCROOT")["SRCROOT"]
            newer = _newer_sources(srcroot, exe)
            if newer:
                raise RuntimeError(f"{len(newer)} source files under {srcroot} are newer than {exe}, "
                                   f"e.g. {newer[0]}")
        if b.get("shared_from"):
            _assert_shared_sourcemods(b.get("source_caseroot"), (spec.get("sourcemods") or {}).get("md5") or {})
    return plan


def _assert_shared_sourcemods(source_caseroot, md5):
    """The case that built the shared exe holds the same SourceMods files ({rel: md5})."""
    for rel, want in sorted(md5.items()):
        f = Path(source_caseroot, "SourceMods", rel)
        if not f.is_file() or _md5(f) != want:
            raise RuntimeError(f"{f} is not the coupled.sourcemods file {rel}: rebuild the spinup case")


def _canonical_path(p):
    return os.path.realpath(os.path.normpath(str(p)))


def _assert_same_srcroot(source_caseroot, caseroot):
    if not source_caseroot:
        raise RuntimeError("build.source_caseroot is required to verify SRCROOT equality")
    # a stored SRCROOT may be spelled '<tree>/cime/..'
    a = xmlquery(source_caseroot, "SRCROOT")["SRCROOT"]
    c = xmlquery(caseroot, "SRCROOT")["SRCROOT"]
    if not a or not c or _canonical_path(a) != _canonical_path(c):
        raise RuntimeError(f"SRCROOT mismatch: {source_caseroot} -> {a}, {caseroot} -> {c}")
