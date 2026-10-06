"""Identify the voltron binary a coupled run launches, and find its restart files."""

import datetime
import glob
import hashlib
import os
import re

# gam files carry the MPI tile in the name; a missing dk.Res is survivable (DragonKing cold-starts).
RES_FAMILIES = ("gam", "volt", "dk", "mix", "raiju", "raiCpl", "gamCpl")
RES_OPTIONAL = ("dk",)

# Printed by voltapp_mpi.F90 when WACCM-X joins the MPMD launch.
WCMX_MARKER = "ranks of W"

_CACHE = {}


def is_elf(path):
    try:
        with open(path, "rb") as f:
            return f.read(4) == b"\x7fELF"
    except OSError:
        return False


def file_sha256(path):
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


def printable_strings(data, n=4):
    """What `strings -a -n n` prints for `data` (bytes)."""
    return [m.decode("ascii") for m in re.findall(rb"[\t\x20-\x7e]{%d,}" % n, data)]


def _file_facts(path):
    try:
        st = os.stat(path)
    except OSError:
        return None
    return {"size": st.st_size, "mtime": datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")}


def _non_elf_hint(path):
    try:
        with open(path, "rb") as f:
            head = f.read(80)
        size = os.path.getsize(path)
    except OSError:
        return "unreadable"
    first = head.split(b"\n", 1)[0].decode("ascii", "replace")
    return f"{size} bytes, first line {first!r}"


def voltron_identity(exe):
    """Facts about voltron_mpi.x `exe` (elf, sha256, wcmx, ...), cached per (realpath, size, mtime)."""
    out = {"path": exe, "realpath": os.path.realpath(exe) if exe else None,
           "exists": bool(exe) and os.path.isfile(exe), "elf": False, "sha256": None, "size": None,
           "mtime": None, "wcmx": False}
    if not out["exists"]:
        return out
    facts = _file_facts(exe) or {}
    key = ("voltron", out["realpath"], facts.get("size"), facts.get("mtime"))
    if key in _CACHE:
        return {**_CACHE[key], "path": exe}
    out.update(facts)
    out["elf"] = is_elf(exe)
    try:
        with open(exe, "rb") as f:
            data = f.read()
    except OSError:
        data = b""
    out["sha256"] = hashlib.sha256(data).hexdigest() if data else file_sha256(exe)
    if out["elf"]:
        out["wcmx"] = WCMX_MARKER in "\n".join(printable_strings(data))
    _CACHE[key] = dict(out)
    return out


def voltron_problem(ident, component):
    """(level, message) for a voltron `ident`, level 'missing' or 'refuse'; (None, None) when fine."""
    path = ident.get("path")
    if not ident.get("exists"):
        return "missing", f"{path} does not exist"
    if not ident.get("elf"):
        return "refuse", (f"{path} is not an ELF executable ({_non_elf_hint(path)}): it cannot be "
                          f"voltron_mpi.x. Point pbs.kaiju_build_directory at a kaiju build whose "
                          f"bin/voltron_mpi.x is the compiled voltron.")
    if component == "grx" and not ident.get("wcmx"):
        return "refuse", (f"{path} has no 'ranks of W': it is not a voltron that couples WACCM-X, so "
                          f"CAM would hang in its first import broadcast. Build voltron from a uidev "
                          f"checkout.")
    return None, None


def res_files(prefix, index, family):
    """Existing restart file(s) of `family` at `index` for run prefix <dir>/<runid>; dangling links skipped."""
    if family == "gam":
        return sorted(p for p in glob.glob(f"{glob.escape(prefix)}*.gam.Res.{int(index):05d}.h5")
                      if os.path.exists(p))
    p = f"{prefix}.{family}.Res.{int(index):05d}.h5"
    return [p] if os.path.exists(p) else []


def job_probe_lines(exe="./voltron_mpi.x", sha256=None):
    """Bash lines that warn at launch when `exe` is not the binary the deck was written for."""
    return ["# ---- voltron check: the binary this deck was written for vs the one this job runs",
            f"want_sha={sha256 or ''}",
            f"have_sha=$(sha256sum {exe} 2>/dev/null | cut -d' ' -f1)",
            'if [ "$have_sha" != "$want_sha" ]; then',
            f'  echo "WARNING: {exe} (sha256 ${{have_sha:-unreadable}}) is not the binary engage recorded '
            f'for this deck (sha256 ${{want_sha:-none, it did not exist yet}})"',
            "fi",
            "# ---- end voltron check"]
