#!/usr/bin/env python


"""engage: configure and generate a coupled MAGE run (GAMERA + RAIJU + TIE-GCM or WACCM-X).

Interactive by default; --mode BASIC|INTERMEDIATE|EXPERT sets how many options are prompted.
"""

import argparse
import contextlib
import copy
import datetime
import glob
import importlib
import io
import json
import os
import re
import shlex
import sys
import math

# No environment is required at import time (so -h works); main() checks what the component needs.
_HERE = os.path.dirname(os.path.abspath(__file__))
KAIJU_ENV_HINT = "source kaiju-private/scripts/setupEnvironment.sh"
TIEGCM_ENV_HINT = "source tiegcm/tiegcmrun/setEnvironment.sh (sets TIEGCMHOME and TIEGCMDATA)"

KAIJUHOME = os.environ.get("KAIJUHOME") or os.path.normpath(os.path.join(_HERE, "..", ".."))

# tiegcmrun itself is imported lazily (grt / tiegcm only); its dir goes on sys.path before makeitso's.
if os.environ.get("TIEGCMHOME"):
    sys.path.append(f'{os.environ["TIEGCMHOME"]}/tiegcmrun')

sys.path.append(f'{KAIJUHOME}/scripts/makeitso')
import makeitso

if _HERE not in sys.path:
    sys.path.append(_HERE)
import grx_options
import binident


def _import_own_jobgen():
    """Return JobGen from scripts/engage/jobgen.py, whichever tool's jobgen is already loaded.

    The tools all `from jobgen import JobGen`; if another copy holds that name, load ours as engage_jobgen."""
    import importlib.util
    path = os.path.join(_HERE, "jobgen.py")
    loaded = sys.modules.get("jobgen")
    if loaded is not None and os.path.abspath(getattr(loaded, "__file__", "") or "") == path:
        return loaded.JobGen
    name = "engage_jobgen" if loaded is not None else "jobgen"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    return sys.modules[name].JobGen


def _import_tiegcmrun():
    """Import tiegcmrun from $TIEGCMHOME/tiegcmrun on first use, after checking the TIE-GCM env."""
    global tiegcmrun
    mod = globals().get("tiegcmrun")
    if mod is not None:
        return mod
    tiegcmhome = makeitso.require_env("TIEGCMHOME", TIEGCM_ENV_HINT)
    makeitso.require_env("TIEGCMDATA", TIEGCM_ENV_HINT)
    tg_dir = os.path.join(tiegcmhome, "tiegcmrun")
    if not os.path.isfile(os.path.join(tg_dir, "tiegcmrun.py")):
        raise SystemExit(f"ERROR: TIEGCMHOME={tiegcmhome} holds no tiegcmrun/tiegcmrun.py; "
                         f"{TIEGCM_ENV_HINT}.")
    if tg_dir not in sys.path:
        sys.path.append(tg_dir)
    tiegcmrun = importlib.import_module("tiegcmrun")
    return tiegcmrun


def __getattr__(name):
    if name == "tiegcmrun":
        return _import_tiegcmrun()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


JobGen = _import_own_jobgen()

DESCRIPTION = "Configure and generate a MAGE run (GRT, GRX or GR) or a standalone TIE-GCM run."

# API_VERSION major engage requires of each tool it drives from another checkout.
REQUIRED_API = {"tiegcmrun": 2, "waxrun": 1}


def check_api_version(name, module):
    """Exit if the tool's API_VERSION major is not REQUIRED_API[name]."""
    have = getattr(module, "API_VERSION", None)
    want = REQUIRED_API[name]
    try:
        major = int(str(have).split(".")[0])
    except ValueError:
        major = None
    if major != want:
        raise SystemExit(
            f"ERROR: engage drives {name} API version {want}.x, but "
            f"{getattr(module, '__file__', name)} reports API_VERSION={have!r}: update the {name} "
            f"checkout (or engage) so that they match.")
    return module

JSON_INDENT = 4

SUPPORT_FILES_DIRECTORY = os.path.join(KAIJUHOME, "scripts", "engage")

MAKEITSO_FILES_DIRECTORY = os.path.join(KAIJUHOME, "scripts", "makeitso")

OPTION_ENGAGE_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_engage_descriptions.json"
)

OPTION_MAKEITSO_DESCRIPTIONS_FILE = os.path.join(
    MAKEITSO_FILES_DIRECTORY, "option_descriptions.json"
)

CONFIG_DIR = os.path.join(SUPPORT_FILES_DIRECTORY, "config")

def layout_help(lead=makeitso.GR_SPINUP_DIR):
    """engage --help epilog: where the run lands and where the layout is documented."""
    return (f"The run is written to pbs.run_directory; a segmented run's GR lead-in runs in {lead}/. "
            f"Layout and replay: docs/source/makeitso/engage.rst and the GR, GRT and GRX pages.")


LAYOUT_HELP = layout_help()


def create_command_line_parser():
    """Return the engage argument parser."""
    parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=LAYOUT_HELP)
    parser.add_argument(
        "--clobber", action="store_true",
        help="Overwrite the run's existing replay files (default: %(default)s)."
    )
    parser.add_argument(
        "--recreate-cases", action="store_true", default=False,
        help="GRX: delete and recreate existing WACCM-X cases (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--mode", default="BASIC", type=str.upper, choices=makeitso.MODES, metavar="MODE",
        help="User mode: BASIC|INTERMEDIATE|EXPERT, case-insensitive (default: %(default)s)."
    )
    parser.add_argument(
        "--component", "-m", default="grt", choices=["grt", "grx", "tiegcm", "gr"],
        help="Run type: coupled grt or grx, standalone tiegcm or gr (default: %(default)s)."
    )
    parser.add_argument(
        "--engage_options_path", "-eo", default=None,
        help="Path to engage JSON file of options (default: %(default)s)."
    )
    parser.add_argument(
        "--makeitso_options_path", "-mo", default=None,
        help="Path to makeitso JSON file of options (default: %(default)s)."
    )
    parser.add_argument(
        "--tiegcm_options_path", "-to", default=None,
        help="Path to tiegcm JSON file of options (default: %(default)s)."
    )
    parser.add_argument(
        "--rederive", default=None, metavar="GROUP[,GROUP]",
        help="With -eo/-mo/-to: re-derive these groups; '--rederive list' names them."
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "--waxrun_options_path", "-wo", default=None,
        help="GRX: path to waxrun JSON file of options (default: %(default)s)."
    )
    parser.add_argument(
        "--dry-run", action="store_true", default=False,
        help="GRX: write the case specs as plain files; run no CIME script."
    )
    parser.add_argument(
        "--skip-build", action="store_true", default=False,
        help="GRX: skip case.build; every cesm.exe must already exist."
    )
    parser.add_argument(
        "--phase", default="both", choices=["spinup", "coupled", "both"],
        help="GRX: which WACCM-X case(s) to generate (default: %(default)s)."
    )
    parser.add_argument(
        "--no-tui", action="store_true",
        help="Ask the options as linear prompts instead of the wizards."
    )
    parser.add_argument(
        "--tui", action="store_true",
        help="The wizards (the default; kept for old command lines)."
    )
    return parser


def _tiegcm_misc():
    """Return tiegcmrun's misc module."""
    _import_tiegcmrun()
    return importlib.import_module("misc")


def build_grt_pbs_options(gr_options: dict, makeitso_options: dict, tiegcm_options: dict) -> dict:
    """Return gr_options with the coupled-GRT task-group layouts and MPMD launch fields added to pbs."""
    options = copy.deepcopy(gr_options)

    # Both models run in one job, so a differing tiegcm hpc_system would give a mixed-machine deck.
    hpc_system = gr_options["simulation"]["hpc_system"]
    tiegcm_hpc = tiegcm_options["simulation"]["hpc_system"]
    if tiegcm_hpc != hpc_system:
        raise ValueError(
            f"Coupled GRT run: tiegcm simulation.hpc_system={tiegcm_hpc!r} differs from engage's "
            f"simulation.hpc_system={hpc_system!r}; both models run in one {hpc_system} job. "
            f"Use the same machine for both (engage's value is authoritative).")

    options["pbs"]["tie_nodes"] = tiegcm_options["job"]["resource"]["select"]
    options["pbs"]["tie_mpiprocs"] = tiegcm_options["job"]["resource"]["mpiprocs"]
    options["pbs"]["tie_mpiranks"] = tiegcm_options["job"]["nprocs"]
    options["pbs"]["tie_exe"] = tiegcm_options["model"]["data"]["coupled_modelexe"]

    options["pbs"]["gamera_nodes"] = makeitso_options["pbs"]["select"]
    options["pbs"]["gamera_mpiprocs"] = makeitso_options["pbs"]["mpiprocs"]
    options["pbs"]["gamera_ompthreads"] = makeitso_options["pbs"]["ompthreads"]
    options["pbs"]["voltron_nodes"] = makeitso_options["pbs"]["select2"]
    options["pbs"]["voltron_mpiprocs"] = makeitso_options["pbs"]["helper_mpiprocs"]
    options["pbs"]["voltron_ompthreads"] = makeitso_options["pbs"]["helper_ompthreads"]
    options["pbs"]["voltron_mpiranks"] = int(options["pbs"]["gamera_nodes"])*int(options["pbs"]["gamera_mpiprocs"])+int( options["pbs"]["voltron_nodes"])*int(options["pbs"]["voltron_mpiprocs"])
    # placer wrapper of makeitso's launcher, e.g. "mpiexec pinCpuCores.sh"
    options["pbs"]["voltron_scripts"] = makeitso.launcher(hpc_system).replace("mpiexec ", "")
    nodes = (int(options["pbs"]["tie_nodes"]) + int(options["pbs"]["gamera_nodes"])
             + int(options["pbs"]["voltron_nodes"]))
    err = grx_options.ncpus_error(hpc_system, options["pbs"].get("queue"), nodes,
                                  "the coupled GRT job (TIE-GCM + GAMERA + voltron nodes)")
    if err:
        raise ValueError(err)

    if hpc_system not in ("derecho", "aitken", "linux", "perlmutter"):
        raise NotImplementedError(
            f"Coupled GRT runs support derecho, aitken, linux and (experimental) perlmutter. "
            f"hpc_system={hpc_system!r} has no coupled launch dispatch in "
            f"build_grt_pbs_options (mpiexec command/option, node-file split). Add a "
            f"branch for it, or run the standalone tools.")
    if hpc_system == "derecho":
        options["pbs"]["mpiexec_command"] = "mpiexec"
        options["pbs"]["mpiexec_option"] = "-n"
    elif hpc_system == "linux":
        options["pbs"]["mpiexec_command"] = "mpiexec"
        options["pbs"]["mpiexec_option"] = "-n"
        options["pbs"]["voltron_scripts"] = ""
    elif hpc_system == "perlmutter":
        # SLURM het-job components: tiegcm = 0, gamera = 1, voltron helpers = 2.
        options["pbs"]["mpiexec_command"] = "srun"
        options["pbs"]["mpiexec_option"] = "-n"
        options["pbs"]["tie_scripts"] = None
        options["pbs"]["voltron_scripts"] = ""
        options["pbs"]["het_groups"] = ["0", "1,2"]
    elif hpc_system == "aitken":
        options["pbs"]["mpiexec_command"] = "mpiexec_mpt"
        options["pbs"]["mpiexec_option"] = "-np"
        options["pbs"]["tie_scripts"] = "correctOMPenvironment.sh $NODEFILE_1 omplace"
        options["pbs"]["voltron_scripts"] = "correctOMPenvironment.sh $NODEFILE_2 omplace"
        options["pbs"]["nodecommand"] = f"""
# This section creates two node files based on the master node file.
# Each node file should include the nodes that contribute to an executable.
# Currently, head corresponds to TIEGCM ranks and tail are the GR ranks
# A more robust system is needed in the future. GR ranks is G + helper R
export NODEFILE_1=${{NODEFILE}}.1
export NODEFILE_2=${{NODEFILE}}.2
head -n {options["pbs"]["tie_mpiranks"]} $NODEFILE > $NODEFILE_1
tail -n {options["pbs"]["voltron_mpiranks"]} $NODEFILE > $NODEFILE_2
wc -l $NODEFILE

export NODEFILE_T=${{NODEFILE}}.T
cat $NODEFILE_1 $NODEFILE_2 > $NODEFILE_T
echo ""
echo "Original Nodes"
cat $PBS_NODEFILE
echo ""

echo ""
echo "New Nodes"
cat $NODEFILE_T
echo ""

export PBS_NODEFILE=$NODEFILE_T

echo ""
echo "PBS_NODEFILE = "
echo $PBS_NODEFILE
echo ""

echo "Running tiegcm and voltron at the same time"
"""
    return options


def render_grt_pbs(options: dict) -> str:
    """Render a coupled GRT job script from the machine/batch config and build_grt_pbs_options()."""
    sim = options["simulation"]
    pbs = options["pbs"]
    jg = JobGen(CONFIG_DIR, machine=sim["hpc_system"], node_type=sim.get("node_type"))

    tie = {"mpi_tasks": int(pbs["tie_nodes"]) * int(pbs["tie_mpiprocs"]),
           "tasks_per_node": int(pbs["tie_mpiprocs"]), "threads": 1}
    gamera = {"mpi_tasks": int(pbs["gamera_nodes"]) * int(pbs["gamera_mpiprocs"]),
              "tasks_per_node": int(pbs["gamera_mpiprocs"]), "threads": int(pbs["gamera_ompthreads"])}
    voltron = {"mpi_tasks": int(pbs["voltron_nodes"]) * int(pbs["voltron_mpiprocs"]),
               "tasks_per_node": int(pbs["voltron_mpiprocs"]), "threads": int(pbs["voltron_ompthreads"])}
    task_groups = [tie, gamera, voltron]

    segment_id = sim["segment_id"]
    jobid = jg.sched("jobid")
    jobname = jg.sched("jobname", job_name=segment_id)
    build_dir = pbs["kaiju_build_directory"]
    # cd to the absolute run directory, not $PBS_O_WORKDIR, so the job can be queued from anywhere.
    run_directory = pbs.get("run_directory")
    workdir_cd = (f"cd {shlex.quote(os.path.abspath(run_directory))}" if run_directory
                  else jg.sched("workdir_cd"))
    job_dir = os.path.abspath(run_directory) if run_directory else os.getcwd()

    def job_spelling(path):
        """Spell a job path relative to the job's directory; paths outside it stay absolute."""
        if path is None or os.path.isabs(str(path)) or "$" in str(path):
            return JobGen.spell(path, job_dir, job_dir)
        return JobGen.spell(os.path.join(job_dir, str(path)), job_dir, job_dir)

    L = [jg.shebang(), "",
         f"# This script was generated to run on {sim['hpc_system']}.", ""]
    pre = jg.directive_lines({
        "job_name": segment_id,
        "account": pbs.get("account_name"),
        "group_list": pbs.get("group_list"),
        "queue": pbs.get("queue"),
    })
    res = jg.resource_lines(task_groups, pbs["walltime"], queue=pbs.get("queue")) + makeitso.constraint_lines(jg)
    post = jg.directive_lines({
        "priority": pbs.get("job_priority"),
        "mail": "abe",
        "join_out": True,
    })
    # Job-wide directives first: on SLURM they would otherwise bind to the last het-job component.
    L += pre + post + res

    L += ["",
          f'echo "Job {jobid} started at `date` on `hostname`."', "",
          "# Specify the ID string for the run.",
          f"export RUNID={jobname}", "",
          "# Run in the job directory.",
          workdir_cd, "",
          "# Load the required modules for kaiju."]
    L += jg.module_lines(stack="grt")
    L += ["module list"]
    L += [m for m in makeitso.as_list(pbs.get("local_modules")) if m not in (None, "None")]

    L += ["",
          "# Define the kaiju installation location.",
          f"export KAIJU_INSTALL_DIR={pbs['kaiju_install_directory']}", "",
          "# Set KAIJUHOME and the other kaiju variables.",
          "source $KAIJU_INSTALL_DIR/scripts/setupEnvironment.sh", "",
          "# Add the kaiju build bin directory to PATH.",
          f"export PATH={build_dir}/bin:$PATH", "",
          "# Set other environment variables.",
          "export MPI_TYPE_DEPTH=32",
          "export OMP_STACKSIZE=100M",
          f"export OMP_NUM_THREADS={pbs.get('voltron_ompthreads', pbs.get('helper_ompthreads'))}", "",
          'echo "The active environment variables are:"',
          "printenv"]

    nodefile_src = jg.sched("nodefile")
    if nodefile_src:
        L += ["",
              "# Copy the node file for use by the placer command.",
              f"NODEFILE=nodefile.{jobid}",
              f"cp {nodefile_src} $NODEFILE"]
        if pbs.get("nodecommand"):
            L.append(pbs["nodecommand"])

    het = pbs.get("het_groups") or [None, None]
    tie_het = f" --het-group={het[0]}" if het[0] else ""
    vol_het = f"--het-group={het[1]} " if het[1] else ""
    tie_cmd = f"{pbs['mpiexec_command']}{tie_het} {pbs['mpiexec_option']} {pbs['tie_mpiranks']}"
    if pbs.get("tie_scripts"):
        tie_cmd += f" {pbs['tie_scripts']}"
    tie_cmd += f" $TIEGCM_EXE {job_spelling(pbs['tie_inp'])}"
    vol_cmd = f"{vol_het}{pbs['mpiexec_option']} {pbs['voltron_mpiranks']}"
    if pbs.get("voltron_scripts"):
        vol_cmd += f" {pbs['voltron_scripts']}"
    vol_cmd += " $VOLTRON_EXE $RUNID.xml"
    run_line = f"{tie_cmd} : {vol_cmd} &>> $OUTPUT_FILE"

    # ./voltron_mpi.x is linked in by link_voltron.
    L += ["",
          "# Run the model. Direct output from the program is saved in a text file.",
          "VOLTRON_EXE=./voltron_mpi.x",
          f"TIEGCM_EXE={shlex.quote(job_spelling(pbs['tie_exe']))}",
          f'OUTPUT_FILE="$VOLTRON_EXE-{segment_id}.out"']
    if pbs.get("handoff"):
        L += list(pbs["handoff"])
    if "voltron_sha256" in pbs:
        L += binident.job_probe_lines("./voltron_mpi.x", pbs["voltron_sha256"])
    last_prim = job_spelling(pbs.get("tie_last_prim"))
    L += ['echo "Running $VOLTRON_EXE on model $RUNID."',
          GRT_OUTPUT_MARK]
    if last_prim:
        L.append(GRT_START_STAMP)
    L += [run_line,
          "rc=$?"]
    L += grt_post_run_lines(last_prim)
    L += ["",
          f'echo "Job {jobid} ended at `date` on `hostname`."',
          "exit $rc"]
    return "\n".join(L) + "\n"


# $OUTPUT_FILE is appended to by every launch; the post-run checks read only this launch's lines.
GRT_OUTPUT_MARK = 'out0=$(cat "$OUTPUT_FILE" 2>/dev/null | wc -l)'
# Without this line voltron runs uncoupled and reports no error.
GRT_COUPLING_REGEX = "V coupling to .* ranks of T"
GRT_START_STAMP = 'touch "$OUTPUT_FILE.start"'


def tiegcm_last_primary(inp_path, run_dir=None):
    """Absolute path of the last primary history the namelist's OUTPUT line writes, or None."""
    if run_dir is not None and not os.path.isabs(inp_path):
        inp_path = os.path.join(run_dir, inp_path)
    try:
        with open(inp_path, encoding="utf-8") as f:
            txt = f.read()
        m = re.search(r"^\s*OUTPUT\s*=\s*(.*?)\s*$", txt, re.M)
        if not m:
            return None
        _import_tiegcmrun()
        import output_solver
        last = output_solver.last_primary_history({"OUTPUT": m.group(1)})
    except (OSError, ImportError, AttributeError):
        return None
    if not last:
        return None
    return os.path.normpath(os.path.join(run_dir, last)) if run_dir is not None else os.path.abspath(last)


def grt_post_run_lines(last_prim=None):
    """Job lines after `rc=$?` that fail the segment if TIE-GCM did not couple or stopped early.

    TIE-GCM's MODEL SHUTDOWN can exit 0 (uninitialised errorcode to mpi_abort), so the output is checked too."""
    this_launch = 'tail -n +$((out0 + 1)) "$OUTPUT_FILE"'
    hist = []
    if last_prim:
        q = shlex.quote(last_prim)
        hist = [
            f"if [ ! -e {q} ]; then",
            f"  echo \"ERROR: TIE-GCM did not write {last_prim}, the last primary history of this segment\"",
            "  [ $rc -ne 0 ] || rc=1",
            f'elif [ -z "$(find {q} -newer "$OUTPUT_FILE.start")" ]; then',
            f"  echo \"ERROR: {last_prim} is older than this job: TIE-GCM stopped before the end of the segment\"",
            "  [ $rc -ne 0 ] || rc=1",
            "fi",
        ]
    return [
        "# ---- post-run checks: fail when TIE-GCM did not couple or stopped early",
        f"if ! {this_launch} | grep -Eq '{GRT_COUPLING_REGEX}'; then",
        "  echo \"ERROR: no 'V coupling to <N> ranks of T' line in $OUTPUT_FILE for this launch: TIE-GCM did "
        "not join the coupling (voltron then runs uncoupled without an error)\"",
        "  [ $rc -ne 0 ] || rc=1",
        "fi",
        f"if {this_launch} | grep -q 'MODEL SHUTDOWN'; then",
        "  echo \"ERROR: TIE-GCM stopped with MODEL SHUTDOWN (see $OUTPUT_FILE); its exit status is not "
        "reliable\"",
        "  [ $rc -ne 0 ] || rc=1",
        "fi",
    ] + hist + [
        'if [ $rc -ne 0 ]; then echo "segment failed rc=$rc"; fi',
        "# ---- end post-run checks",
    ]


def deck_count_error(makeitso_options, tiegcm_inp_scripts):
    """Error message if makeitso's pbs.num_segments differs from the TIE-GCM segment count, else None."""
    n_mage = makeitso_options["pbs"].get("num_segments")
    if str(n_mage) != str(len(tiegcm_inp_scripts)):
        return (f"ERROR: the coupled run has {len(tiegcm_inp_scripts)} TIE-GCM segments but "
                f"{n_mage} MAGE coupled decks (pbs.num_segments).")
    return None


def remove_gr_only_scripts(run_dir, job_name):
    """Remove makeitso's GR-only <job>-NN.pbs and <job>_pbs.sh; return the names removed."""
    paths = sorted(glob.glob(os.path.join(run_dir, f"{job_name}-[0-9][0-9].pbs")))
    paths.append(os.path.join(run_dir, f"{job_name}_pbs.sh"))
    removed = []
    for p in paths:
        if os.path.isfile(p):
            os.remove(p)
            removed.append(os.path.basename(p))
    return removed


def create_pbs_scripts(gr_options: dict, makeitso_options: dict, makeitso_pbs_scripts: list,
                       tiegcm_options: dict, tiegcm_inp_scripts: list, tiegcm_pbs_scripts: list,
                       run_dir: str = None, voltron: dict = None):
    """Write one coupled GRT job script per segment and the <job>_pbs.sh submit script.

    Returns (coupled job script paths, submit script path).
    """
    err = deck_count_error(makeitso_options, tiegcm_inp_scripts)
    if err:
        raise SystemExit(err)
    options = build_grt_pbs_options(gr_options, makeitso_options, tiegcm_options)
    if voltron is None:
        voltron = binident.voltron_identity(
            os.path.join(makeitso_options["pbs"]["kaiju_build_directory"], "bin", "voltron_mpi.x"))
    rd = os.path.abspath(run_dir or os.getcwd())
    handoff = makeitso.handoff_plan(makeitso_options, {"coupling": gr_options.get("coupling") or {}})
    runid = makeitso.restart_runid(makeitso_options)

    pbs_scripts = []
    n_coupled = len(tiegcm_inp_scripts)
    print(f'Creating {n_coupled} PBS scripts')
    for job in range(1, n_coupled + 1):

        opt = copy.deepcopy(options)
        runid = opt["simulation"]["job_name"]
        segment_id = f"{runid}-{job:02d}"
        opt["simulation"]["segment_id"] = segment_id

        opt["pbs"]["tie_inp"] = os.path.join(rd, tiegcm_inp_scripts[job-1])
        opt["pbs"]["tie_last_prim"] = tiegcm_last_primary(tiegcm_inp_scripts[job-1], rd)
        opt["pbs"]["voltron_sha256"] = voltron.get("sha256")
        if segment_id in handoff:
            opt["pbs"]["handoff"] = makeitso.handoff_lines(runid, handoff[segment_id])

        pbs_script = makeitso.run_directory_path(opt, f"{opt['simulation']['segment_id']}.pbs")
        # absolute: a relative value would resolve one level too deep from the run directory
        opt["pbs"]["run_directory"] = rd
        pbs_content = render_grt_pbs(opt)
        pbs_scripts.append(pbs_script)
        with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)

    submit_all_jobs_script = f"{gr_options['simulation']['job_name']}_pbs.sh"
    jg = JobGen(CONFIG_DIR, machine=gr_options["simulation"]["hpc_system"],
                node_type=gr_options["simulation"].get("node_type"))
    # Empty lead-in chains are skipped; segment 01 waits only on the chains submitted.
    lead_ins = [("GR lead-in", "makeitso_job_id", makeitso.submit_entries(makeitso_pbs_scripts, rd)),
                ("TIE-GCM spinup", "tiegcm_job_id", makeitso.submit_entries(tiegcm_pbs_scripts, rd))]
    grt_jobs = makeitso.submit_entries(pbs_scripts, rd)
    chains = [(c, v, s, []) for c, v, s in lead_ins if s]
    chains.append(("coupled segments", "job_id", grt_jobs, [f"${v}" for _, v, s in lead_ins if s]))
    with open(submit_all_jobs_script, "w", encoding="utf-8") as f:
        f.write(jg.submit_chain_script(
            chains, header=f"Submit the {gr_options['simulation']['job_name']} GRT run. Run from anywhere.",
            jobids=f"{gr_options['simulation']['job_name']}.jobids"))
    os.chmod(submit_all_jobs_script, 0o755)

    return pbs_scripts, submit_all_jobs_script


# GRX (GAMERA + RAIJU + WACCM-X) launcher: WACCM-X replaces TIE-GCM as the third MPMD task group.

GRX_SEGMENT_DURATION = grx_options.GRX_SEGMENT_DURATION
GRX_GCM_SPIN_UP_TIME = grx_options.GRX_GCM_SPIN_UP_TIME
GRX_DATE_FORMAT = "%Y-%m-%dT%H:%M:%S"
GRX_ARTIFACTS = ("engage_parameters.json", "makeitso_parameters.json", "waxrun_parameters.json")


def _import_waxrun_module(name):
    """Import one of waxrun's modules from scripts/waxrun."""
    wax_dir = os.path.join(KAIJUHOME, "scripts", "waxrun")
    if wax_dir not in sys.path:
        sys.path.append(wax_dir)
    return importlib.import_module(name)


def _import_segment_solver():
    return _import_waxrun_module("segment_solver")


def _import_waxrun():
    return _import_waxrun_module("waxrun")


def _grx_datetime(v):
    if isinstance(v, datetime.datetime):
        return v
    return datetime.datetime.fromisoformat(str(v))


def apply_grx_defaults(engage_options):
    """Fill blank GRX segment/spin-up defaults in a replayed or partial options dict."""
    sim = engage_options["simulation"]
    cpl = engage_options.setdefault("coupling", {})
    # a 0 is kept (and refused by the range checks)
    if sim.get("segment_duration") is None or str(sim["segment_duration"]).strip() == "":
        sim["segment_duration"] = GRX_SEGMENT_DURATION
    if cpl.get("gcm_spin_up_time") is None or str(cpl["gcm_spin_up_time"]).strip() == "":
        cpl["gcm_spin_up_time"] = GRX_GCM_SPIN_UP_TIME
    return engage_options


def _module_name_version(mod):
    name, _, version = str(mod).partition("/")
    return name, version


def _waxrun_runtime_modules(waxrun, hpc_system):
    """waxrun's stacks.<hpc>.runtime_modules default, or None."""
    path = os.path.join(os.path.dirname(os.path.abspath(waxrun.__file__)), "config", "options_description.json")
    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        desc = json.load(f)
    leaf = ((desc.get("stacks") or {}).get(hpc_system) or {}).get("runtime_modules")
    if not leaf:
        return None
    return leaf.get("default") if isinstance(leaf, dict) else leaf


def assert_grx_stack(hpc_system, wax_modules, source):
    """Assert machines.yaml <hpc>.stacks.grx matches waxrun's runtime_modules (cesm.exe is built under the latter).

    A version-less pin on either side matches any version."""
    if not wax_modules:
        raise AssertionError(f"{source} not available; the machines.yaml {hpc_system}.stacks.grx "
                             "cannot be cross-checked against waxrun's runtime_modules")
    jg = JobGen(CONFIG_DIR, machine=hpc_system)
    engage_mods = [ln[len("module load "):] for ln in jg.module_lines(stack="grx")
                   if ln.startswith("module load ")]
    eng = dict(_module_name_version(m) for m in engage_mods)
    wax = dict(_module_name_version(m) for m in wax_modules)
    if set(eng) != set(wax):
        raise AssertionError(
            f"machines.yaml {hpc_system}.stacks.grx != {source}: only in engage "
            f"{sorted(set(eng) - set(wax))}, only in waxrun {sorted(set(wax) - set(eng))}")
    mismatched = {n: (eng[n], wax[n]) for n in eng if eng[n] and wax[n] and eng[n] != wax[n]}
    if mismatched:
        raise AssertionError(
            f"machines.yaml {hpc_system}.stacks.grx != {source}: version pins differ {mismatched}")


def _grx_layout(waxrun, job_name, run_dir):
    return waxrun.layout(job_name, run_dir)


def _voltron_exe(pbs):
    build_dir = pbs.get("kaiju_build_directory") or makeitso.kaiju_build_default("grx")
    return os.path.join(build_dir, "bin", "voltron_mpi.x")


def build_waxrun_payload(engage_options, makeitso_options, lay, coupling_start, stop, delay_h, step_h,
                         args, answers=None, asked=None):
    """Return the --engage payload for waxrun: engage's options plus makeitso's voltron section and the GRX values."""
    sim, pbs, cpl = engage_options["simulation"], engage_options["pbs"], engage_options["coupling"]
    payload = copy.deepcopy(engage_options)
    payload["pbs"]["run_directory"] = lay["rundir"]
    payload["voltron"] = makeitso_options.get("voltron", {})
    payload["layout"] = dict(lay)
    payload["waccmx"] = {"kaijuhome": pbs.get("kaiju_install_directory") or KAIJUHOME,
                         "voltron_exe": _voltron_exe(pbs)}
    spin_s = int(float(cpl["gcm_spin_up_time"]))
    # Raise waxrun's IC-pick bound when the requested spin-up exceeds its default.
    max_h = grx_spinup_max_hours(spin_s)
    if max_h > GRX_SPINUP_MAX_HOURS:
        payload["waccmx"]["spinup"] = {"max_hours": max_h}
    payload["grx"] = {
        "coupling_start": coupling_start.strftime(GRX_DATE_FORMAT),
        "stop": stop.strftime(GRX_DATE_FORMAT),
        "t0": (coupling_start - datetime.timedelta(hours=delay_h)).strftime(GRX_DATE_FORMAT),
        "delaycouplinghrs": delay_h,
        "stephours": step_h,
        "target_start": (coupling_start - datetime.timedelta(seconds=spin_s)).strftime(GRX_DATE_FORMAT),
        "project": pbs.get("account_name"),
        "job_priority": pbs.get("job_priority"),
        "run_directory": lay["rundir"],
    }
    payload["engage_args"] = {"mode": args.mode, "clobber": bool(args.clobber),
                              "recreate_cases": bool(getattr(args, "recreate_cases", False)),
                              "debug": bool(args.debug), "verbose": bool(args.verbose),
                              "tui": bool(getattr(args, "tui", False))}
    if answers is not None:
        payload["prompt_answers"] = answers
    if asked:
        payload["asked_answers"] = asked
    return payload


# waxrun's spinup.max_hours default
GRX_SPINUP_MAX_HOURS = 480


def grx_spinup_max_hours(spin_s):
    """max(480, the requested spin-up in hours, rounded up)."""
    return max(GRX_SPINUP_MAX_HOURS, int(math.ceil(float(spin_s) / 3600.0)))


def grx_spinup_warning(requested_s, spinup_spec):
    """Warning text if waxrun's planned spin-up (set by the IC restart it picked) differs from the request."""
    if not isinstance(spinup_spec, dict) or spinup_spec.get("stop_seconds") in (None, ""):
        return None
    planned = int(spinup_spec["stop_seconds"])
    requested = int(round(float(requested_s)))
    if planned == requested:
        return None
    ic = spinup_spec.get("ic") or {}
    m = re.fullmatch(r"(\d{4}-\d{2}-\d{2})-(\d{5})", str(ic.get("pick_stamp") or ""))
    pick = (f" (start file {m.group(1)}T{int(m.group(2)) // 3600:02d}:{int(m.group(2)) % 3600 // 60:02d}:"
            f"{int(m.group(2)) % 60:02d})" if m else "")
    return (f"WARNING: the WACCM-X spin-up is {planned / 3600:g} h{pick}, not the requested "
            f"{requested / 3600:g} h (coupling.gcm_spin_up_time).")


def build_waxrun_argv(engage_options, makeitso_options, lay, coupling_start, stop, delay_h, step_h,
                      args, answers=None, asked=None):
    """Return waxrun's argv; --interactive only without -wo, and --recreate-cases becomes waxrun's --clobber."""
    payload = build_waxrun_payload(engage_options, makeitso_options, lay, coupling_start, stop,
                                   delay_h, step_h, args, answers, asked)
    argv = ["--engage", json.dumps(payload, default=str), "--mode", args.mode]
    if getattr(args, "waxrun_options_path", None):
        argv += ["-wo", args.waxrun_options_path]
    else:
        argv.append("--interactive")
        if answers is None:
            argv.append("--no-tui")     # no wizard answers: waxrun asks linearly
    if getattr(args, "recreate_cases", False):
        argv.append("--clobber")
    if getattr(args, "dry_run", False):
        argv.append("--dry-run")
    if getattr(args, "skip_build", False):
        argv.append("--skip-build")
    argv += ["--phase", getattr(args, "phase", None) or "both"]
    return argv


def engage_owned_paths(engage_options):
    """Key paths of the makeitso options engage overrides in a coupled run."""
    paths = [("simulation", k) for k in engage_options.get("simulation", {})]
    paths += [("pbs", k) for k in engage_options.get("pbs", {}) if k != "num_segments"]
    return paths + [("gamera", "sim", "runid"), ("voltron", "time", "tFin"), ("gamera", "sim", "H5Grid"),
                    ("gamera", "iPdir", "N"), ("gamera", "jPdir", "N"), ("gamera", "kPdir", "N"),
                    ("voltron", "helpers", "useHelpers"), ("voltron", "helpers", "numHelpers")]


def _get_path(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return None, False
        d = d[k]
    return d, True


def note_engage_owned(before, after, engage_options, source, derived=()):
    """Print one NOTE of the engage-owned -mo values engage replaced or dropped; return their names."""
    noted, changes = [], []
    for path in engage_owned_paths(engage_options):
        name = ".".join(path)
        if name in derived:
            continue
        old, had = _get_path(before, path)
        new, _ = _get_path(after, path)
        if not had or str(old) == str(new):
            continue
        noted.append(name)
        changes.append(f"{name} {old!r} -> {new!r}")
    dropped = [f"pbs.{k}" for k in (before.get("pbs") or {})
               if k not in (after.get("pbs") or {}) and k not in ("select", "select2", "num_segments")]
    if changes or dropped:
        makeitso._note(f"the makeitso options JSON ({source}) takes engage's values for the keys engage "
                       f"owns (change them in the engage options): "
                       + "; ".join(changes + [f"{k} dropped" for k in dropped]) + ".")
    return noted + dropped


def assign_engage_options(makeitso_options, engage_options):
    """Copy engage's simulation values and pbs block (but num_segments) into makeitso options."""
    for k, v in engage_options["simulation"].items():
        makeitso_options["simulation"][k] = v
    makeitso_options["pbs"] = {k: copy.deepcopy(v) for k, v in engage_options["pbs"].items()
                               if k != "num_segments"}
    makeitso_options["gamera"]["sim"]["runid"] = engage_options["simulation"]["job_name"]
    return makeitso_options


def load_makeitso_json(path):
    """Load a -mo JSON; exit if it is in an older makeitso format."""
    with open(path, "r", encoding="utf-8") as f:
        options = json.load(f)
    problem = makeitso.old_format_problem(options)
    if problem:
        raise SystemExit(f"ERROR: {path}: {problem}")
    try:
        makeitso.ensure_sections(options)
    except ValueError as e:
        raise SystemExit(f"ERROR: {path}: {e}") from None
    return options


def patch_makeitso_options(makeitso_options, engage_options, source="-mo"):
    """GRX -mo: apply engage's owned values to loaded makeitso options, as the GRT path does."""
    with open(OPTION_MAKEITSO_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        od = json.load(f)
    before = copy.deepcopy(makeitso_options)
    assign_engage_options(makeitso_options, engage_options)
    makeitso.apply_grid_layout(makeitso_options, od)
    # makeitso's grx branch re-derives these, so they are not edits to note
    note_engage_owned(before, makeitso_options, engage_options, source,
                      derived=("simulation.start_date", "simulation.stop_date", "pbs.run_directory"))
    return makeitso_options


def write_temp_json(directory, name, options):
    """Write options to <directory>/.<name>.<pid>.tmp and return its path; the caller removes it."""
    path = os.path.join(os.path.abspath(directory), f".{name}.{os.getpid()}.tmp")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)
    return path


def link_voltron(run_dir, voltron_exe, warn=True):
    """Symlink <run_dir>/voltron_mpi.x to the absolute voltron executable."""
    voltron_exe = os.path.abspath(voltron_exe)
    link = os.path.join(run_dir, "voltron_mpi.x")
    if os.path.islink(link):
        os.remove(link)
    elif os.path.exists(link):
        print(f"WARNING: {link} is a regular file; left in place (engage would have linked {voltron_exe})")
        return link
    os.symlink(voltron_exe, link)
    if warn and not os.path.exists(voltron_exe):
        print(f"WARNING: voltron_mpi.x -> {voltron_exe} is dangling (build voltron before submitting)")
    return link


def link_lead_in_voltron(run_dir):
    """Link gr_spinup/voltron_mpi.x to ../voltron_mpi.x so the lead-in runs the coupled binary; return it or None."""
    sub = os.path.join(run_dir, makeitso.GR_SPINUP_DIR)
    if not os.path.isdir(sub) or os.path.islink(sub):
        return None
    link = os.path.join(sub, "voltron_mpi.x")
    target = os.path.join(os.pardir, "voltron_mpi.x")
    if os.path.islink(link):
        if os.readlink(link) == target:
            return link
        os.remove(link)
    elif os.path.exists(link):
        print(f"WARNING: {link} is a regular file; left in place (engage would have linked it to "
              f"{target}; the lead-in jobs copy the build's voltron over it)")
        return link
    os.symlink(target, link)
    return link


# A missing executable is a warning (decks may precede the build); a wrong one is refused.
_REPORTED = set()

MISSING_EXE_HINT = {
    "spinup": "the TIE-GCM spin-up executable (the tiegcm_standalone jobs run it): build it with "
              "tiegcmrun.py --onlycompile",
    "coupled": "the coupled TIE-GCM executable (the coupled jobs run it): build it with "
               "tiegcmrun.py --onlycompile --coupling",
    "voltron": "voltron_mpi.x (the jobs run the run directory's link to it): build it in "
               "pbs.kaiju_build_directory",
}


def _report_once(key, message):
    if key in _REPORTED:
        return
    _REPORTED.add(key)
    print(message)


TIEGCM_RUN_ROOT_KEYS = ("parentdir", "execdir", "workdir", "histdir", "modelexe", "coupled_modelexe",
                        "log_file")


def tiegcm_json_relative(options, run_dir):
    """Copy of the tiegcmrun options with run-root paths inside run_dir made relative to it."""
    out = copy.deepcopy(options if isinstance(options, dict) else dict(options))
    data = (out.get("model") or {}).get("data")
    if not isinstance(data, dict):
        return out
    rd = os.path.abspath(run_dir)
    for key in TIEGCM_RUN_ROOT_KEYS:
        v = data.get(key)
        if isinstance(v, str) and os.path.isabs(v):
            data[key] = JobGen.run_relative(v, rd, rd, dot=key in ("modelexe", "coupled_modelexe"))
    return out


def dump_json(path, options):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)


def voltron_preflight(exe, component):
    """Identify voltron_mpi.x by content; warn if missing, exit if it cannot be the right binary."""
    ident = binident.voltron_identity(exe)
    level, msg = binident.voltron_problem(ident, component)
    if level == "refuse":
        raise SystemExit(f"ERROR: {msg}")
    if level == "missing":
        _report_once(("missing", exe), f"WARNING: {exe} does not exist - {MISSING_EXE_HINT['voltron']}")
    return ident


def _grt_expected_grids(engage_options):
    """{'spinup'|'coupled': (horires, vertres, nres_grid)} for the GAMERA grid; Nones if unknown."""
    _import_tiegcmrun()
    try:
        spin, cpl = importlib.import_module("engage_solver").gamres_to_res(
            engage_options["simulation"]["gamera_grid_type"])
    except (KeyError, ValueError):
        return {"spinup": (None, None, None), "coupled": (None, None, None)}
    solver = _tiegcm_misc().resolution_solver
    return {role: (h, solver(h)[0], solver(h)[2]) for role, h in (("spinup", spin), ("coupled", cpl))}


def _tiegcm_exes_of(tiegcm_options_path):
    """{modelexe, coupled_modelexe, zitop} from a -to JSON, or None."""
    if not tiegcm_options_path:
        return None
    try:
        with open(tiegcm_options_path, "r", encoding="utf-8") as f:
            to = json.load(f)
        data = to["model"]["data"]
        exe, cexe = data.get("modelexe"), data.get("coupled_modelexe")
        zitop = (to["model"].get("specification") or {}).get("zitop")
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None
    # Relative paths are relative to the JSON; a bare name is resolved later by tiegcmrun (execdir).
    jdir = os.path.dirname(os.path.abspath(tiegcm_options_path))

    def resolved(p):
        if not isinstance(p, str) or not p:
            return None
        if os.path.isabs(p):
            return p
        if os.path.dirname(p):
            return os.path.normpath(os.path.join(jdir, p))
        return None
    return {"modelexe": resolved(exe), "coupled_modelexe": resolved(cexe), "zitop": zitop}


def grt_to_format_preflight(tiegcm_options_path):
    """Exit if the -to JSON is in an older tiegcmrun format."""
    if not tiegcm_options_path:
        return None
    try:
        with open(tiegcm_options_path, "r", encoding="utf-8") as f:
            options = json.load(f)
    except (OSError, ValueError):
        return None                       # tiegcmrun reports it
    _import_tiegcmrun()
    problem = importlib.import_module("replay_solver").old_format_problem(options)
    if problem:
        raise SystemExit(f"ERROR: {tiegcm_options_path}: {problem}")
    return None


def grt_to_resource_preflight(tiegcm_options_path, engage_options):
    """Exit if the spin-up task count pinned in a -to JSON breaks tiegcmrun's rule (misc.ntask_error)."""
    if not tiegcm_options_path:
        return None
    try:
        with open(tiegcm_options_path, "r", encoding="utf-8") as f:
            res = (json.load(f).get("job") or {}).get("resource") or {}
        n = int(res["select"]) * int(res["mpiprocs"])
    except (OSError, ValueError, AttributeError, KeyError, TypeError):
        return None
    horires = _grt_expected_grids(engage_options)["spinup"][0]
    if horires is None:
        return None
    err = _tiegcm_misc().ntask_error(n, horires, None,
                                     where=f"job.resource of the TIE-GCM options JSON (the spin-up ranks "
                                           f"at {horires} deg)")
    if err:
        raise SystemExit("ERROR: the TIE-GCM jobs of this coupled run would be rejected (nothing was "
                         "written): " + ANSI_RE.sub("", err))
    return None


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

# Errors a tool raises that engage prints as one ERROR line.
TOOL_ERRORS = (ValueError, FileExistsError, RuntimeError)


def tool_error(tool, path, exc):
    """Format a tool's error as one line: ERROR: <tool> options (<file>): <message>."""
    where = f" options ({os.path.basename(path)})" if path else ""
    return f"ERROR: {tool}{where}: {ANSI_RE.sub('', str(exc)).strip()}"


def grt_executable_preflight(modelexe, coupled_modelexe, engage_options, zitop=None):
    """Check the spin-up, coupled TIE-GCM and voltron executables by content; return their identities.

    modelexe must be a standalone build at the spin-up grid, coupled_modelexe a -DGAMERA build at the coupled grid."""
    misc = _tiegcm_misc()
    grids = _grt_expected_grids(engage_options)
    tiegcm, refusals = {}, []
    for role, key, exe in (("spinup", "modelexe", modelexe), ("coupled", "coupled_modelexe", coupled_modelexe)):
        ident = misc.exe_identity(exe, zitop)
        tiegcm[role] = ident
        if not exe:
            continue
        if not ident["exists"]:
            _report_once(("missing", exe), f"WARNING: {exe} does not exist - {MISSING_EXE_HINT[role]}")
            continue
        not_elf, problems = misc.exe_problems(ident, f"model.data.{key}", role == "coupled",
                                              *grids[role], zitop)
        refusals += ([not_elf] if not_elf else []) + problems
    if refusals:
        raise SystemExit("ERROR: " + "\nERROR: ".join(refusals))
    voltron = voltron_preflight(_voltron_exe(engage_options["pbs"]), "grt")
    return {"voltron": voltron, "tiegcm": tiegcm}


def _grx_waccmx_pe(coupled_spec):
    """(nwnodes, ntasks_atm) of the coupled WACCM-X case: from the spec, the caseroot, or defaults."""
    pe = coupled_spec.get("pe") or {}
    ntasks = coupled_spec.get("ntasks_atm") or pe.get("ntasks_atm") or pe.get("NTASKS_ATM")
    nwnodes = coupled_spec.get("nwnodes") or pe.get("nodes") or pe.get("nwnodes")
    tasks_per_node = pe.get("tasks_per_node")
    caseroot = coupled_spec.get("caseroot")
    if ntasks is None and caseroot and os.path.isfile(os.path.join(caseroot, "env_mach_pes.xml")):
        case_solver = _import_waxrun_module("case_solver")
        q = case_solver.xmlquery(caseroot, ["NTASKS_ATM", "MAX_MPITASKS_PER_NODE"])
        ntasks, tasks_per_node = q.get("NTASKS_ATM"), q.get("MAX_MPITASKS_PER_NODE")
    if ntasks is None:
        seg = _import_segment_solver()
        nwnodes, tasks_per_node = seg.DEFAULTS["nwnodes"], seg.DEFAULTS["nwcpus"]
        ntasks = nwnodes * tasks_per_node
        print(f"WARNING: coupled spec carries no NTASKS_ATM; using the segment_solver default "
              f"{nwnodes} x {tasks_per_node} = {ntasks} WACCM-X ranks")
    ntasks = int(ntasks)
    if nwnodes is None:
        if tasks_per_node is None:
            tasks_per_node = _import_segment_solver().DEFAULTS["nwcpus"]
        nwnodes = math.ceil(ntasks / int(tasks_per_node))
    nwnodes = int(nwnodes)
    if ntasks % nwnodes:
        raise ValueError(f"NTASKS_ATM={ntasks} is not a multiple of nwnodes={nwnodes}")
    return nwnodes, ntasks


def build_grx_pbs_options(gr_options: dict, makeitso_options: dict, coupled_spec: dict) -> dict:
    """GRX analog of build_grt_pbs_options: GR layout from makeitso, WACCM-X from waxrun's coupled spec."""
    options = copy.deepcopy(gr_options)
    p = options["pbs"]

    p["gamera_nodes"] = makeitso_options["pbs"]["select"]
    p["gamera_mpiprocs"] = makeitso_options["pbs"]["mpiprocs"]
    p["gamera_ompthreads"] = makeitso_options["pbs"]["ompthreads"]
    p["voltron_nodes"] = makeitso_options["pbs"]["select2"]
    p["voltron_mpiprocs"] = makeitso_options["pbs"]["helper_mpiprocs"]
    p["voltron_ompthreads"] = makeitso_options["pbs"]["helper_ompthreads"]
    p["voltron_mpiranks"] = (int(p["gamera_nodes"]) * int(p["gamera_mpiprocs"])
                             + int(p["voltron_nodes"]) * int(p["voltron_mpiprocs"]))

    hpc = options["simulation"]["hpc_system"]
    if hpc != "derecho":
        raise NotImplementedError(
            f"Coupled GRX runs currently support only derecho (CESM/WACCM-X MPMD launch). "
            f"hpc_system={hpc!r} has no GRX launch dispatch — add a branch or run standalone.")
    nwnodes, ntasks = _grx_waccmx_pe(coupled_spec)
    p["waccmx_nodes"] = nwnodes
    p["waccmx_mpiprocs"] = ntasks // nwnodes
    p["waccmx_mpiranks"] = ntasks

    cesm_exe = coupled_spec["cesm_exe"]
    if not os.path.isabs(cesm_exe):
        raise ValueError(f"coupled spec cesm_exe must be absolute: {cesm_exe!r}")
    p["cesm_exe"] = cesm_exe
    p["kaiju_install_directory"] = (p.get("kaiju_install_directory") or coupled_spec.get("kaijuhome")
                                    or KAIJUHOME)
    p["pin"] = "${KAIJUHOME}/scripts/preproc/pinCpuCores.sh"
    p["product_path"] = coupled_spec.get("ncdata")
    p["cesm_rundir"] = coupled_spec.get("rundir")
    return options


def _grx_header(job_name, seg_name, phase):
    return [f"# {seg_name}: GRX {phase} segment of {job_name}."]


def _grx_preflight(pbs, job_name, nres_str, phase="both"):
    """Preflight lines of <job>-01: lead-in products exist and postdate this attempt, executables resolve,
    the GR restarts are linked in from gr_spinup/, and nuopc.runconfig is reset to startup
    (<job>-02 flips it to continue, so a resubmitted <job>-01 would otherwise start in continue mode)."""
    run_dir = pbs["cesm_rundir"]

    def here(path):
        return JobGen.spell(path, run_dir, run_dir)

    product = here(pbs["product_path"])
    cesm_exe = here(pbs["cesm_exe"])
    stamp = here(os.path.join(run_dir, f"{job_name}_grx.stamp"))
    jobids = here(os.path.join(run_dir, f"{job_name}_grx.jobids"))
    lead = makeitso.GR_SPINUP_DIR
    required = [f"{job_name}.{fam}.Res.{nres_str}.h5" for fam in binident.RES_FAMILIES
                if fam != "gam" and fam not in binident.RES_OPTIONAL]
    if phase == "coupled":
        product_age = [f"# --phase coupled: the spinup product comes from the earlier --phase spinup attempt"]
    else:
        product_age = [f'[ "{product}" -nt "$stamp" ] || fail "spinup product predates this attempt: {product}"']
    return [
        "# ---- preflight: the lead-in output must exist and be from this attempt",
        f"stamp={stamp}",
        f"g1=${{GR_WARMUP_JOBID:-$(sed -n 's/^gr_warmup=//p' {jobids} 2>/dev/null)}}",
        'fail() { echo "GATE: $1"; exit 90; }',
        f'[ -e "{product}" ] || fail "spinup product missing: {product}"',
        *product_age,
        f'[ -n "$g1" ] || fail "{job_name}-WARMUP-01 job id unknown (GR_WARMUP_JOBID unset and no gr_warmup= line in {jobids})"',
        f'[ -e "{lead}/WARMUP_DONE.$g1" ] || fail "{lead}/WARMUP_DONE.$g1 missing (GR warm-up did not finish in this attempt)"',
        # a missing dk.Res only cold-starts DragonKing, so it is a warning
        f"# ---- GR restart hand-off: link restart index {nres_str} from {lead}/",
        f'pull() {{ s="{lead}/$1"; d="./$1"',
        '  if [ -e "$s" ]; then [ "$s" -ef "$d" ] || ln -f "$s" "$d" 2>/dev/null || ln -sf "$s" "$d"',
        '  elif [ -e "$d" ]; then :',
        '  else return 1; fi; }',
        f'src() {{ if [ -e "{lead}/$1" ]; then echo "{lead}/$1"; else echo "./$1"; fi; }}',
        f'gam=$(for s in {lead}/{job_name}*.gam.Res.{nres_str}.h5 ./{job_name}*.gam.Res.{nres_str}.h5; do '
        'if [ -e "$s" ]; then echo "${s##*/}"; fi; done | sort -u)',
        f'[ -n "$gam" ] || fail "missing {job_name}*.gam.Res.{nres_str}.h5 (neither in {lead}/ nor here)"',
        "for f in $gam " + " ".join(required) + "; do",
        '  s=$(src "$f")',
        f'  [ -e "$s" ] || fail "missing $f (neither in {lead}/ nor here)"',
        '  [ "$s" -nt "$stamp" ] || fail "$s predates this attempt"',
        '  pull "$f" || fail "cannot link $s here"',
        "done",
        f"f={job_name}.dk.Res.{nres_str}.h5",
        's=$(src "$f")',
        'if [ -e "$s" ]; then [ "$s" -nt "$stamp" ] || fail "$s predates this attempt"; '
        'pull "$f" || fail "cannot link $s here"; '
        'else echo "WARNING: $f missing: DragonKing cold-starts"; fi',
        "# ---- end GR restart hand-off",
        f'[ -x "{cesm_exe}" ] || fail "cesm.exe missing or not executable: {cesm_exe}"',
        f"[ \"$(ldd {cesm_exe} | grep -c 'not found')\" -eq 0 ] || "
        "{ echo \"GATE: cesm.exe has unresolved libs under the runtime modules\"; exit 91; }",
        '[ -x ./voltron_mpi.x ] || fail "voltron_mpi.x missing or not executable: ./voltron_mpi.x"',
        "[ \"$(ldd ./voltron_mpi.x | grep -c 'not found')\" -eq 0 ] || "
        "{ echo \"GATE: voltron_mpi.x has unresolved libs under the runtime modules\"; exit 91; }",
        "# back to startup mode (a resubmitted first segment)",
        r"sed -i 's/^\(\s*start_type\s*=\s*\)continue\b/\1startup/' nuopc.runconfig",
        "grep -Eq '^\\s*start_type\\s*=\\s*startup' nuopc.runconfig || "
        "fail \"nuopc.runconfig is not in startup mode; run ./preview_namelists in the coupled caseroot\"",
        "# ---- end preflight",
    ]


def render_grx_pbs(options: dict, ii: int, ic: int, iu: int, coupling_start, run_phase="both") -> str:
    """Render the job of GRX segment ii: 0 = SPINUP, 1..ic = WARMUP (voltron only), > ic = coupled MPMD.

    Coupled segments self-chain on rc==0; the lead-in jobs cd to gr_spinup/ by absolute path."""
    seg = _import_segment_solver()
    sim, pbs, mage = options["simulation"], options["pbs"], options["mage"]
    step, delay = int(mage["stephours"]), int(mage["delaycouplinghrs"])
    jg = JobGen(CONFIG_DIR, machine=sim["hpc_system"], node_type=sim.get("node_type"))
    coupled = ii > ic
    lead_in = not coupled

    job_name = sim["job_name"]
    seg_name = seg.segment_id(job_name, ii, ic)
    volt_xml = f"{seg_name}.xml"
    nv = pbs["voltron_mpiranks"]
    pin = pbs["pin"]
    phase = "fully coupled" if coupled else ("spinup" if ii == 0 else "coupling-delay")

    gamera = {"mpi_tasks": int(pbs["gamera_nodes"]) * int(pbs["gamera_mpiprocs"]),
              "tasks_per_node": int(pbs["gamera_mpiprocs"]), "threads": int(pbs["gamera_ompthreads"])}
    voltron = {"mpi_tasks": int(pbs["voltron_nodes"]) * int(pbs["voltron_mpiprocs"]),
               "tasks_per_node": int(pbs["voltron_mpiprocs"]), "threads": int(pbs["voltron_ompthreads"])}
    if coupled:
        waccmx = {"mpi_tasks": int(pbs["waccmx_mpiranks"]),
                  "tasks_per_node": int(pbs["waccmx_mpiprocs"]), "threads": 1}
        task_groups = [waccmx, gamera, voltron]
    else:
        task_groups = [gamera, voltron]

    L = [jg.shebang()]
    L += _grx_header(job_name, seg_name, phase)
    L += jg.directive_lines({"job_name": seg_name})
    L += jg.directive_lines({"join_out": True})
    L += jg.directive_lines({"queue": pbs.get("queue")})
    L += jg.directive_lines({"account": pbs.get("account_name")})
    resources = jg.resource_lines(task_groups, pbs["walltime"], queue=pbs.get("queue"))
    L += resources[:1]
    L += jg.directive_lines({"priority": pbs.get("job_priority")})
    L += resources[1:]

    L += ["", f"export KAIJUHOME={pbs['kaiju_install_directory']}", ""]
    tmpdir = jg.stack_option("grx", "tmpdir")
    if tmpdir:
        L += [f"export TMPDIR={tmpdir}", "mkdir -p $TMPDIR"]
    L += jg.module_lines(stack="grx")
    run_dir = pbs["cesm_rundir"]
    if lead_in:
        workdir_cd = f"cd {shlex.quote(os.path.join(os.path.abspath(run_dir), makeitso.GR_SPINUP_DIR))}"
    else:
        workdir_cd = jg.sched("workdir_cd")
    L += [f"export OMP_NUM_THREADS={pbs['voltron_ompthreads']}",
          workdir_cd,
          "export KMP_STACKSIZE=512M",
          "export JNUM=${PBS_ARRAY_INDEX:-0}",
          "ulimit -s unlimited",
          "module list",
          "hostname",
          "date"]

    log = f"voltron_mpi.x-{seg_name}.${{PBS_JOBID}}.out"
    out = f"> {log}"
    probe_lines = []
    if "voltron_sha256" in pbs and ii >= 1:
        probe_lines = binident.job_probe_lines("./voltron_mpi.x", pbs["voltron_sha256"])
    if coupled:
        if ii == ic + 1:
            nres = seg.nres(ii, ic, delay, step, dtRes=mage.get("dtRes", 1800.0))
            L += _grx_preflight(pbs, job_name, f"{nres:05d}", run_phase)
        if ii == ic + 2:
            L += [r"sed -i 's/\bstartup\b/continue/g' nuopc.runconfig",
                  r"sed -i 's/rpointer\.cpl\..*/rpointer.cpl/' nuopc.runconfig"]
        if ii >= ic + 2:
            # the first coupled segment is a CESM startup and reads no rpointer
            suf = seg.rpointer_suffix(ii, ic, coupling_start, step)
            L += [f"ln -sf rpointer.ice.{suf} rpointer.ice",
                  f"ln -sf rpointer.cpl.{suf} rpointer.cpl"]
        L += grx_log_name_lines(ii - ic)
        L += probe_lines
        cesm_exe = JobGen.spell(pbs["cesm_exe"], run_dir, run_dir)
        L += [f"mpiexec --label --line-buffer -n {pbs['waccmx_mpiranks']} {cesm_exe} : "
              f"-n {nv} {pin} ./voltron_mpi.x {volt_xml} {out}"]
    else:
        L += probe_lines
        L += [f'echo "Running {seg_name}"',
              f"mpiexec -n {nv} {pin} ./voltron_mpi.x {volt_xml} {out}"]

    if ii == 0:
        L += ["rc=$?",
              'if [ $rc -ne 0 ]; then echo "segment failed rc=$rc"; fi',
              "exit $rc"]
    elif ii == ic:
        L += ["rc=$?",
              'if [ $rc -eq 0 ]; then touch WARMUP_DONE.${PBS_JOBID}; '
              f'echo "GR warm-up done; {seg.segment_id(job_name, ic + 1, ic)} waits on this marker and the '
              f'WACCM-X spinup"; '
              'else echo "segment failed rc=$rc; NOT chaining"; fi',
              "exit $rc"]
    else:
        L += ["rc=$?"]
        L += grx_mixcpl_lines(log)
        if ii < iu - 1:
            nxt = f"{seg.segment_id(job_name, ii + 1, ic)}.pbs"
            jobids = JobGen.spell(os.path.join(run_dir, f"{job_name}_grx.jobids"), run_dir, run_dir)
            L += [f'if [ $rc -eq 0 ]; then id=$(qsub {nxt}) && {{ echo "$id"; echo "{nxt}=$id" >> {jobids}; }}; '
                  'else echo "segment failed rc=$rc; NOT chaining"; fi']
        else:
            L += ['if [ $rc -ne 0 ]; then echo "segment failed rc=$rc"; fi']
        L += ["exit $rc"]
    return "\n".join(L) + "\n"


def grx_log_name_lines(k):
    """Job lines that rename CESM logs to <comp>.log.seg<k>.<jobid> in nuopc.runconfig.

    buildnml fixes one log name per case, so every segment would otherwise append to the same logs."""
    return [f"# ---- CESM log names: <comp>.log.seg{k:02d}.<jobid>",
            r'sed -i -E "s/^(\s*logfile\s*=\s*[A-Za-z0-9_]+\.log)\b.*/\1.seg'
            + f'{k:02d}' + r'.${PBS_JOBID}/" nuopc.runconfig',
            f'echo "CESM logs of this segment: *.log.seg{k:02d}.${{PBS_JOBID}}"',
            "# ---- end CESM log names"]


def grx_mixcpl_lines(out):
    """Job lines reporting voltron's MIXCPL WARNING count (an unknown WACCM-X variable skipped); rc unchanged."""
    return [f"n=$(grep -c 'MIXCPL WARNING' {out} 2>/dev/null); "
            'if [ "${n:-0}" -gt 0 ]; then echo "MIXCPL WARNING seen $n times: voltron skipped a WACCM-X '
            'variable it does not know (see the voltron log)"; fi']


def create_grx_pbs_scripts(gr_options: dict, makeitso_options: dict, coupled_spec: dict,
                           t0, t1, step_h, delay_h, run_dir=None, phase="both", voltron=None):
    """Write one job per GRX segment (lead-in jobs in gr_spinup/, coupled ones in the run directory); return the paths."""
    seg = _import_segment_solver()
    options = build_grx_pbs_options(gr_options, makeitso_options, coupled_spec)
    if voltron is None:
        voltron = binident.voltron_identity(coupled_spec.get("voltron_exe") or _voltron_exe(options["pbs"]))
    options["pbs"]["voltron_sha256"] = voltron.get("sha256")
    step, delay = int(step_h), int(delay_h)
    # must match the decks' cadence so the preflight checks the restart <job>-01 resumes from
    dtRes = makeitso_options.get("voltron", {}).get("restart", {}).get("dtRes", 1800.0)
    options["mage"] = {"stephours": step, "delaycouplinghrs": delay, "dtRes": float(dtRes)}
    t0 = _grx_datetime(t0)
    ic, iu, maxhours = seg.compute_segments(t0, _grx_datetime(t1), step, delay)
    coupling_start = t0 + datetime.timedelta(hours=delay)

    run_dir = os.path.abspath(run_dir or options["pbs"].get("cesm_rundir")
                              or options["pbs"].get("run_directory", "."))
    options["pbs"]["cesm_rundir"] = run_dir
    os.makedirs(run_dir, exist_ok=True)
    pbs_scripts = []
    print(f"Creating {iu} GRX segment PBS scripts (ic={ic}, maxhours={maxhours})")
    for ii in range(0, iu):
        content = render_grx_pbs(copy.deepcopy(options), ii, ic, iu, coupling_start, phase)
        where = makeitso.lead_in_directory(run_dir) if ii <= ic else run_dir
        path = os.path.join(where, f"{seg.segment_id(options['simulation']['job_name'], ii, ic)}.pbs")
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        pbs_scripts.append(path)
    return pbs_scripts


def _nudging_scripts(coupled_spec):
    """The nudging generator job scripts waxrun wrote (each touches <script>.OK on success)."""
    return [str(p) for p in coupled_spec.get("nudging_scripts") or [] if p]


def write_grx_dag(run_dir, job_name, coupled_spec, spinup_submit_scripts, pbs_files, lay,
                  clobber=False, phase="both"):
    """Write <run_dir>/<job>_grx.sh, the GRX submit script; return its path.

    nudging -> WACCM-X spinup (case.submit chain) in parallel with SPINUP -> WARMUP-01; <job>-01 afterok both.
    """
    run_dir = os.path.abspath(run_dir)
    if phase not in ("spinup", "coupled", "both"):
        raise ValueError(f"write_grx_dag: unknown phase {phase!r}")
    if len(pbs_files) < 3:
        raise ValueError(f"a GRX run needs at least 3 jobs (SPINUP, WARMUP-01, 01); got {len(pbs_files)}")

    def here(path):
        return JobGen.spell(os.path.join(run_dir, path) if not os.path.isabs(path) else path,
                            run_dir, run_dir)

    def job_dir(p):
        """(script directory relative to run_dir, or None for run_dir itself; script name)."""
        ap = os.path.join(run_dir, p)
        sub = JobGen.run_relative(os.path.dirname(ap), run_dir, run_dir, dot=False)
        return (None if sub == os.curdir else sub), os.path.basename(p)

    (s_spin, p_spin), (s_warm, p_warm), (s_first, p_first) = (job_dir(p) for p in pbs_files[:3])
    lead = makeitso.GR_SPINUP_DIR
    gen_scripts = [JobGen.run_relative(os.path.join(run_dir, g), run_dir, run_dir, dot=False)
                   for g in _nudging_scripts(coupled_spec)]
    chunks = len(spinup_submit_scripts) if spinup_submit_scripts else 1
    capture = r"grep -oP 'Submitted job case\.run with id \K\d+'"
    stamp, jobids = f"{job_name}_grx.stamp", f"{job_name}_grx.jobids"
    do_spinup, do_coupled = phase != "coupled", phase != "spinup"
    cesm_exe = coupled_spec.get("cesm_exe")
    cesm_exe = here(cesm_exe) if cesm_exe else cesm_exe
    product = coupled_spec.get("ncdata") or coupled_spec.get("spinup_product")
    product = here(product) if product else product

    def qsub_from(sub, script, opts=""):
        """qsub the script from its own directory, so $PBS_O_WORKDIR and its output land there."""
        opts = f"{opts} " if opts else ""
        if not sub:
            return f"qsub {opts}{script}"
        return f"cd {shlex.quote(JobGen.dot_relative(sub))} && qsub {opts}{shlex.quote(JobGen.dot_relative(script))}"

    L = ["#!/usr/bin/env bash",
         f"# Submit the {job_name} GRX run" + (f" (--phase {phase})" if phase != "both" else "")
         + ". Run from anywhere.",
         "set -euo pipefail",
         "unset CDPATH",                 # an exported CDPATH would redirect the relative cds below
         'cd "$(dirname "$0")"', "",
         "# checks"]
    if do_coupled:
        L += [f'[ -x "{cesm_exe}" ] || {{ echo "REFUSING: cesm.exe missing; build the coupled case first: '
              f'{cesm_exe}"; exit 3; }}']
        if not do_spinup:
            L += [f'[ -e "{product}" ] || {{ echo "REFUSING: --phase coupled needs the spinup product: {product}"; exit 3; }}']
    restart_globs = ([f"{job_name}*.gam.Res.*.h5"]
                     + [f"{job_name}.{fam}.Res.*.h5" for fam in binident.RES_FAMILIES if fam != "gam"])
    L += [f': "${{CLOBBER:={1 if clobber else 0}}}"',
          f"stale=$(ls -1 WARMUP_DONE.* " + " ".join(restart_globs)
          + f" {lead}/WARMUP_DONE.* " + " ".join(f"{lead}/{g}" for g in restart_globs)
          + f" rpointer.* {job_name}.*.r.*.nc {job_name}.cam.r[sh]*.nc 2>/dev/null || true)",
          'if [ -n "$stale" ]; then',
          '  if [ "$CLOBBER" = "1" ]; then',
          '    echo "Removing leftovers of a previous attempt."; rm -f $stale',
          "  else",
          '    echo "REFUSING: leftovers of a previous attempt (run CLOBBER=1 $0 to remove them):"; echo "$stale"; exit 2',
          "  fi",
          "fi",
          f"touch {stamp}", ""]
    nud_ids = []
    nudged = bool((coupled_spec.get("nudging") or {}).get("enabled", True))
    if do_spinup and nudged:
        L += ["# nudging data", 'nud=""']
    if do_spinup:
        for k, g in enumerate(gen_scripts, start=1):
            var = "nud" if len(gen_scripts) == 1 else f"nud{k}"
            nud_ids.append(var)
            gsub, gname = os.path.split(g)
            L += [f"rm -f {g}.OK", f"{var}=$({qsub_from(gsub, gname)})", f'echo "{var}=${var}"']
        if len(gen_scripts) > 1:
            L += ['nud="' + ":".join(f"${v}" for v in nud_ids) + '"']
        L += ([""] if nudged else []) + [
              "# WACCM-X spinup",
              f"cd {shlex.quote(here(lay['caseroot_spinup']))}"]
        if not gen_scripts:
            L += [f"sp1=$(./case.submit | {capture})"]
        else:
            L += ['if [ -n "$nud" ]; then',
                  f'  sp1=$(./case.submit --prereq "$nud" | {capture})',
                  "else",
                  f"  sp1=$(./case.submit | {capture})",
                  "fi"]
        L += ["sp_last=$sp1"]
        for k in range(2, chunks + 1):
            L += [f'sp{k}=$(./case.submit --prereq "$sp_last" | {capture})',
                  f"sp_last=$sp{k}"]
        L += [f"# resubmit one chunk by hand: bash {here(s_)}" for s_ in spinup_submit_scripts or []]
        L += ["cd - >/dev/null", ""]
    if do_coupled:
        L += ["# GR lead-in",
              f"g0=$({qsub_from(s_spin, p_spin)})",
              f"g1=$({qsub_from(s_warm, p_warm, '-W depend=afterok:$g0')})", "",
              "# coupled segments"]
        if do_spinup:
            L += [f"c=$({qsub_from(s_first, p_first, '-v GR_WARMUP_JOBID=$g1 -W depend=afterok:$sp_last:$g1')})", ""]
        else:
            L += [f"c=$({qsub_from(s_first, p_first, '-v GR_WARMUP_JOBID=$g1 -W depend=afterok:$g1')})", ""]
    ids = [("nudging", '"$nud"' if do_spinup and nudged else '""'),
           ("spinup", '"$sp_last"' if do_spinup else '""'),
           ("gr_spinup", '"$g0"' if do_coupled else '""'),
           ("gr_warmup", '"$g1"' if do_coupled else '""'),
           ("coupled", '"$c"' if do_coupled else '""')]
    L += ["printf '" + "".join(f"{k}=%s\\n" for k, _ in ids) + "' \\",
          "       " + " ".join(v for _, v in ids) + f" > {jobids}",
          f"cat {jobids}"]
    # The heredocs below are unquoted for ${here} only; their text must hold no other $, ` or \.
    L += ["here=$(printf '%q' \"$PWD\")"]
    if do_coupled:
        L += ["cat <<EOF",
              "To resubmit the coupled chain after a fix (remove CESM leftovers first, or rerun with CLOBBER=1):",
              f"    cd ${{here}} && qsub {p_first if not s_first else os.path.join(s_first, p_first)}",
              "EOF"]
    else:
        L += ["cat <<EOF",
              "Spinup only. Once the spinup product exists, build the coupled case and run:",
              "    cd ${here} && engage.py --component grx --phase coupled --clobber -eo engage_parameters.json "
              "-mo makeitso_parameters.json -wo waxrun_parameters.json",
              "EOF"]
    dag = os.path.join(run_dir, f"{job_name}_grx.sh")
    with open(dag, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    os.chmod(dag, 0o755)
    return dag


def grt_machines(option_engage_descriptions=None):
    """Machines from machines.yaml that have an engage pbs block, in file order."""
    if option_engage_descriptions is None:
        option_engage_descriptions = grx_options.load_engage_descriptions(OPTION_ENGAGE_DESCRIPTIONS_FILE)
    return [m for m in makeitso.machine_names() if m in option_engage_descriptions["pbs"]]


# Other machines are offered but labelled experimental.
GRT_VERIFIED_MACHINES = ("derecho", "aitken")


def grt_machine_warning(hpc_system):
    """Warning text for a coupled GRT run on an experimental machine, else None."""
    if hpc_system in GRT_VERIFIED_MACHINES:
        return None
    return (f"{makeitso.YELLOW}WARNING: a coupled GRT run on {hpc_system} is EXPERIMENTAL: it has never "
            f"run there, and the job scripts (from config/machines.yaml) are unchecked on that "
            f"machine.{makeitso.RESET}")


def set_grt_machine_choices(hpc_description, option_engage_descriptions=None):
    """Restrict the hpc_system description to grt_machines(); keep the host default if valid."""
    machines = grt_machines(option_engage_descriptions)
    hpc_description["valids"] = machines
    if hpc_description.get("default") not in machines:
        hpc_description["default"] = machines[0] if machines else None
    experimental = [m for m in machines if m not in GRT_VERIFIED_MACHINES]
    if experimental:
        hpc_description["NOTE"] = (f"coupled GRT: derecho and aitken are the supported machines; "
                                   f"{', '.join(experimental)}: experimental (never run coupled there).")
    return hpc_description


def set_grid_choices(grid_description, component="grt"):
    """Restrict the grid-type description to makeitso.RUN_GRID_TYPES[component]."""
    grids = list(makeitso.RUN_GRID_TYPES[component])
    grid_description["valids"] = grids
    if grid_description.get("default") not in grids:
        grid_description["default"] = grids[0]
    return grid_description


def account_required(hpc_system):
    """True if the machine's jobs need an account line (machines.yaml account_required)."""
    return makeitso.account_required(hpc_system, CONFIG_DIR)


def settle_engage_account(engage_options, label="the engage options JSON"):
    """Exit if pbs.account_name is missing where required; drop it where the machine takes none."""
    try:
        makeitso.settle_account(engage_options, label, config_dir=CONFIG_DIR)
    except ValueError as e:
        raise SystemExit(f"ERROR: {e}") from None
    return engage_options


def engage_machine_hint(engage_options, platform_od, label="the engage options JSON"):
    """NOTE (suggesting --rederive job) if a replayed pbs block looks made for another machine; return the evidence."""
    sim, pbs = engage_options["simulation"], engage_options["pbs"]
    hpc = str(sim["hpc_system"])
    try:
        evidence = makeitso.machine_block_evidence(pbs, platform_od, hpc)
    except ValueError as e:
        raise SystemExit(f"ERROR: {label}: {e}") from None
    if evidence:
        makeitso._note(f"{label}: the pbs block looks made for another machine than simulation.hpc_system="
                       f"{hpc!r} ({'; '.join(evidence)}): pass --rederive job to follow the machine.")
    return evidence


def _raw_makeitso_descriptions():
    with open(OPTION_MAKEITSO_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        return json.load(f)


def grt_min_nodes(engage_options, option_makeitso_descriptions):
    """GR node count from the grid's default decomposition, plus 1 + num_helpers; None if unknown."""
    sim, pbs = engage_options["simulation"], engage_options["pbs"]
    grid = str(sim.get("gamera_grid_type"))
    try:
        g = option_makeitso_descriptions["gamera"]
        ranks = 1
        for axis in ("iPdir", "jPdir", "kPdir"):
            ranks *= int(g[axis]["N"]["default"][grid])
        return math.ceil(ranks / int(pbs["mpiprocs"])) + 1 + int(pbs["num_helpers"])
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return None


def coupled_segment_seconds(sim, component):
    """Longest coupled segment in seconds (GRT: clamped to the run window)."""
    seg = sim.get("segment_duration")
    if component != "grt":
        return seg
    try:
        window = (makeitso.parse_run_datetime(sim["stop_date"], "simulation.stop_date")
                  - makeitso.parse_run_datetime(sim["start_date"], "simulation.start_date")).total_seconds()
        if str(sim.get("use_segments", "Y")).upper() != "Y":
            return window
        return min(float(seg), window)
    except (KeyError, TypeError, ValueError):
        return seg


def coupled_queue_error(engage_options, component, min_nodes=None):
    """Error if the coupled segments do not fit the queue's walltime or ncpus limits, else None."""
    sim, pbs = engage_options["simulation"], engage_options["pbs"]
    hpc, queue = sim.get("hpc_system"), pbs.get("queue")
    grid = sim.get("gamera_grid_type")
    grx_options.announce(True)          # the wizard (if any) has returned
    err = grx_options.walltime_error(None, pbs.get("walltime"), hpc, queue, component, grid)
    if err:
        return f"ERROR: {err}"
    if component == "grt":
        err = grx_options.grt_window_error(sim.get("start_date"), sim.get("stop_date"))
        if not err and makeitso.yn(sim.get("use_segments", "Y")) == "Y":
            err = grx_options.segment_hours_error(
                "grt", hpc, queue, pbs.get("walltime"), grid, sim.get("start_date"), sim.get("stop_date"),
                coupled_segment_seconds(sim, component))
        if err:
            return f"ERROR: {err}"
    seg = coupled_segment_seconds(sim, component)
    err = grx_options.walltime_error(seg, pbs.get("walltime"), hpc, queue, component, grid)
    if err:
        return f"ERROR: {err}"
    if min_nodes:
        label = "GRX" if component == "grx" else "GRT"
        err = grx_options.ncpus_error(hpc, queue, min_nodes,
                                      f"the coupled {label} job (at least its {min_nodes} GAMERA + voltron nodes)")
        if err:
            return f"ERROR: {err}"
    return None


def makeitso_rules_error(makeitso_options, component, coupling=None):
    """makeitso.option_rules_error for a replayed -mo JSON, or None."""
    args = {"grx": True} if component == "grx" else {"coupling": coupling or {}}
    err = makeitso.option_rules_error(makeitso_options, args)
    return f"ERROR: the makeitso options JSON: {err}" if err else None


# Keys from older engage JSONs, refused with what replaces them.
OLD_FORMAT_KEYS = {
    ("provenance",): "no longer written",
    ("coupling", "doGCM"): "every coupled GRT deck runs doGCM=T; voltron without a GCM is --component gr",
    ("coupling", "root_directory"): "the run lives in pbs.run_directory",
    ("pbs", "select"): "derived from the GAMERA decomposition",
    ("pbs", "ncpus"): "from config/machines.yaml",
    ("pbs", "other"): "use simulation.node_type",
    ("pbs", "modules"): "from config/machines.yaml",
    ("pbs", "moduledir"): "from config/machines.yaml",
    ("pbs", "environment_variables"): "from config/machines.yaml",
    ("pbs", "mpiexec_command"): "from config/machines.yaml",
}


def old_format_problem(engage_options):
    """Error message if the engage JSON carries OLD_FORMAT_KEYS, else None."""
    old = []
    for path, what in OLD_FORMAT_KEYS.items():
        node = engage_options
        for k in path[:-1]:
            node = node.get(k) if isinstance(node, dict) else None
        if isinstance(node, dict) and path[-1] in node:
            old.append(f"{'.'.join(path)} ({what})")
    if not old:
        return None
    return (f"the engage options JSON was written by an older engage: {'; '.join(old)}. Remove "
            f"{'it' if len(old) == 1 else 'them'} from the JSON, or generate a new one.")


# -eo keys with no default (GRT also requires the grid type).
EO_REQUIRED = {"simulation": ("start_date", "stop_date")}


def fill_engage_options(engage_options, option_makeitso_descriptions, option_engage_descriptions,
                        component="grt", filled=None):
    """Fill omitted -eo entries with the prompt defaults and NOTE them; exit on missing required keys.

    With a `filled` list, the filled names are appended to it and the NOTE is left to the caller."""
    sim = engage_options.setdefault("simulation", {})
    pbs = engage_options.setdefault("pbs", {})
    cpl = engage_options.setdefault("coupling", {})
    required = [(sec, k) for sec, keys in EO_REQUIRED.items() for k in keys]
    if component != "grx":
        required.append(("simulation", "gamera_grid_type"))
    missing = [f"{sec}.{k}" for sec, k in required if engage_options[sec].get(k) in (None, "")]
    names = []
    sod = option_makeitso_descriptions["simulation"]
    if sim.get("hpc_system") in (None, ""):
        if component == "grx":
            sim["hpc_system"] = grx_options.grx_defaults()["simulation.hpc_system"]
        elif sod["hpc_system"].get("default") in grt_machines(option_engage_descriptions):
            sim["hpc_system"] = sod["hpc_system"]["default"]
        else:
            missing.append("simulation.hpc_system")
        if sim.get("hpc_system") not in (None, ""):
            names.append("simulation.hpc_system")
    if missing:
        raise SystemExit(f"ERROR: the engage options JSON has no value for {', '.join(missing)}: "
                         f"add {'it' if len(missing) == 1 else 'them'} to the JSON.")
    for k, default in (("job_name", sod["job_name"].get("default")), ("use_segments", "Y"),
                       ("gamera_grid_type", "Q" if component == "grx" else None)):
        if sim.get(k) in (None, "") and default is not None:
            sim[k] = default
            names.append(f"simulation.{k}")
    err = makeitso.grid_type_error(sim["gamera_grid_type"], component)
    if err:
        raise SystemExit(f"ERROR: the engage options JSON: {err}")
    if sim.get("node_type") in (None, "") and sim["hpc_system"] in makeitso.machine_names():
        node_type_od = makeitso.node_type_description(option_makeitso_descriptions, sim["hpc_system"])
        if node_type_od is not None and node_type_od.get("default") not in (None, ""):
            sim["node_type"] = node_type_od["default"]
            names.append("simulation.node_type")
    defaults = {"account_name": makeitso.account_default(sim["hpc_system"], CONFIG_DIR)
                if sim["hpc_system"] in makeitso.machine_names() else None,
                "run_directory": option_makeitso_descriptions["pbs"]["_common"]["run_directory"].get("default"),
                "kaiju_install_directory": KAIJUHOME,
                "kaiju_build_directory": makeitso.kaiju_build_default(component)}
    for k, v in defaults.items():
        if pbs.get(k) in (None, "") and v is not None:
            pbs[k] = v
            names.append(f"pbs.{k}")
        pbs.setdefault(k, v)
    for k, desc in option_engage_descriptions["coupling"].items():
        if k in cpl or (component == "grx" and k == "gcm_spin_up_time"):
            continue                            # set by apply_grx_defaults
        cpl[k] = os.environ.get("CONDA_DEFAULT_ENV") if k == "conda_env" else desc.get("default")
        names.append(f"coupling.{k}")
    if filled is not None:
        filled.extend(names)
    else:
        note_engage_filled(names)
    return engage_options


def note_engage_filled(names):
    if names:
        makeitso._note(makeitso.fill_note_text("engage", names, "engage_parameters.json"))


def _blank(value):
    return value is None or str(value).strip() in ("", "None", "null")


def validate_engage_options(engage_options, option_makeitso_descriptions, option_engage_descriptions,
                            component="grt"):
    """Check and normalize every engage option against its description; exit listing all problems."""
    problems = []
    pbs = engage_options.get("pbs") or {}
    run_dir = os.path.abspath(os.path.expanduser(str(pbs.get("run_directory") or ".")))

    def walk(descs, values, prefix):
        if not isinstance(values, dict):
            return
        for k, desc in descs.items():
            if not isinstance(desc, dict) or "LEVEL" not in desc or k not in values:
                continue
            path = f"{prefix}.{k}"
            if path == "simulation.hpc_system":
                continue                    # checked later against engage's machines
            try:
                # bcwind_file may already be in the run directory; makeitso checks it there
                values[k] = makeitso.normalize_option(path, desc, values[k],
                                                      None if path == "pbs.run_directory" else run_dir,
                                                      check_files=path != "simulation.bcwind_file")
            except ValueError as e:
                problems.append(str(e))

    walk(option_makeitso_descriptions["simulation"], engage_options.get("simulation"), "simulation")
    hpc = str((engage_options.get("simulation") or {}).get("hpc_system"))
    pbs_desc = dict(option_makeitso_descriptions["pbs"]["_common"])
    if hpc in makeitso.machine_names():
        pbs_desc.update(makeitso.get_pbs_descriptions(option_makeitso_descriptions, hpc))
        pbs_desc.update(option_engage_descriptions["pbs"].get(hpc) or {})
    walk(pbs_desc, pbs, "pbs")
    cpl = engage_options.setdefault("coupling", {})
    walk(option_engage_descriptions["coupling"], cpl, "coupling")
    grid = (engage_options.get("simulation") or {}).get("gamera_grid_type")
    if component == "grx":
        if _blank(cpl.get("conda_env")):
            cpl["conda_env"] = os.environ.get("CONDA_DEFAULT_ENV")
    elif component == "grt" and str(grid) != "D" and _blank(cpl.get("conda_env")):
        problems.append(f"coupling.conda_env is required: the grid-{grid} TIE-GCM hand-off regrids in it")
    if problems:
        raise SystemExit("ERROR: " + "\nERROR: ".join(problems))
    return engage_options


def check_grt_spinup_window(engage_options):
    """Exit if the TIE-GCM spin-up window (from the midnight before start - lead-in, to start) is empty."""
    sim, cpl = engage_options["simulation"], engage_options["coupling"]
    start = makeitso.parse_run_datetime(sim["start_date"], "simulation.start_date")
    try:
        lead = float(cpl["gr_warm_up_time"]) + float(cpl["gcm_spin_up_time"])
    except (TypeError, ValueError):
        raise SystemExit(f"ERROR: coupling.gr_warm_up_time={cpl.get('gr_warm_up_time')!r} and "
                         f"coupling.gcm_spin_up_time={cpl.get('gcm_spin_up_time')!r} must be numbers "
                         f"of seconds.") from None
    spin_start = start - datetime.timedelta(seconds=lead)
    spin_start = datetime.datetime.combine(spin_start.date(), datetime.time())
    if spin_start >= start:
        raise SystemExit(f"ERROR: the TIE-GCM spin-up window is empty: coupling.gr_warm_up_time + "
                         f"coupling.gcm_spin_up_time = {lead:g} s before a midnight coupled start "
                         f"({sim['start_date']}). Set coupling.gcm_spin_up_time > 0.")


def _asked(description, mode):
    """True if makeitso.get_run_option prompts for this option in `mode`."""
    level = description.get("LEVEL", "BASIC")
    if description.get("prompt") is None:
        return False
    return level not in makeitso.MODES or makeitso.MODES.index(level) <= makeitso.MODES.index(mode)


def _ask_checked(name, description, mode, check):
    """makeitso.get_run_option, re-asked while check(value) returns an error; exit if it was not asked."""
    while True:
        value = makeitso.get_run_option(name, description, mode)
        err = check(value)
        if not err:
            return value
        if not _asked(description, mode):
            raise SystemExit(err)
        print(f"{makeitso.RED}{err}{makeitso.RESET}")


def prompt_user_for_run_options(args):
    """Prompt for engage's options and return the nested options dict."""
    mode = args.mode
    grx_options.announce(True)

    option_makeitso_descriptions = makeitso.load_option_descriptions()
    option_engage_descriptions = grx_options.load_engage_descriptions(OPTION_ENGAGE_DESCRIPTIONS_FILE)

    options = {}

    o = options["simulation"] = {}
    od = option_makeitso_descriptions["simulation"]
    set_grt_machine_choices(od["hpc_system"], option_engage_descriptions)
    set_grid_choices(od["gamera_grid_type"], getattr(args, "component", "grt"))
    grx = getattr(args, "component", "grt") == "grx"
    gd = grx_options.grx_defaults() if grx else {}
    if grx:
        od["start_date"]["default"] = gd["simulation.start_date"]
        od["stop_date"]["default"] = gd["simulation.stop_date"]
        od["hpc_system"]["default"] = gd["simulation.hpc_system"]
        od["hpc_system"]["valids"] = grx_options.grx_machines()

    for on in ["job_name"]:
        o[on] = makeitso.get_run_option(on, od[on], mode)


    # The run dates are also the date range of the CDAWeb boundary-condition data.
    for on in ["start_date", "stop_date"]:
        if grx and on == "start_date":
            o[on] = _ask_checked(on, od[on], mode, grx_options.coupling_start_error)
        elif not grx and on == "stop_date":
            o[on] = _ask_checked(on, od[on], mode,
                                 lambda v: grx_options.grt_window_error(o["start_date"], v))
        else:
            o[on] = makeitso.get_run_option(on, od[on], mode)

    t1 = makeitso.parse_run_datetime(o["start_date"], "simulation.start_date")
    t2 = makeitso.parse_run_datetime(o["stop_date"], "simulation.stop_date")
    o["start_date"] = start_date = t1.strftime(makeitso.RUN_DATETIME_FORMAT)
    o["stop_date"] = stop_date = t2.strftime(makeitso.RUN_DATETIME_FORMAT)
    if t2 <= t1:
        raise ValueError(
            f"stop_date ({stop_date}) must be strictly after start_date ({start_date}); "
            f"the simulation duration would be non-positive.")
    simulation_duration = (t2 - t1).total_seconds()

    # The segment is asked after the walltime, which sets its default.
    for on in ["gamera_grid_type", "hpc_system"]:
        o[on] = makeitso.get_run_option(on, od[on], mode)
    node_type_od = makeitso.node_type_description(option_makeitso_descriptions, o["hpc_system"])
    if node_type_od is not None:
        o["node_type"] = makeitso.get_run_option("node_type", node_type_od, mode)

    o = options["pbs"] = {}

    od = option_makeitso_descriptions["pbs"]["_common"]
    hpc_for_account = options["simulation"]["hpc_system"]
    makeitso.seed_account_prompt(od["account_name"], hpc_for_account, CONFIG_DIR)
    if not account_required(hpc_for_account):
        od["account_name"]["prompt"] = None
        od["account_name"]["default"] = None
    od["kaiju_install_directory"]["default"] = os.environ["KAIJUHOME"]
    od["kaiju_build_directory"]["default"] = makeitso.kaiju_build_default("grx" if grx else "grt")
    for on in od:
        if on == "account_name" and od[on].get("prompt") is not None:
            check = makeitso.account_prompt_check(hpc_for_account)
            o[on] = _ask_checked(on, od[on], mode, lambda v: (lambda e: e and f"ERROR: {e}")(check(v)))
        elif on == "num_segments":
            o[on] = None                  # derived from the segment below
        else:
            # relative kaiju directories are anchored to the run directory asked before them
            run_dir = (os.path.abspath(os.path.expanduser(str(o["run_directory"])))
                       if o.get("run_directory") not in (None, "") else None)
            o[on] = makeitso.get_run_option(on, od[on], mode, path=f"pbs.{on}", run_dir=run_dir)

    hpc_platform = options["simulation"]["hpc_system"]
    gamera_grid_type = options["simulation"]["gamera_grid_type"]
    od = makeitso.get_pbs_descriptions(option_makeitso_descriptions, hpc_platform)
    if hpc_platform not in option_engage_descriptions["pbs"]:
        raise ValueError(
            f"Coupled GRT runs support {', '.join(option_engage_descriptions['pbs'])}; hpc_system="
            f"{hpc_platform!r} has no engage PBS defaults — pick a supported machine.")
    oed = option_engage_descriptions["pbs"][hpc_platform]
    od["num_helpers"]["default"] = (
        od["num_helpers"]["default"][gamera_grid_type]
    )
    if "local_modules" in od and "local_modules" in oed:
        od["local_modules"] = oed["local_modules"]
    if grx:
        od["walltime"]["default"] = grx_options.grx_defaults(hpc_platform)["pbs.walltime"]
    node_type = options["simulation"].get("node_type")
    component = "grx" if grx else "grt"
    for on in od:
        makeitso.set_thread_default(od, on, o, hpc_platform, node_type)
        if on == "walltime":
            # must fit the queue and hold at least a 1 h coupled segment
            o[on] = _ask_checked(on, od[on], mode, lambda v: (lambda e: e and f"ERROR: {e}")(
                grx_options.walltime_error(None, v, hpc_platform, o.get("queue"), component,
                                           gamera_grid_type)))
        else:
            o[on] = makeitso.get_run_option(on, od[on], mode)

    # Segment length (seconds, whole hours): default is the largest that fits the walltime.
    o = options["simulation"]
    sod = option_makeitso_descriptions["simulation"]
    wall, queue = options["pbs"]["walltime"], options["pbs"].get("queue")
    cod = option_engage_descriptions["coupling"]["gr_warm_up_time"]
    # GRX: only an unasked (default) warm-up constrains the segment choices
    warm0 = cod.get("default") if grx and not _asked(cod, mode) else None
    valid, info, top = grx_options.segment_choices(component, hpc_platform, queue, wall,
                                                   gamera_grid_type, start_date, stop_date, warm0)
    for on in ["use_segments"]:
        sod[on]["default"] = "Y"
        o[on] = "Y" if grx else makeitso.get_run_option(on, sod[on], mode)
    if makeitso.yn(o["use_segments"]) == "Y":
        if not valid:
            raise SystemExit("ERROR: " + grx_options.segment_hours_error(
                component, hpc_platform, queue, wall, gamera_grid_type, start_date, stop_date,
                3600 * max(top, 1), warm0))
        print(grx_options.segment_suggestion(component, hpc_platform, queue, wall, gamera_grid_type,
                                             valid, info, top))
        hod = dict(sod["segment_duration"])
        hod["default"] = grx_options.suggested_seconds(valid, info, top)

        def check_segment(v):
            try:
                secs = float(str(v).strip())
            except (TypeError, ValueError):
                return (f"ERROR: segment_duration {v!r} is not a number of seconds. Valid values "
                        f"(seconds): {grx_options.seconds_list(valid)}")
            err = grx_options.segment_hours_error(component, hpc_platform, queue, wall, gamera_grid_type,
                                                  start_date, stop_date, secs, warm0)
            if not err and grx and not _asked(cod, mode):
                err = grx_options.warmup_error(warm0, secs)
            return err and f"ERROR: {err}"
        secs = _ask_checked("segment_duration", hod, mode, check_segment)
        o["segment_duration"] = str(int(round(float(secs))))
    else:
        o["segment_duration"] = str(simulation_duration)
        err = grx_options.walltime_error(simulation_duration, wall, hpc_platform, queue, component,
                                         gamera_grid_type)
        if err:
            raise ValueError(f"{err} (use_segments=N runs the whole {simulation_duration / 3600:g} h "
                             f"window as one segment: answer Y and pick a segment length)")
    num_segments = math.ceil(simulation_duration / float(o["segment_duration"]))
    options["pbs"]["num_segments"] = str(num_segments)

    options["coupling"] = {}
    o = options["coupling"]
    od = option_engage_descriptions["coupling"]

    od["conda_env"]["default"] = os.environ.get('CONDA_DEFAULT_ENV')
    if grx:
        od["gcm_spin_up_time"]["default"] = gd["coupling.gcm_spin_up_time"]

    sim = options["simulation"]
    for on in ["gr_warm_up_time", "gcm_spin_up_time", "conda_env", "tfin_delta"]:
        if grx and on == "conda_env":
            o[on] = os.environ.get("CONDA_DEFAULT_ENV")     # GRX has no hand-off regrid
        elif grx and on == "gr_warm_up_time":
            o[on] = _ask_checked(on, od[on], mode, lambda v: grx_options.warmup_error(
                v, sim["segment_duration"], sim["start_date"], sim["stop_date"], options["pbs"]["walltime"],
                hpc_platform, options["pbs"].get("queue")))
        else:
            o[on] = makeitso.get_run_option(on, od[on], mode)
    return options


def run_standalone(component, args):
    """Run standalone tiegcmrun (tiegcm) or makeitso (gr) with engage's arguments."""
    tui = getattr(args, "tui", False)
    if component == "tiegcm":
        tg_args = ["--mode", args.mode]
        if getattr(args, "tiegcm_options_path", None):
            tg_args += ["--options_path", args.tiegcm_options_path]
            groups = forwarded_rederive(args, _tiegcm_rederive_groups("tiegcm"))
            if groups:
                tg_args += ["--rederive", groups]
        if args.clobber:
            tg_args.append("--clobber")
        if not tui:
            tg_args.append("--no-tui")
        if args.debug:
            tg_args.append("--debug")
        print("engage: launching standalone tiegcmrun ...")
        try:
            return _import_tiegcmrun().tiegcmrun(tg_args)
        except TOOL_ERRORS as e:
            if args.debug:
                raise
            raise SystemExit(tool_error("tiegcmrun", getattr(args, "tiegcm_options_path", None), e)) from None
    if component == "gr":
        mk_args = {"mode": args.mode, "tui": tui, "clobber": args.clobber,
                   "debug": args.debug, "verbose": args.verbose}
        if getattr(args, "makeitso_options_path", None):
            mk_args["options_path"] = args.makeitso_options_path
            mk_args["rederive"] = forwarded_rederive(args, makeitso.REDERIVE_GROUPS)
        print("engage: launching standalone makeitso (GR: GAMERA+RAIJU) ...")
        try:
            return makeitso.makeitso(mk_args)
        except TOOL_ERRORS as e:
            if args.debug:
                raise
            raise SystemExit(tool_error("makeitso", getattr(args, "makeitso_options_path", None), e)) from None
    raise ValueError(f"run_standalone: unknown component {component!r}")


def _save_waxrun_parameters(coupled_spec, run_dir=None):
    """Copy waxrun's replay payload to <run_dir>/waxrun_parameters.json with relocated paths; return its path."""
    src = coupled_spec.get("waxrun_parameters")
    if not src or not os.path.isfile(src):
        raise RuntimeError(f"waxrun did not report its replay payload (coupled_spec['waxrun_parameters']={src!r})")
    waxrun = _import_waxrun()
    run_dir = os.path.abspath(run_dir or os.getcwd())
    target = os.path.join(run_dir, "waxrun_parameters.json")
    with open(src, "r", encoding="utf-8") as f:
        payload = waxrun.relocate_payload(json.load(f), os.path.dirname(os.path.abspath(src)), run_dir, run_dir)
    with open(target, "w", encoding="utf-8") as f:      # waxrun's own layout
        f.write(json.dumps(payload, indent=2, default=str) + "\n")
    return target


def _recreate_hint(exc):
    """Point a waxrun --clobber refusal at engage's --recreate-cases."""
    return (f"{exc}\n[engage] use engage --recreate-cases to recreate the WACCM-X cases (engage --clobber "
            f"only overwrites engage's replay files).")


def grx_preconditions(engage_options, voltron=None):
    """Exit on any GRX option that run_grx cannot use; return (coupling_start, stop, delay_h, step_h)."""
    sim, pbs, cpl = engage_options["simulation"], engage_options["pbs"], engage_options["coupling"]
    err = grx_options.machine_error(sim.get("hpc_system"))
    if err:
        raise SystemExit(err)
    coupling_start = makeitso.parse_run_datetime(sim["start_date"], "simulation.start_date")
    stop = makeitso.parse_run_datetime(sim["stop_date"], "simulation.stop_date")
    vcpl = (voltron or {}).get("coupling") or {}
    dtc = vcpl.get("dtCouple")
    dtc = None if dtc in (None, "") else dtc
    # Refuse before makeitso's expensive preprocessing.
    if float(cpl["gcm_spin_up_time"]) <= 0:
        raise SystemExit(f"ERROR: coupling.gcm_spin_up_time={cpl['gcm_spin_up_time']!r} must be > 0 (the WACCM-X "
                         f"spin-up target)")
    for err in (grx_options.hours_error("gr_warm_up_time", cpl["gr_warm_up_time"]),
                grx_options.segment_error(sim["start_date"], sim["stop_date"], sim["segment_duration"],
                                          cpl["gr_warm_up_time"], pbs["walltime"], sim["hpc_system"],
                                          pbs.get("queue")),
                grx_options.warmup_error(cpl["gr_warm_up_time"], sim["segment_duration"],
                                         sim["start_date"], sim["stop_date"], pbs["walltime"],
                                         sim["hpc_system"], pbs.get("queue")),
                grx_options.coupling_start_error(sim["start_date"])):
        if err:
            raise SystemExit(err)
    if dtc is not None:
        tspin = ((voltron or {}).get("spinup") or {}).get("tSpin")
        err = grx_options.cadence_error(dtc, tSpin=None if tspin in (None, "") else tspin,
                                        gr_warm_up_time=cpl["gr_warm_up_time"],
                                        segment_duration=sim["segment_duration"],
                                        start_date=sim["start_date"])
        if err:
            raise SystemExit(err)
    err = coupled_queue_error(engage_options, "grx",
                              grt_min_nodes(engage_options, _raw_makeitso_descriptions()))
    if err:
        raise SystemExit(err)
    delay_h = int(float(cpl["gr_warm_up_time"])) // 3600        # gcm_spin_up_time is not added
    step_h = int(float(sim["segment_duration"])) // 3600
    return coupling_start, stop, delay_h, step_h


def grx_srcroot(waxrun, payload, mode, run_dir):
    """Linear GRX: check waxrun's simulation.srcroot before makeitso writes; ask until valid if allowed.

    Returns the asked value (waxrun then does not ask again), or None."""
    import tui_form_waxrun
    opt = waxrun.getparams(mode, engage=payload)
    sim = opt["simulation"]
    check = lambda v: tui_form_waxrun.srcroot_error(v, run_dir)     # noqa: E731
    err = check(sim.get("srcroot"))
    if not err:
        return None
    spec = waxrun.load_descriptions()["simulation"]["srcroot"]
    if "simulation.srcroot" in set(opt.get("_engage_filled", ())) \
            or not waxrun._at_or_below(spec["LEVEL"], mode):
        name = "" if err.startswith("simulation.srcroot") else "simulation.srcroot "
        waxrun._refuse(f"{name}{err}", "simulation.srcroot", opt)
    print(f"{waxrun.RED}{err}{waxrun.RESET}")
    waxrun._prompt_leaf("srcroot", spec, sim, "simulation.srcroot", check=check, run_dir=run_dir)
    return os.path.abspath(os.path.expanduser(str(sim["srcroot"])))


def check_waxrun_options(waxrun, engage_options, makeitso_options, lay, coupling_start, stop, delay_h,
                         step_h, args):
    """Run waxrun.plan_check silently on the options makeitso is about to write (GRX -wo before_write hook)."""
    argv = build_waxrun_argv(engage_options, makeitso_options, lay, coupling_start, stop, delay_h,
                             step_h, args)
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            waxrun.plan_check(argv)
    except FileExistsError as e:
        raise SystemExit(f"ERROR: {_recreate_hint(e)}") from None
    except SystemExit as e:
        if isinstance(e.code, str) and "--clobber" in e.code:
            raise SystemExit(_recreate_hint(e.code)) from None
        raise
    except TOOL_ERRORS as e:
        if args.debug:
            raise
        raise SystemExit(tool_error("waxrun", args.waxrun_options_path, e)) from None
    finally:
        # so the real waxrun pass prints its NOTEs
        waxrun._NOTED.clear()
        waxrun._CONSTITUENT_WARNED.clear()


def run_grx(engage_options, args):
    """Coupled GRX run: makeitso, then waxrun (WACCM-X cases), then engage's jobs and submit script.

    Returns (pbs_files, dag)."""
    waxrun = check_api_version("waxrun", _import_waxrun())
    apply_grx_defaults(engage_options)
    sim, pbs, cpl = engage_options["simulation"], engage_options["pbs"], engage_options["coupling"]
    job_name = sim["job_name"]
    loaded = None
    if getattr(args, "makeitso_options_path", None):
        loaded = load_makeitso_json(args.makeitso_options_path)
    coupling_start, stop, delay_h, step_h = grx_preconditions(engage_options,
                                                              (loaded or {}).get("voltron"))
    t0 = coupling_start - datetime.timedelta(hours=delay_h)
    t1 = stop
    run_dir = os.path.abspath(os.path.expanduser(str(pbs["run_directory"])))
    _refuse_existing_artifacts(run_dir, GRX_ARTIFACTS, args.clobber)

    assert_grx_stack(sim["hpc_system"], _waxrun_runtime_modules(waxrun, sim["hpc_system"]),
                     f"waxrun options_description stacks.{sim['hpc_system']}.runtime_modules")

    lay = _grx_layout(waxrun, job_name, run_dir)
    phase = getattr(args, "phase", None) or "both"
    replay = getattr(args, "waxrun_options_path", None)

    # A voltron without 'ranks of W' cannot couple WACCM-X and leaves CAM hanging.
    voltron_path = _voltron_exe(pbs)
    volt_ident = voltron_preflight(voltron_path, "grx")

    # waxrun's wizard runs before anything is written; mage.dtcouple/dtout come from makeitso later.
    answers = None
    if getattr(args, "tui", False) and not replay:
        payload = build_waxrun_payload(engage_options, {"voltron": {}}, lay, coupling_start, stop,
                                       delay_h, step_h, args)
        answers = waxrun.wizard_answers(payload, args.mode, hidden=("mage.dtcouple", "mage.dtout"))
    # Linear prompts: waxrun asks after makeitso, so settle the CAM tree now.
    srcroot = None
    if answers is None and not replay:
        payload = build_waxrun_payload(engage_options, {"voltron": {}}, lay, coupling_start, stop,
                                       delay_h, step_h, args)
        srcroot = grx_srcroot(waxrun, payload, args.mode, lay["rundir"])

    if loaded is not None:
        patch_makeitso_options(loaded, engage_options, os.path.basename(args.makeitso_options_path))
        err = makeitso_rules_error(loaded, "grx")
        if err:
            raise SystemExit(err)

    # --skip-build replay: refuse a missing cesm.exe before the run directory is created.
    if replay and getattr(args, "skip_build", False) and not getattr(args, "dry_run", False) \
            and phase != "spinup":
        check_waxrun_options(waxrun, engage_options, loaded or {"voltron": {}}, lay, coupling_start,
                             stop, delay_h, step_h, args)

    os.makedirs(lay["rundir"], exist_ok=True)

    # Replay JSONs are written only after makeitso returns, so a refused run needs no --clobber to retry.
    makeitso_args = {'clobber': args.clobber, 'debug': args.debug, 'verbose': args.verbose,
                     'mode': args.mode, 'tui': getattr(args, "tui", False), 'save_parameters': False,
                     'grx': {'run_dir': lay['rundir'], 't0': t0.strftime(GRX_DATE_FORMAT),
                             't1': t1.strftime(GRX_DATE_FORMAT),
                             'stephours': step_h, 'delaycouplinghrs': delay_h}}
    cwd = os.getcwd()
    tmp_mo = None
    if loaded is not None:
        # absolute run_directory: makeitso reads a relative one from the JSON's directory
        tmp_mo = write_temp_json(lay["rundir"], "makeitso_parameters.json",
                                 dict(loaded, pbs=dict(loaded.get("pbs") or {}, run_directory=lay["rundir"])))
        makeitso_args['options_path'] = tmp_mo
        makeitso_args['rederive'] = forwarded_rederive(args, makeitso.REDERIVE_GROUPS)
    makeitso_args.update(engage_options)
    if replay:
        makeitso_args["before_write"] = lambda options: check_waxrun_options(
            waxrun, engage_options, options, lay, coupling_start, stop, delay_h, step_h, args)
    try:
        makeitso_options, volt_files, _ = makeitso.makeitso(makeitso_args)
    except TOOL_ERRORS as e:
        if args.debug:
            raise
        raise SystemExit(tool_error("makeitso", getattr(args, "makeitso_options_path", None), e)) from None
    finally:
        os.chdir(cwd)   # makeitso's GRX mode cds into the run dir
        if tmp_mo and os.path.exists(tmp_mo):
            os.remove(tmp_mo)
    rundir = lay["rundir"]
    with open(os.path.join(rundir, "engage_parameters.json"), "w", encoding="utf-8") as f:
        json.dump(makeitso.run_options_for_json(engage_options, rundir, rundir), f, indent=JSON_INDENT)
    with open(os.path.join(rundir, "makeitso_parameters.json"), "w", encoding="utf-8") as f:
        json.dump(makeitso.run_options_for_json(makeitso_options, rundir, rundir), f, indent=JSON_INDENT)

    argv = build_waxrun_argv(engage_options, makeitso_options, lay, coupling_start, stop,
                             delay_h, step_h, args, answers,
                             {"simulation": {"srcroot": srcroot}} if srcroot else None)
    try:
        coupled_spec, spinup_submit_scripts, coupled_case_files = waxrun.waxrun(argv)
    except FileExistsError as e:
        raise FileExistsError(_recreate_hint(e)) from None
    except SystemExit as e:
        if isinstance(e.code, str) and "--clobber" in e.code:
            raise SystemExit(_recreate_hint(e.code)) from None
        raise
    os.chdir(cwd)
    cesm_exe = coupled_spec.get("cesm_exe")
    if phase != "spinup" and getattr(args, "skip_build", False) and not getattr(args, "dry_run", False) \
            and not (cesm_exe and os.path.isfile(cesm_exe)):
        raise SystemExit(f"--skip-build but the coupled cesm.exe does not exist: {cesm_exe} "
                         "(build the coupled case, or drop --skip-build)")
    # coupled_case_files is [caseroot, user_nl_cam], not one per segment.
    if coupled_spec.get("runtime_modules"):
        assert_grx_stack(sim["hpc_system"], coupled_spec["runtime_modules"], "coupled_spec['runtime_modules']")
    exe = coupled_spec.get("voltron_exe") or voltron_path
    if os.path.abspath(exe) != os.path.abspath(voltron_path):
        volt_ident = voltron_preflight(exe, "grx")
    link_voltron(lay["rundir"], exe, warn=False)        # voltron_preflight warned already
    link_lead_in_voltron(lay["rundir"])

    pbs_files = create_grx_pbs_scripts(engage_options, makeitso_options, coupled_spec,
                                       t0, t1, step_h, delay_h, run_dir=lay["rundir"], phase=phase,
                                       voltron=volt_ident)
    dag = write_grx_dag(lay["rundir"], job_name, coupled_spec, spinup_submit_scripts, pbs_files,
                        lay, clobber=args.clobber or getattr(args, "recreate_cases", False),
                        phase=phase)

    # re-dump the replay JSONs with the segment count actually used
    engage_options["pbs"]["num_segments"] = str(len(pbs_files))
    try:
        with open(coupled_spec.get("waxrun_parameters"), "r", encoding="utf-8") as f:
            warn = grx_spinup_warning(cpl["gcm_spin_up_time"], json.load(f).get("spinup"))
    except (OSError, ValueError, TypeError, KeyError):
        warn = None
    if warn:
        print(warn)
    dump_json(os.path.join(rundir, "engage_parameters.json"),
              makeitso.run_options_for_json(engage_options, rundir, rundir))
    dump_json(os.path.join(rundir, "makeitso_parameters.json"),
              makeitso.run_options_for_json(makeitso_options, rundir, rundir))
    _save_waxrun_parameters(coupled_spec, rundir)

    print(f"GRX voltron configs = {volt_files}")
    print(f"GRX spinup submits  = {spinup_submit_scripts}")
    print(f"GRX coupled case    = {coupled_case_files}")
    print(f"GRX segment PBS     = {pbs_files}")
    print(f"GRX submit script   = {dag}")
    return pbs_files, dag


ENGAGE_REDERIVE_GROUPS = {
    "job": "engage's pbs machine settings (queue, walltime, job_priority / group_list / "
           "local_modules, ranks and threads per node), from config/machines.yaml for "
           "simulation.hpc_system; also forwarded to makeitso and tiegcmrun",
}


def _tiegcm_rederive_groups(component):
    """tiegcmrun's --rederive groups, or {} without a TIE-GCM part."""
    if component not in ("grt", "tiegcm"):
        return {}
    _import_tiegcmrun()
    return dict(getattr(importlib.import_module("replay_solver"), "REDERIVE_GROUPS", {}))


def parse_rederive(value, component="grt"):
    """Parse --rederive into a tuple of groups ('all' expands, 'list' stays); ValueError on unknown ones."""
    if value in (None, ""):
        return ()
    known = list(dict.fromkeys(list(ENGAGE_REDERIVE_GROUPS) + list(makeitso.REDERIVE_GROUPS)
                               + list(_tiegcm_rederive_groups(component))))
    groups = [g.strip().lower() for g in str(value).split(",") if g.strip()]
    if groups == ["list"]:
        return ("list",)
    if "all" in groups:
        return tuple(known)
    unknown = [g for g in groups if g not in known]
    if unknown or not groups:
        raise ValueError(f"unknown --rederive group {', '.join(unknown) or repr(value)} "
                         f"(groups: {', '.join(known)}, all, list)")
    return tuple(dict.fromkeys(groups))


def rederive_groups_text(component="grt"):
    """Text for --rederive list."""
    lines = []
    for tool, groups in (("engage", ENGAGE_REDERIVE_GROUPS), ("makeitso", makeitso.REDERIVE_GROUPS),
                         ("tiegcmrun", _tiegcm_rederive_groups(component))):
        lines += [f"{tool:<9} {name:<11} {what}" for name, what in groups.items()]
    return "\n".join(lines)


def forwarded_rederive(args, groups):
    """The requested --rederive groups that a tool knows, as 'g1,g2', or None."""
    mine = [g for g in (getattr(args, "rederive", None) or ()) if g in groups]
    return ",".join(mine) or None


# The option files each --component takes: all of them, or none.
MODE_FILES = {"grt": ("-eo", "-mo", "-to"), "grx": ("-eo", "-mo", "-wo"), "gr": ("-mo",), "tiegcm": ("-to",)}
OPTION_FILE_ARGS = {"-eo": "engage_options_path", "-mo": "makeitso_options_path",
                    "-to": "tiegcm_options_path", "-wo": "waxrun_options_path"}


def _and_list(items):
    return items[0] if len(items) == 1 else f"{', '.join(items[:-1])} and {items[-1]}"


def option_files_error(args):
    """Error if the option files given do not match --component (MODE_FILES), else None."""
    component = getattr(args, "component", "grt")
    wanted = MODE_FILES[component]
    given = [f for f, attr in OPTION_FILE_ARGS.items() if getattr(args, attr, None)]
    for f in given:
        if f not in wanted:
            return f"ERROR: {f} is not used by --component {component} (it takes {', '.join(wanted)})."
    missing = [f for f in wanted if f not in given]
    if given and missing:
        return (f"ERROR: --component {component} takes {_and_list(list(wanted))} together; missing: "
                f"{', '.join(missing)}. Give all of them, or none to be prompted.")
    return None


EO_SECTIONS = ("simulation", "pbs", "coupling", "waxrun")


def engage_key_schema(component):
    """The keys an -eo JSON may hold, as {key: sub-schema or None}."""
    od = makeitso.options_key_schema()
    with open(OPTION_ENGAGE_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        oed = json.load(f)
    pbs = dict(od["pbs"])
    for block in oed["pbs"].values():
        pbs.update(dict.fromkeys(block))
    schema = {"simulation": od["simulation"], "pbs": pbs, "coupling": dict.fromkeys(oed["coupling"])}
    if component == "grx":
        schema["waxrun"] = None
    return makeitso.add_schema_paths(schema, list(EO_PREVIOUS_KEYS) + [".".join(p) for p in OLD_FORMAT_KEYS])


# Keys from older engage_parameters.json files that are still accepted.
EO_PREVIOUS_KEYS = ("coupling.dtOut", "coupling.gamera_spin_up_time")


def option_file_problem(flag, path, shown, component):
    """One-line problem with an option file (unreadable, bad shape, old format, unknown keys), or None."""
    try:
        data = makeitso.read_options_json(path, shown)
    except ValueError as e:
        return str(e)
    problem = None
    if flag == "-eo":
        problem = next((f"{s} must be a JSON object" for s in EO_SECTIONS
                        if s in data and not isinstance(data[s], dict)), None)
        problem = problem or old_format_problem(data) or makeitso.unknown_keys_problem(
            data, engage_key_schema(component))
    elif flag == "-mo":
        problem = (makeitso.section_problem(data) or makeitso.old_format_problem(data)
                   or makeitso.unknown_keys_problem(data, makeitso.options_key_schema()))
    elif flag == "-to":
        check = getattr(_import_tiegcmrun(), "options_file_problem", None)
        problem = check(data, unknown=True) if check else None
    elif flag == "-wo":
        check = getattr(_import_waxrun(), "options_shape_problem", None)
        problem = check(data) if check else None
    return f"{shown}: {problem}" if problem else None


def check_option_files(args, shown):
    """Exit on the first option file with a problem."""
    component = getattr(args, "component", "grt")
    for flag, attr in OPTION_FILE_ARGS.items():
        path = getattr(args, attr, None)
        if path:
            problem = option_file_problem(flag, path, shown.get(attr, path), component)
            if problem:
                raise SystemExit(f"ERROR: {problem}")


def main():
    """engage's command-line entry point."""
    parser = create_command_line_parser()

    args = parser.parse_args()
    args.tui = makeitso.use_tui(args.no_tui)
    if args.debug:
        print(f"args = {args}")
    _REPORTED.clear()
    # Absolute now: the run later works inside the run directory.
    shown = {}      # as the user spelled them, for error messages
    for _opt in ("engage_options_path", "makeitso_options_path", "tiegcm_options_path",
                 "waxrun_options_path"):
        if getattr(args, _opt, None):
            shown[_opt] = getattr(args, _opt)
            setattr(args, _opt, os.path.abspath(getattr(args, _opt)))
    err = option_files_error(args)
    if err:
        raise SystemExit(err)
    clobber = args.clobber
    debug = args.debug
    engage_options_path = args.engage_options_path

    component = getattr(args, "component", "grt")
    if debug or args.verbose:
        print(f"makeitso from {makeitso.__file__}")
    makeitso.require_env("KAIJUHOME", KAIJU_ENV_HINT)
    if component == "grx":
        _import_waxrun_module("ic_solver").require_ncks()
    if component in ("grt", "tiegcm"):
        check_api_version("tiegcmrun", _import_tiegcmrun())
        if debug or args.verbose:
            print(f"tiegcmrun from {_import_tiegcmrun().__file__}")
    check_option_files(args, shown)
    try:
        args.rederive = parse_rederive(getattr(args, "rederive", None), component)
    except ValueError as e:
        raise SystemExit(f"ERROR: {e}") from None
    if args.rederive == ("list",):
        print(rederive_groups_text(component))
        return None
    if component in ("tiegcm", "gr"):
        return run_standalone(component, args)

    # GRT replay files gated by --clobber (GRX: GRX_ARTIFACTS in run_grx).
    _artifacts = ["engage_parameters.json", "makeitso_parameters.json", "tiegcmrun_parameters.json"]

    if engage_options_path:
        option_makeitso_descriptions = makeitso.load_option_descriptions()
        option_engage_descriptions = grx_options.load_engage_descriptions(OPTION_ENGAGE_DESCRIPTIONS_FILE)
        with open(engage_options_path, "r", encoding="utf-8") as f:
            engage_options = json.load(f)
        problem = old_format_problem(engage_options)
        if problem:
            raise SystemExit(f"ERROR: {engage_options_path}: {problem}")
        eo_run_directory_note(engage_options, engage_options_path)
        # the NOTE on filled entries is printed once the segment is known (below)
        filled = []
        fill_engage_options(engage_options, option_makeitso_descriptions, option_engage_descriptions,
                            component, filled)
        if component == "grx":
            filled += [f"{sec}.{k}" for sec, k in (("simulation", "segment_duration"),
                                                    ("coupling", "gcm_spin_up_time"))
                       if _blank(engage_options[sec].get(k))]
            apply_grx_defaults(engage_options)
        validate_engage_options(engage_options, option_makeitso_descriptions, option_engage_descriptions,
                                component)
        o = engage_options["simulation"]
        t1 = makeitso.parse_run_datetime(o["start_date"], "simulation.start_date")
        t2 = makeitso.parse_run_datetime(o["stop_date"], "simulation.stop_date")
        o["start_date"] = start_date = t1.strftime(makeitso.RUN_DATETIME_FORMAT)
        o["stop_date"] = stop_date = t2.strftime(makeitso.RUN_DATETIME_FORMAT)
        if t2 <= t1:
            raise ValueError(
                f"stop_date ({stop_date}) must be strictly after start_date ({start_date}); "
                f"the simulation duration would be non-positive.")
        simulation_duration = float((t2 - t1).total_seconds())

        makeitso.check_node_type(engage_options)
        o = engage_options["pbs"]

        hpc_platform = engage_options["simulation"]["hpc_system"]
        gamera_grid_type = engage_options["simulation"]["gamera_grid_type"]
        od = makeitso.get_pbs_descriptions(option_makeitso_descriptions, hpc_platform)
        if hpc_platform not in option_engage_descriptions["pbs"]:
            raise ValueError(
                f"Coupled GRT runs support {', '.join(option_engage_descriptions['pbs'])}; hpc_system="
                f"{hpc_platform!r} has no engage PBS defaults — pick a supported machine.")
        oed = option_engage_descriptions["pbs"][hpc_platform]
        od["num_helpers"]["default"] = (
            od["num_helpers"]["default"][gamera_grid_type]
        )
        if "local_modules" in od and "local_modules" in oed:
            od["local_modules"] = oed["local_modules"]
        node_type = engage_options["simulation"].get("node_type")
        rederived = set()
        if "job" in args.rederive:
            rederived = {f for f in makeitso.MACHINE_PBS_FIELDS if f in o}
            for f in makeitso.MACHINE_PBS_FIELDS:
                o.pop(f, None)
        else:
            engage_machine_hint(engage_options, od)
        for on in od:
            if on not in o:
                makeitso.set_thread_default(od, on, o, hpc_platform, node_type)
                o[on] = od[on]["default"]
                if on not in rederived:
                    filled.append(f"pbs.{on}")
        settle_engage_account(engage_options)

        sim = engage_options["simulation"]
        if sim.get("segment_duration") in (None, ""):
            if sim["use_segments"] != "Y":
                sim["segment_duration"] = str(simulation_duration)
            else:
                valid, info, top = grx_options.segment_choices(
                    component, hpc_platform, o.get("queue"), o.get("walltime"), gamera_grid_type,
                    start_date, stop_date, None)
                if not valid:
                    raise SystemExit("ERROR: " + grx_options.segment_hours_error(
                        component, hpc_platform, o.get("queue"), o.get("walltime"), gamera_grid_type,
                        start_date, stop_date, 3600 * max(top, 1), None))
                sim["segment_duration"] = str(int(round(float(grx_options.suggested_seconds(valid, info, top)))))
            filled.append("simulation.segment_duration")
        note_engage_filled(filled)
        num_segments = 1
        if sim["use_segments"] == "Y":
            num_segments = math.ceil(simulation_duration / float(sim["segment_duration"]))
        o["num_segments"] = str(num_segments)

        o = engage_options["coupling"]
        if "conda_env" not in engage_options["coupling"]:
            o["conda_env"] = os.environ.get('CONDA_DEFAULT_ENV')
    else:
        # The wizard falls back to the linear prompts without a TTY or prompt_toolkit.
        engage_options = None
        if getattr(args, "tui", False):
            import tui_form_engage
            try:
                engage_options = tui_form_engage.run_engage_wizard(args)
            except tui_form_engage._FormCancelled:
                print("engage wizard cancelled (Ctrl-C); nothing was written.")
                raise SystemExit(130) from None
            except tui_form_engage._FormUnavailable as e:
                print(f"Engage wizard unavailable ({e}); using the linear prompts.")
                if debug:
                    import traceback
                    traceback.print_exc()
        if engage_options is None:
            engage_options = prompt_user_for_run_options(args)
        settle_engage_account(engage_options, "the engage options")
        validate_engage_options(engage_options, makeitso.load_option_descriptions(),
                                grx_options.load_engage_descriptions(OPTION_ENGAGE_DESCRIPTIONS_FILE),
                                component)
    if debug:
        print(f"engage_options = {engage_options}")

    if component == "grx":
        return run_grx(engage_options, args)

    check_grt_spinup_window(engage_options)
    err = coupled_queue_error(engage_options, "grt",
                              grt_min_nodes(engage_options, _raw_makeitso_descriptions()))
    if err:
        raise SystemExit(err)
    launch_dir = os.getcwd()
    run_dir = makeitso.prepare_run_directory(engage_options["pbs"]["run_directory"])
    _refuse_existing_artifacts(run_dir, _artifacts, clobber)
    try:
        return _run_grt(engage_options, args, run_dir, launch_dir)
    finally:
        os.chdir(launch_dir)


def eo_run_directory_note(engage_options, eo_path):
    """Print and return a NOTE if a relative pbs.run_directory resolves differently from the cwd than from the JSON."""
    rd = (engage_options.get("pbs") or {}).get("run_directory")
    if rd in (None, "") or os.path.isabs(os.path.expanduser(str(rd))):
        return None
    jdir = os.path.dirname(os.path.abspath(eo_path))
    here = os.path.abspath(os.path.expanduser(str(rd)))
    there = os.path.normpath(os.path.join(jdir, os.path.expanduser(str(rd))))
    if here == there:
        return None
    note = (f"NOTE: pbs.run_directory {rd!r} is relative to the current directory, not to the directory "
            f"of {os.path.basename(eo_path)}.")
    print(note)
    return note


def _refuse_existing_artifacts(directory, artifacts, clobber):
    """Raise FileExistsError if a replay file exists in `directory` and clobber is off."""
    for artifact in artifacts:
        path = os.path.normpath(os.path.join(directory, artifact))
        if os.path.exists(path) and not clobber:
            raise FileExistsError(f"Options file {path} exists! (pass --clobber to overwrite)")


def _api_at_least(module, version):
    """True if the tool's API_VERSION is at least (major, minor)."""
    try:
        have = tuple(int(x) for x in str(getattr(module, "API_VERSION", "0")).split(".")[:2])
    except ValueError:
        return False
    return have >= tuple(version)


def grt_tiegcm_args(engage_options, voltron, run_dir, args):
    """Return tiegcmrun's argv for a coupled GRT run (engage's options plus makeitso's voltron section)."""
    coupled_options = copy.deepcopy(engage_options)
    coupled_options["voltron"] = voltron
    coupled_options["pbs"]["run_directory"] = run_dir
    node_type = makeitso.node_type_name(engage_options["simulation"]["hpc_system"],
                                        engage_options["simulation"].get("node_type"))
    if node_type:
        coupled_options["simulation"]["node_type"] = node_type
    tiegcm_options_path = args.tiegcm_options_path
    tiegcm_args = [
        "--coupling",
        "--engage", json.dumps(coupled_options),
        "--options_path", tiegcm_options_path,
        "--mode", args.mode,
    ]
    if tiegcm_options_path:
        groups = forwarded_rederive(args, _tiegcm_rederive_groups("grt"))
        if groups:
            tiegcm_args += ["--rederive", groups]
    if args.debug:
        tiegcm_args.append("--debug")
    if not getattr(args, "tui", False) and not tiegcm_options_path:
        tiegcm_args.append("--no-tui")
    return tiegcm_args


def check_tiegcm_options(engage_options, makeitso_options, run_dir, args):
    """Run tiegcmrun --check silently on the options makeitso is about to write (GRT -to before_write hook)."""
    argv = grt_tiegcm_args(engage_options, makeitso_options.get("voltron", {}), run_dir, args) + ["--check"]
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            _import_tiegcmrun().tiegcmrun(argv)
    except TOOL_ERRORS as e:
        if args.debug:
            raise
        raise SystemExit(tool_error("tiegcmrun", args.tiegcm_options_path, e)) from None
    finally:
        os.chdir(cwd)


def _run_grt(engage_options, args, run_dir, launch_dir):
    """Coupled GRT run: makeitso (from launch_dir), then tiegcmrun and the coupled jobs in run_dir.

    tiegcmrun runs with run_dir as cwd because TIE-GCM's coupled files land in its cwd."""
    clobber, debug, verbose = args.clobber, args.debug, args.verbose
    makeitso_options_path = args.makeitso_options_path
    tiegcm_options_path = args.tiegcm_options_path
    warn = grt_machine_warning(engage_options["simulation"]["hpc_system"])
    if warn:
        print(warn)

    # A prompted run names the TIE-GCM executables inside tiegcmrun; they are checked after it.
    grt_to_format_preflight(tiegcm_options_path)
    grt_to_resource_preflight(tiegcm_options_path, engage_options)
    to_exes = _tiegcm_exes_of(tiegcm_options_path)
    if to_exes is not None:
        grt_executable_preflight(to_exes["modelexe"], to_exes["coupled_modelexe"], engage_options,
                                 zitop=to_exes["zitop"])
    else:
        voltron_preflight(_voltron_exe(engage_options["pbs"]), "grt")

    # Replay JSONs are written only after makeitso returns, so a refused run needs no --clobber to retry.
    makeitso_args = {'clobber': clobber, 'debug': debug, 'verbose': verbose, 'mode': args.mode,
                     'tui': getattr(args, "tui", False), 'save_parameters': False}
    tmp_mo = None
    if makeitso_options_path:
        with open(OPTION_MAKEITSO_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
            option_makeitso_descriptions = json.load(f)
        makeitso_options = load_makeitso_json(makeitso_options_path)
        loaded_makeitso_options = copy.deepcopy(makeitso_options)
        assign_engage_options(makeitso_options, engage_options)
        # The GR run starts gr_warm_up_time before engage's coupled start.
        t0 = makeitso.parse_run_datetime(engage_options["simulation"]["start_date"], "simulation.start_date")
        t0 -= datetime.timedelta(seconds=float(engage_options["coupling"]["gr_warm_up_time"]))
        t1 = makeitso.parse_run_datetime(engage_options["simulation"]["stop_date"], "simulation.stop_date")
        makeitso_options["simulation"]["start_date"] = datetime.datetime.isoformat(t0)
        segment_duration = float(engage_options["simulation"]["segment_duration"])
        # written as a fresh run writes them, so an unchanged replay rewrites the same JSON
        makeitso_options["voltron"]["time"]["tFin"] = (t1-t0).total_seconds()
        makeitso_options["pbs"]["num_segments"] = str(math.ceil((t1-t0).total_seconds()/segment_duration))
        makeitso_options["pbs"] = {k: makeitso_options["pbs"][k] for k in engage_options["pbs"]
                                   if k in makeitso_options["pbs"]} | makeitso_options["pbs"]
        makeitso.apply_grid_layout(makeitso_options, option_makeitso_descriptions)
        note_engage_owned(loaded_makeitso_options, makeitso_options, engage_options,
                          os.path.basename(makeitso_options_path))
        err = makeitso_rules_error(makeitso_options, "grt", engage_options["coupling"])
        if err:
            raise SystemExit(err)
        # absolute run_directory: makeitso reads a relative one from the JSON's directory
        tmp_mo = write_temp_json(run_dir, "makeitso_parameters.json",
                                 dict(makeitso_options, pbs=dict(makeitso_options["pbs"], run_directory=run_dir)))
        makeitso_args['options_path'] = tmp_mo
        makeitso_args['rederive'] = forwarded_rederive(args, makeitso.REDERIVE_GROUPS)
    makeitso_args.update(engage_options)
    if tiegcm_options_path and _api_at_least(_import_tiegcmrun(), (2, 1)):
        makeitso_args["before_write"] = lambda options: check_tiegcm_options(
            engage_options, options, run_dir, args)

    # makeitso resolves pbs.run_directory against the launch directory and moves into it.
    os.chdir(launch_dir)
    try:
        makeitso_options, makeitso_spinup_pbs_scripts, makeitso_warmup_pbs_scripts = makeitso.makeitso(makeitso_args)
    except TOOL_ERRORS as e:
        if debug:
            raise
        raise SystemExit(tool_error("makeitso", makeitso_options_path, e)) from None
    finally:
        if tmp_mo and os.path.exists(tmp_mo):
            os.remove(tmp_mo)
    makeitso_pbs_scripts = makeitso_spinup_pbs_scripts + makeitso_warmup_pbs_scripts
    os.chdir(run_dir)
    # engage records makeitso's num_segments (over the warm-up-extended span)
    dump_json(os.path.join(run_dir, "makeitso_parameters.json"),
              makeitso.run_options_for_json(makeitso_options, run_dir, run_dir))
    engage_options["pbs"]["num_segments"] = makeitso_options["pbs"]["num_segments"]
    dump_json(os.path.join(run_dir, "engage_parameters.json"),
              makeitso.run_options_for_json(engage_options, run_dir, run_dir))

    tiegcm_args = grt_tiegcm_args(engage_options, makeitso_options.get("voltron", {}), run_dir, args)

    try:
        tiegcm_options, tiegcm_pbs_scripts, tiegcm_inp_scripts = _import_tiegcmrun().tiegcmrun(tiegcm_args)
    except TOOL_ERRORS as e:
        if debug:
            raise
        raise SystemExit(tool_error("tiegcmrun", tiegcm_options_path, e)) from None
    except SystemExit as e:
        if e.code == 130:
            print(f"tiegcmrun wizard cancelled (Ctrl-C). Already written in {run_dir}: "
                  f"makeitso_parameters.json, engage_parameters.json and makeitso's GR decks, job "
                  f"scripts and inputs ({', '.join(os.path.relpath(p, run_dir) for p in makeitso_pbs_scripts)}"
                  f"{', ...' if makeitso_pbs_scripts else ''}). Nothing was removed; rerun with --clobber "
                  f"to regenerate the run there.")
        raise

    # The replay JSON holds the spin-up options; tiegcm_options are the coupled ones.
    tiegcm_replay_options = getattr(tiegcm_options, "replay_options", None) or tiegcm_options
    dump_json(os.path.join(run_dir, 'tiegcmrun_parameters.json'),
              tiegcm_json_relative(tiegcm_replay_options, run_dir))

    # On refusal, remove makeitso's GR-only job scripts that the coupled scripts would have replaced.
    try:
        err = deck_count_error(makeitso_options, tiegcm_inp_scripts)
        if err:
            raise SystemExit(err)
        identities = grt_executable_preflight(
            tiegcm_replay_options["model"]["data"].get("modelexe"),
            tiegcm_options["model"]["data"].get("coupled_modelexe"), engage_options,
            zitop=tiegcm_options["model"].get("specification", {}).get("zitop"))
    except SystemExit as e:
        removed = remove_gr_only_scripts(run_dir, engage_options["simulation"]["job_name"])
        if removed:
            raise SystemExit(f"{e.code}\n(no coupled job script was written; makeitso's GR-only "
                             f"{', '.join(removed)} {'was' if len(removed) == 1 else 'were'} "
                             f"removed)") from None
        raise
    link_voltron(run_dir, os.path.join(makeitso_options["pbs"]["kaiju_build_directory"], "bin",
                                       "voltron_mpi.x"), warn=False)
    link_lead_in_voltron(run_dir)

    pbs_scripts, submit_all_jobs_script = create_pbs_scripts(
        engage_options, makeitso_options, makeitso_pbs_scripts, tiegcm_options, tiegcm_inp_scripts,
        tiegcm_pbs_scripts, run_dir=run_dir, voltron=identities["voltron"])
    print(f"GR_pbs_scripts = {makeitso_pbs_scripts}")
    print(f"Tiegcm_pbs_scripts = {tiegcm_pbs_scripts}")
    print(f"GRT_pbs_scripts = {pbs_scripts}")
    print(f"submit_all_jobs_script = {submit_all_jobs_script}")


def cli():
    """Run main(), printing ValueError / FileExistsError as one ERROR line unless --debug."""
    try:
        return main()
    except (ValueError, FileExistsError) as e:
        if "--debug" in sys.argv[1:] or "-d" in sys.argv[1:]:
            raise
        raise SystemExit(f"ERROR: {ANSI_RE.sub('', str(e)).strip()}") from None


if __name__ == "__main__":
    cli()