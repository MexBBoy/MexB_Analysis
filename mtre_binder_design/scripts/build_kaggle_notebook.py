#!/usr/bin/env python3
"""Generate the Kaggle notebook for running the MtrE campaigns on a free GPU.

    python3 scripts/build_kaggle_notebook.py

Writes notebooks/mtre_kaggle.ipynb. The notebook is kept in generated form so
the cell sources stay reviewable as ordinary Python here rather than as JSON
string arrays.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "notebooks" / "mtre_kaggle.ipynb"

REPO = "https://github.com/MexBBoy/MexB_Analysis.git"
BRANCH = "claude/upbeat-allen-akc6ua"

cells: list[dict] = []


def md(text: str) -> None:
    cells.append({"cell_type": "markdown", "metadata": {},
                  "source": text.strip("\n").splitlines(keepends=True)})


def code(text: str) -> None:
    cells.append({"cell_type": "code", "metadata": {}, "execution_count": None,
                  "outputs": [], "source": text.strip("\n").splitlines(keepends=True)})


# ---------------------------------------------------------------------------
md("""
# MtrE binder design on Kaggle

Runs the MtrE campaigns on Kaggle's free GPU (30 h/week, 16 GB P100 or 2×T4,
12 h sessions).

**Before you run anything**, in the notebook sidebar:

1. **Session options → Accelerator → GPU** (P100 or T4 ×2; either works).
2. **Session options → Internet → On.** Needs a phone-verified Kaggle account.
   Without it every install below fails.
3. *(after your first run)* **Input → attach your `mtre-weights` dataset**, so
   the model weights are not re-downloaded every session. Cell 3 tells you how
   to create it.

**What this notebook will and will not do.** It will prove the pipeline end to
end and produce real designs. It will not run a full campaign: BindCraft needs
several hundred trajectories per target and a session is capped at 12 hours.
Use this to validate and to generate a first handful of designs; scale on a lab
GPU or an HPC allocation.

**This notebook has not been executed.** There is no GPU in the environment it
was written in, so every cell is written to verify itself and fail loudly
rather than silently. If a cell fails, paste the output of cell 1 together with
the error.
""")

# ---------------------------------------------------------------------------
md("""
## 1 · Preflight

Checks the things that actually break: no GPU, no internet, wrong CUDA, no
disk. If anything here is wrong, stop and fix it before going on.
""")

code("""
import os, shutil, subprocess, sys, textwrap

report = {}
report["python"] = sys.version.split()[0]

try:
    import torch
    report["torch"] = torch.__version__
    report["torch_cuda"] = torch.version.cuda
    report["cuda_available"] = torch.cuda.is_available()
    if torch.cuda.is_available():
        report["gpu"] = torch.cuda.get_device_name(0)
        cap = torch.cuda.get_device_capability(0)
        report["compute_capability"] = f"sm_{cap[0]}{cap[1]}"
        report["vram_gb"] = round(
            torch.cuda.get_device_properties(0).total_memory / 1024**3, 1)
except ImportError:
    report["torch"] = "NOT INSTALLED"

free_gb = round(shutil.disk_usage("/kaggle/working").free / 1024**3, 1)
report["working_free_gb"] = free_gb

try:
    subprocess.run(["curl", "-sS", "--max-time", "15", "-o", "/dev/null",
                    "https://github.com"], check=True, capture_output=True)
    report["internet"] = True
except Exception:
    report["internet"] = False

print("=" * 60)
for k, v in report.items():
    print(f"{k:22} {v}")
print("=" * 60)

problems = []
if not report.get("cuda_available"):
    problems.append("No GPU. Sidebar -> Session options -> Accelerator -> GPU.")
if not report["internet"]:
    problems.append("No internet. Sidebar -> Session options -> Internet -> On "
                    "(needs a phone-verified account).")
if free_gb < 12:
    problems.append(f"Only {free_gb} GB free in /kaggle/working; need ~12 GB "
                    "unless weights come from an attached dataset.")

if problems:
    print("\\nSTOP:")
    for p in problems:
        print(" - " + p)
    raise SystemExit("preflight failed")
print("\\npreflight OK")
""")

# ---------------------------------------------------------------------------
md("""
## 2 · Get the project

Clones the MtrE repository (targets, configs, filter scripts) and RFdiffusion.
""")

code(f"""
import os, subprocess, pathlib

WORK = pathlib.Path("/kaggle/working")
PROJ = WORK / "mtre"          # our repo
RFD  = WORK / "RFdiffusion"   # the tool

def run(cmd, **kw):
    print("$", " ".join(cmd) if isinstance(cmd, list) else cmd)
    subprocess.run(cmd, shell=isinstance(cmd, str), check=True, **kw)

if not PROJ.exists():
    run(["git", "clone", "--depth", "1", "-b", "{BRANCH}", "{REPO}", str(PROJ)])
if not RFD.exists():
    run(["git", "clone", "--depth", "1",
         "https://github.com/RosettaCommons/RFdiffusion", str(RFD)])

DESIGN = PROJ / "mtre_binder_design"
assert (DESIGN / "targets" / "target_loop2.pdb").exists(), "project layout unexpected"
print("\\nproject:", DESIGN)
print("targets:", sorted(p.name for p in (DESIGN / "targets").glob("target_*.pdb")))
""")

# ---------------------------------------------------------------------------
md("""
## 3 · Model weights

**Only one RFdiffusion checkpoint is needed.** All three of our campaigns set
`ppi.hotspot_res`, and RFdiffusion selects `Complex_base_ckpt.pt` for that case
(`rfdiffusion/inference/model_runners.py`). The other six checkpoints are for
motif scaffolding, inpainting and symmetry, none of which we use — so this
downloads one file instead of 5 GB.

### Make this permanent (do it once)

Re-downloading every session wastes your weekly quota. After this cell has run:

1. Kaggle sidebar → **Output** → find `/kaggle/working/weights`
2. **Create Dataset** from it, name it `mtre-weights`, keep it **Private**
3. In later sessions, attach it via **Input → Add Input** and skip the download

The cell checks `/kaggle/input/**` first and only downloads what is missing.
""")

code("""
import pathlib, subprocess, os

WORK = pathlib.Path("/kaggle/working")
WEIGHTS = WORK / "weights"
WEIGHTS.mkdir(exist_ok=True)

# Only what hotspot-driven binder design actually loads.
NEEDED = {
    "Complex_base_ckpt.pt":
        "http://files.ipd.uw.edu/pub/RFdiffusion/"
        "e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt",
}
# Uncomment if you later run the scaffold-guided variant in run_lumen_plug.sh:
# NEEDED["Complex_Fold_base_ckpt.pt"] = (
#     "http://files.ipd.uw.edu/pub/RFdiffusion/"
#     "60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt")

def find_in_inputs(name):
    root = pathlib.Path("/kaggle/input")
    if not root.exists():
        return None
    for p in root.rglob(name):
        return p
    return None

for name, url in NEEDED.items():
    dest = WEIGHTS / name
    if dest.exists():
        print(f"already here: {name}")
        continue
    found = find_in_inputs(name)
    if found:
        dest.symlink_to(found)
        print(f"from attached dataset: {name}  <- {found}")
        continue
    print(f"downloading {name} (this is the slow bit) ...")
    # NOTE: files.ipd.uw.edu is plain HTTP, not HTTPS. That is upstream's doing.
    subprocess.run(["wget", "-q", "--show-progress", "-O", str(dest), url], check=True)

for name in NEEDED:
    p = WEIGHTS / name
    size = p.stat().st_size if p.exists() else 0
    assert size > 1_000_000, f"{name} missing or truncated ({size} bytes)"
    print(f"{name}: {size/1024**2:.0f} MB")
print("\\nweights ready at", WEIGHTS)
""")

# ---------------------------------------------------------------------------
md("""
## 4 · Install RFdiffusion

Kaggle already ships PyTorch built against CUDA, so we install into the
existing environment rather than building a conda env — much faster, and it
sidesteps the stale `env/SE3nv.yml` pins (PyTorch 1.9 / CUDA 11.1) that break
on modern cards.

**`dgl` is the fragile dependency.** It must match the installed PyTorch and
CUDA exactly. The cell tries the matching wheel index and falls back; if it
fails, the error tells you which combination was attempted.

Kaggle's P100 is `sm_60` and the T4 is `sm_75`, both old enough that the
CUDA-version problem that bites RTX 4090 and H100 users does not apply here.
""")

code("""
import subprocess, sys, torch, re

def pip(*args):
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", *args], check=True)

torch_ver = torch.__version__.split("+")[0]
cuda_ver = (torch.version.cuda or "").replace(".", "")   # "12.1" -> "121"
print(f"torch {torch_ver}, cuda {torch.version.cuda}")

pip("e3nn==0.5.1", "hydra-core", "pyrsistent", "opt_einsum", "icecream",
    "pynvml", "dllogger @ git+https://github.com/NVIDIA/dllogger")

try:
    import dgl  # noqa: F401
    print("dgl already present:", dgl.__version__)
except ImportError:
    mm = ".".join(torch_ver.split(".")[:2])
    indexes = [
        f"https://data.dgl.ai/wheels/torch-{mm}/cu{cuda_ver}/repo.html",
        f"https://data.dgl.ai/wheels/cu{cuda_ver}/repo.html",
        "https://data.dgl.ai/wheels/repo.html",
    ]
    for idx in indexes:
        print("trying", idx)
        try:
            pip("dgl", "-f", idx)
            import dgl  # noqa: F401
            print("dgl installed:", dgl.__version__)
            break
        except Exception as exc:
            print("  failed:", str(exc)[:200])
    else:
        raise SystemExit(
            f"Could not install dgl for torch {torch_ver} / cu{cuda_ver}.\\n"
            "Paste this cell's output back and I'll pin a working combination."
        )
""")

code("""
import os, subprocess, sys, pathlib
RFD = pathlib.Path("/kaggle/working/RFdiffusion")

# NVIDIA's SE(3)-Transformer, vendored in the RFdiffusion repo.
subprocess.run([sys.executable, "setup.py", "install"],
               cwd=RFD / "env" / "SE3Transformer", check=True,
               stdout=subprocess.DEVNULL)
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-e", "."],
               cwd=RFD, check=True)

# Verify the import chain the inference script actually needs.
check = subprocess.run(
    [sys.executable, "-c",
     "import torch, dgl, e3nn, hydra, se3_transformer, rfdiffusion; "
     "print('imports OK; cuda:', torch.cuda.is_available())"],
    capture_output=True, text=True)
print(check.stdout or check.stderr)
assert check.returncode == 0, "RFdiffusion import chain is broken, see above"
""")

# ---------------------------------------------------------------------------
md("""
## 5 · Smoke test

Two designs against `target_loop2.pdb` (42 residues — the smallest target).
This should take a couple of minutes and proves the install, the weights, the
contig string and the target file all work together.

If this passes, everything downstream is a matter of scale.
""")

code("""
import subprocess, sys, pathlib

PROJ = pathlib.Path("/kaggle/working/mtre/mtre_binder_design")
RFD = pathlib.Path("/kaggle/working/RFdiffusion")
OUT = pathlib.Path("/kaggle/working/designs/smoke")
OUT.parent.mkdir(exist_ok=True)

cmd = [
    sys.executable, str(RFD / "scripts" / "run_inference.py"),
    f"inference.model_directory_path=/kaggle/working/weights",
    f"inference.output_prefix={OUT}/test",
    f"inference.input_pdb={PROJ / 'targets' / 'target_loop2.pdb'}",
    "inference.num_designs=2",
    "contigmap.contigs=[A315-335/0 B315-335/0 55-80]",
    "ppi.hotspot_res=[A326,A322,A327]",
    "denoiser.noise_scale_ca=0",
    "denoiser.noise_scale_frame=0",
]
print(" ".join(cmd), "\\n")
subprocess.run(cmd, check=True)

made = sorted(OUT.glob("*.pdb"))
print("\\ndesigns produced:", [p.name for p in made])
assert made, "no PDBs written - check the traceback above"
""")

code("""
# Geometry filter: rejects designs docked into the membrane or the periplasm,
# before any scoring compute is spent on them.
import subprocess, sys, pathlib
PROJ = pathlib.Path("/kaggle/working/mtre/mtre_binder_design")
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "biopython"], check=True)
subprocess.run(
    [sys.executable, str(PROJ / "scripts" / "filter_designs.py"),
     "--mode", "loop2"] + [str(p) for p in
                           sorted(pathlib.Path("/kaggle/working/designs/smoke").glob("*.pdb"))],
    check=True, cwd=PROJ)
""")

# ---------------------------------------------------------------------------
md("""
## 6 · A real RFdiffusion run

Scale up once the smoke test is green. Pick a campaign and a design count that
fits the time you have left in the session — roughly 30–60 s per design on a
P100 for targets this size, so 200 designs is about 2–3 hours.

Outputs land in `/kaggle/working/designs/`, which Kaggle saves as notebook
output when the session ends. Save before you run out of session.
""")

code("""
CAMPAIGN = "crown"      # "crown" | "lumen" | "macrocycle"
NUM_DESIGNS = 50        # raise once you know the per-design time

import subprocess, sys, pathlib, json
PROJ = pathlib.Path("/kaggle/working/mtre/mtre_binder_design")
RFD = pathlib.Path("/kaggle/working/RFdiffusion")
contigs = json.loads((PROJ / "targets" / "mtre_sites.json").read_text())["contigs"]

SPECS = {
    "crown": dict(
        pdb="target_crown.pdb",
        contigs=f"[{contigs['crown']} 55-95]",
        hotspots="[A326,A333,C114,C119]", extra=[]),
    "lumen": dict(
        pdb="target_lumen.pdb",
        contigs=f"[{contigs['lumen']} 60-100]",
        hotspots="[A123,A312,A311,A132]", extra=[]),
    # Macrocycle: peptide contig FIRST, because cyc_chains='a' cyclises the
    # first output chain. Reversing this silently cyclises part of the target.
    "macrocycle": dict(
        pdb="target_lumen.pdb",
        contigs=f"[12-12 {contigs['lumen']}]",
        hotspots="[A111,A120,A123]",
        extra=["--config-name", "base", "inference.cyclic=True",
               "inference.cyc_chains=a", "diffuser.T=50"]),
}
spec = SPECS[CAMPAIGN]
OUT = pathlib.Path(f"/kaggle/working/designs/{CAMPAIGN}")
OUT.parent.mkdir(parents=True, exist_ok=True)

cmd = [sys.executable, str(RFD / "scripts" / "run_inference.py")] + spec["extra"] + [
    "inference.model_directory_path=/kaggle/working/weights",
    f"inference.output_prefix={OUT}/{CAMPAIGN}",
    f"inference.input_pdb={PROJ / 'targets' / spec['pdb']}",
    f"inference.num_designs={NUM_DESIGNS}",
    f"contigmap.contigs={spec['contigs']}",
    f"ppi.hotspot_res={spec['hotspots']}",
    "denoiser.noise_scale_ca=0", "denoiser.noise_scale_frame=0",
]
print(" ".join(cmd), "\\n")
subprocess.run(cmd, check=True)
print("\\ndone:", len(list(OUT.glob("*.pdb"))), "designs")
""")

# ---------------------------------------------------------------------------
md("""
## 7 · BindCraft (optional, heavy)

BindCraft is the primary generator — 50% experimental hit rate against
RFdiffusion's 24% in the benchmark this project follows — but it is the
expensive half of the setup: PyRosetta plus 5.3 GB of AlphaFold2 weights, and
it wants several hundred trajectories before designs start passing filters.

On a free 12-hour session that is enough to confirm it runs and to collect a
few accepted designs, not to finish a campaign. Set `RUN_BINDCRAFT = True` only
once the RFdiffusion half is working and you have hours of session left.

The install follows BindCraft's own Colab notebook (pip, not conda).
""")

code("""
RUN_BINDCRAFT = False   # flip to True when you have the session time

if RUN_BINDCRAFT:
    import subprocess, sys, pathlib, os
    WORK = pathlib.Path("/kaggle/working")
    BC = WORK / "BindCraft"
    if not BC.exists():
        subprocess.run(["git", "clone", "--depth", "1",
                        "https://github.com/martinpacesa/BindCraft", str(BC)], check=True)
    os.chmod(BC / "functions" / "dssp", 0o755)
    os.chmod(BC / "functions" / "DAlphaBall.gcc", 0o755)

    def pip(*a):
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", *a], check=True)

    # JAX must be the CUDA build or BindCraft silently runs on CPU.
    pip("-U", "jax[cuda12]")
    pip("git+https://github.com/sokrypton/ColabDesign.git", "--no-deps")
    pip("pyrosetta", "--find-links",
        "https://west.rosettacommons.org/pyrosetta/quarterly/release.cxx11thread.serialization")

    import jax
    print("jax devices:", jax.devices())
    assert any("cuda" in str(d).lower() or "gpu" in str(d).lower() for d in jax.devices()), \\
        "JAX is on CPU - BindCraft would run but take days. Fix before continuing."

    # AlphaFold2 weights, 5.3 GB. Reuse from an attached dataset if present.
    params = BC / "params"
    params.mkdir(exist_ok=True)
    found = next(pathlib.Path("/kaggle/input").rglob("params_model_5_ptm.npz"), None) \\
        if pathlib.Path("/kaggle/input").exists() else None
    if found:
        for f in found.parent.iterdir():
            tgt = params / f.name
            if not tgt.exists():
                tgt.symlink_to(f)
        print("AF2 params from attached dataset:", found.parent)
    elif not (params / "params_model_5_ptm.npz").exists():
        print("downloading AF2 params (5.3 GB) ...")
        subprocess.run("apt-get install -y -qq aria2", shell=True, check=False)
        subprocess.run(
            "aria2c -q -x 16 -d /kaggle/working "
            "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar",
            shell=True, check=True)
        subprocess.run(f"tar -xf /kaggle/working/alphafold_params_2022-12-06.tar -C {params}",
                       shell=True, check=True)
        os.remove("/kaggle/working/alphafold_params_2022-12-06.tar")
    assert (params / "params_model_5_ptm.npz").exists(), "AF2 weights missing"
    print("BindCraft ready")
else:
    print("skipped - set RUN_BINDCRAFT = True to install")
""")

code("""
if RUN_BINDCRAFT:
    import subprocess, sys, pathlib, json, tempfile
    PROJ = pathlib.Path("/kaggle/working/mtre/mtre_binder_design")
    BC = pathlib.Path("/kaggle/working/BindCraft")

    # The repo's config uses relative paths; rewrite them absolute for Kaggle
    # and lower the target count so the session ends with something saved.
    cfg = json.loads((PROJ / "configs" / "bindcraft" / "target_crown_loops.json").read_text())
    cfg["starting_pdb"] = str(PROJ / "targets" / "target_crown.pdb")
    cfg["design_path"] = "/kaggle/working/designs/bindcraft_crown/"
    cfg["number_of_final_designs"] = 5
    tmp = pathlib.Path("/kaggle/working/bc_target.json")
    tmp.write_text(json.dumps({k: v for k, v in cfg.items() if not k.startswith("_")}, indent=2))

    subprocess.run([
        sys.executable, "-u", str(BC / "bindcraft.py"),
        "--settings", str(tmp),
        "--filters", str(PROJ / "configs" / "bindcraft" / "filters_mtre_strict.json"),
        "--advanced", str(BC / "settings_advanced" / "default_4stage_multimer.json"),
    ], check=True)
""")

# ---------------------------------------------------------------------------
md("""
## 8 · Before the session ends

Kaggle saves `/kaggle/working` as notebook output when the session finishes,
but only up to 20 GB and only if the session ends cleanly. Anything large or
important is better pulled out explicitly.
""")

code("""
import pathlib, shutil, subprocess
D = pathlib.Path("/kaggle/working/designs")
if D.exists():
    n = len(list(D.rglob("*.pdb")))
    print(f"{n} design PDBs in {D}")
    for sub in sorted(p for p in D.iterdir() if p.is_dir()):
        print(f"  {sub.name}: {len(list(sub.glob('*.pdb')))}")
    # One archive is easier to download than hundreds of files.
    shutil.make_archive("/kaggle/working/mtre_designs", "zip", D)
    size = pathlib.Path("/kaggle/working/mtre_designs.zip").stat().st_size
    print(f"\\nwrote /kaggle/working/mtre_designs.zip ({size/1024**2:.0f} MB)")
    print("Download it from the Output tab, or commit the notebook to keep it.")
else:
    print("no designs yet")

# Keep the weights out of the saved output unless you are building the dataset.
print("\\nIf you are NOT creating the weights dataset, delete them before saving:")
print("  shutil.rmtree('/kaggle/working/weights')")
""")

# ---------------------------------------------------------------------------
md("""
## Troubleshooting

| Symptom | Cause | Fix |
| --- | --- | --- |
| Cell 1 says no internet | Internet toggle off | Sidebar → Session options → Internet → On (phone-verified account required) |
| `dgl` install fails | no wheel for this torch/CUDA pair | paste cell 4's output back; the pair needs pinning |
| `CUDA out of memory` | target too large for 16 GB | use `target_loop2.pdb`, or reduce the binder length range |
| JAX reports CPU only | CUDA JAX build missing | re-run the BindCraft cell; `jax[cuda12]` must install cleanly |
| Weights download very slow | `files.ipd.uw.edu` is plain HTTP and far from Kaggle | create the `mtre-weights` dataset once and attach it thereafter |
| Session died mid-run | 12 h cap, or idle timeout | lower `NUM_DESIGNS`; Kaggle does not resume interrupted runs |
| BindCraft accepts nothing | normal | several hundred trajectories is expected; free tier will not finish a campaign |

## Where this is going

Free compute is for **validating the pipeline and producing first designs**.
The campaign numbers in `PROTOCOL_MtrE.md` — 3,000 RFdiffusion backbones, 48–96
BindCraft designs ordered after filtering — need a lab GPU or an HPC
allocation. Getting designs out of this notebook is what makes that a
straightforward ask rather than a speculative one.
""")

notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python"},
        "accelerator": "GPU",
        "kaggle": {"accelerator": "nvidiaTeslaT4", "dataSources": [],
                   "isInternetEnabled": True, "language": "python",
                   "sourceType": "notebook"},
    },
    "nbformat": 4,
    "nbformat_minor": 4,
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(notebook, indent=1))
print(f"wrote {OUT} ({len(cells)} cells)")
