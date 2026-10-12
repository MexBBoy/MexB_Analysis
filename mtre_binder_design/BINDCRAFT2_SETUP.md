# BindCraft2 — step-by-step install

From opening a terminal to a first design run. Verified against
[PacesaLab/BindCraft2](https://github.com/PacesaLab/BindCraft2) `c0ee638`
(10 October 2026).

Allow about 45 minutes, nearly all of it the 5.3 GB weights download.

---

## Before you start

| | |
| --- | --- |
| **Operating system** | Linux. On Windows use WSL2 with an NVIDIA driver installed on the Windows side. macOS is not supported, Apple Silicon included. |
| **GPU** | An NVIDIA card. **There is no CPU install** — `install.sh cpu` is refused outright, not quietly built. |
| **Python** | 3.12 or newer. If you have none, the installer brings its own. |
| **Disk** | ~20 GB free. |

---

## Step 1 — Open a terminal

On the Linux machine with the GPU. If you use a cluster, open a terminal on the
machine you will actually run jobs from.

## Step 2 — Check the GPU

```bash
nvidia-smi
```

You want a table listing your card. Note two things from it:

- **"CUDA Version"** in the top-right corner — the newest CUDA your driver serves.
- **Your card's name**, which tells you its generation.

If this says `command not found`, there is no NVIDIA driver on this machine and
nothing below will work. Stop here.

## Step 3 — Check git and disk space

```bash
git --version
df -h ~
```

If git is missing: `sudo apt install git` on Ubuntu/Debian. You need ~20 GB
free in the "Avail" column.

## Step 4 — Deactivate any Python environment

**This step is easy to skip and causes confusing failures.** The installer
installs into whatever environment is active. If a conda env or venv is active
and it is older than Python 3.12, the install stops; if it is a different
project's environment, BindCraft2 lands inside it.

```bash
conda deactivate 2>/dev/null
deactivate 2>/dev/null
echo "VIRTUAL_ENV=$VIRTUAL_ENV  CONDA_PREFIX=$CONDA_PREFIX"
```

Both should print empty. Run `conda deactivate` more than once if you have
stacked environments — some setups auto-activate `base` on login.

*(If you deliberately want BindCraft2 inside an environment you already use,
activate that one instead and skip to step 6. It must be Python 3.12+.)*

## Step 5 — Choose a permanent location

**The install is editable, which means the clone IS the installation.** Moving,
renaming or deleting the folder afterwards breaks the `bindcraft` command. Do
not clone into `/tmp` or any scratch space that gets purged.

```bash
mkdir -p ~/software
cd ~/software
```

## Step 6 — Clone the repository

```bash
git clone https://github.com/PacesaLab/BindCraft2.git
cd BindCraft2
```

Note the org is **PacesaLab**, not `martinpacesa`. The old repository is
BindCraft v1 and stops at v1.5.3.

## Step 7 — Run the installer

For most machines, no argument — it reads the accelerator off the driver:

```bash
bash install.sh
```

**Pass an argument instead if any of these apply:**

| Your situation | Command |
| --- | --- |
| V100, P100, or any card older than a T4 | `bash install.sh cuda12` |
| AMD card | `bash install.sh rocm` |
| Cluster login node, compute nodes have older cards | `bash install.sh cuda12` |
| You want to skip the 5.3 GB download for now | `bash install.sh --no-weights` |

Why: CUDA 13 dropped every GPU below compute capability 7.5, so a V100 or P100
needs the CUDA 12 wheels however new its driver is. The installer detects this
from the card — but on a login node with no visible GPU it cannot, and assumes
CUDA 13. **Match the choice to the compute nodes, not the login node.**

This takes 20–40 minutes. It builds a `.venv`, installs the wheels, downloads
the AlphaFold parameters to `~/.cache/bindcraft/alphafold`, and self-checks.

A successful run ends with something like:

```
bindcraft: jax runs here on gpu, 1 device(s)
bindcraft: installed, designing on cuda13.
  source .venv/bin/activate            # once per terminal
```

## Step 8 — Activate the environment

```bash
source .venv/bin/activate
```

Your prompt should now show `(.venv)`.

## Step 9 — Verify

Three checks. Run all three.

```bash
bindcraft --help
```
Prints usage. If `command not found`, the environment is not active — go back
to step 8.

```bash
python -c "import jax; print(jax.devices())"
```
Want: `[CudaDevice(id=0)]`, one per card. **If this says `[CpuDevice(id=0)]`,
stop.** The install works but a campaign would run on the CPU and take days
instead of hours. Re-run step 7 naming your accelerator explicitly.

```bash
python -m bindcraft.selfcheck cuda13     # or cuda12 / rocm, whichever you installed
```
Names every module and checkpoint it needs and reports anything missing. A pip
install can resolve a wheel and unpack nothing while reporting success; this is
what catches that.

## Step 10 — First design run

A short bounded run, to confirm the whole thing works before committing hours:

```bash
bindcraft design examples/pdl1.json \
  --set number_of_final_designs=1 \
  --set max_trajectories=3 \
  --set 'project_folder=results/smoke_test'
```

This designs against the shipped human PD-L1 target and stops after three
attempts whether or not anything passes. It should take a few minutes and write
to `results/smoke_test/`.

Before it starts you will see a `campaign preflight:` line naming the models,
target and output folder — read it and check it matches what you intended. If
something is wrong you get `campaign refused:` listing **every** problem at
once, rather than one at a time.

Once that works, the real thing:

```bash
bindcraft design examples/pdl1.json
```

which runs until it has 10 accepted binders, with no cap on attempts.

To read the results afterwards:

```bash
bindcraft rank results/pdl1              # best designs first
bindcraft filter results/pdl1 --list     # what you can filter on
```

## Step 11 — Every new terminal

```bash
cd ~/software/BindCraft2
source .venv/bin/activate
```

The activation does not persist between terminals. This is the single most
common "it stopped working" after a successful install.

---

## Running against your own target

Put a `design.json` next to your structure:

```json
{
  "targets": [{"name": "my_target",
               "target_path": "target.pdb",
               "chains": "A",
               "hotspots": "54,56,66-70"}],
  "modality": "binder",
  "binder_lengths": [60, 100],
  "number_of_final_designs": 10,
  "project_folder": "results/my_target"
}
```

```bash
bindcraft design design.json
```

- `binder_lengths` `[80, 80]` pins one length; `[60, 100]` samples the range.
- `hotspots` uses residue numbers from your own structure; chain-prefixed
  (`A54,B12-16`) when you select several chains.
- Each experiment gets its own `project_folder`.

Useful listings: `bindcraft design --list-targets`, `--list-modalities`,
`--list-properties`, `--list-settings`.

---

## For the MtrE project

Our three targets translate directly. The crown campaign:

```json
{
  "targets": [{"name": "MtrE_crown",
               "target_path": "targets/target_crown.pdb",
               "chains": "A,B,C",
               "hotspots": "A326,A333,C114,C119"}],
  "modality": "binder",
  "binder_lengths": [55, 100],
  "number_of_final_designs": 100,
  "project_folder": "results/mtre_crown"
}
```

Two BindCraft2 features worth using here that v1 does not have:

```bash
# Cyclic peptides: replaces the separate RFpeptides arm. 7-16mers, head to
# tail, cysteine biased out of the sequence.
bindcraft design mtre_mouth.json --modality cyclic_peptide

# Negative design against TolC: selectivity as a design objective rather than
# a post-hoc decoy filter. Add E. coli TolC under "targets" with
# "objective": "detarget".
```

See `examples/pdl1_cyclic_peptide.json` and `examples/pdl1_detarget_pd1.json`
for the full shape of each.

---

## Troubleshooting

| Symptom | Cause | Fix |
| --- | --- | --- |
| `bindcraft: command not found` | environment not active | `source .venv/bin/activate` |
| `[CpuDevice(id=0)]` | wheels do not match the driver | re-run `bash install.sh` naming the accelerator |
| Campaign runs but is very slow | JAX silently took the CPU | same device check; it warns only once |
| `... is active and BindCraft2 needs 3.12 or newer` | an old environment is active | step 4, then re-run the installer |
| `module not installed: <name>` | wheel resolved, unpacked nothing | `bash install.sh` again in the same environment |
| `checkpoints missing or unfinished` | interrupted download | `rm -rf ~/.cache/bindcraft/alphafold && bindcraft fetch-weights` |
| `AlphaFold parameters are not on this machine` | installed with `--no-weights` | `bindcraft fetch-weights` |
| Works on login node, fails on compute nodes | installer read the login node's driver | `bash install.sh cuda12` |
| Out of GPU memory | too many workers packed onto the card | `export BINDCRAFT_WORKERS_PER_GPU=1` |
| `already holds campaign output` | `resume: false` on a non-empty folder | use a fresh `project_folder` |
| `bindcraft: this environment has no ...` | wrong interpreter | `source .venv/bin/activate` |

**GPU memory per worker**, by padded complex size: 7.1 GB at 64 residues, 8.0
at 128, 11.8 at 256, 18.0 at 384, 26.7 at 512. A campaign packs up to seven
workers per card automatically. Our MtrE targets are 42–180 residues, so one
worker fits a 16 GB card comfortably.

**Reporting a problem upstream:** include `campaign_metadata.json` from the
results folder, the `campaign refused:` or worker log text, and the output of
`python -m bindcraft.selfcheck <accelerator>`.
