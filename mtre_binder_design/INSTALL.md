# Installing RFdiffusion and BindCraft locally

Verified against the repositories as of 2 October 2026:
RosettaCommons/RFdiffusion `86507b6`, martinpacesa/BindCraft `7713aa0`.

Budget about 2 hours, most of it downloads.

> **No GPU yet?** `notebooks/mtre_kaggle.ipynb` runs the whole pipeline on
> Kaggle's free tier (30 GPU hours a week, 16 GB card, 12-hour sessions). It is
> enough to validate the setup and produce first designs, not to finish a
> campaign. Start there if you are still arranging compute.

---

## 0. Before you start

**You need a Linux machine with an NVIDIA GPU.** There is no way around this:

| | |
| --- | --- |
| **Linux x86_64 + NVIDIA GPU** | Supported. Do this. |
| **Windows** | Only via WSL2 with an NVIDIA driver on the Windows side. Not native. |
| **macOS** | Not supported, Apple Silicon included. No CUDA, and neither pipeline has a working CPU path. |

**GPU memory.** BindCraft's authors recommend ≥32 GB for large targets. Our
trimmed MtrE targets are 42–180 residues, which runs comfortably on 16–24 GB
(RTX 4090, A5000, A100). RFdiffusion needs roughly 8–12 GB at these sizes.

**Disk.** About 30 GB: AlphaFold2 weights 5.3 GB, RFdiffusion weights ~5 GB,
and the two conda environments together ~20 GB.

**Check your GPU and driver first.** The CUDA version in the top-right of this
output is the maximum your driver supports, and you will need it in step 2:

```bash
nvidia-smi
```

**Install Miniforge** if you have no conda. Prefer it over Anaconda: it
defaults to conda-forge and ships `mamba`, which resolves these environments in
minutes rather than tens of minutes.

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# reopen your shell afterwards
```

---

## 1a. BindCraft2 — the successor, and probably what you want

Verified against PacesaLab/BindCraft2 `c0ee638` (10 Oct 2026).

> **Installing BindCraft2 on its own?** `BINDCRAFT2_SETUP.md` is a clean
> step-by-step walkthrough from opening a terminal, with verification at each
> stage. This section is the condensed version for people also installing v1
> and RFdiffusion.

**BindCraft2 is a separate repository under a new org**
([PacesaLab/BindCraft2](https://github.com/PacesaLab/BindCraft2)), not a branch
or release of `martinpacesa/BindCraft` — that repo stops at v1.5.3. BC2 is a
rewrite with a different CLI, a different input schema, and a far better
installer. Read §1b before deciding; our configs currently target v1.

```bash
git clone https://github.com/PacesaLab/BindCraft2.git ~/software/BindCraft2
cd ~/software/BindCraft2
bash install.sh
source .venv/bin/activate
```

That is the whole install. No conda, no `--cuda` guess: `install.sh` reads the
accelerator off the driver, builds a `.venv` with whichever of `uv`, `venv` or
conda the machine has (fetching `uv` if there is no usable Python at all),
installs the matching wheels, downloads the AlphaFold parameters, and then
verifies the result.

**Two things that will catch you out:**

- **The install is editable, so the clone *is* the installation.** Moving,
  renaming or deleting the directory breaks the `bindcraft` command. Do not
  clone it into scratch space that gets purged.
- **`source .venv/bin/activate` once per terminal.**

Allow ~20 GB free: the 5.3 GB parameter archive and its unpacked contents sit
side by side during the transfer, and the accelerator wheels are several GB
more. ProteinMPNN weights ship inside the package — nothing to fetch.

### Accelerator selection

Usually automatic. Override by passing an argument:

| Argument | When |
| --- | --- |
| *(none)* | read from the driver |
| `cuda13` | compute capability 7.5+ (T4, A100, L40S, 4090, H100) |
| `cuda12` | **Volta, Pascal, Maxwell** — V100, P100, and anything below 7.5 |
| `rocm` | AMD, against a local ROCm 7 install |
| `--no-weights` | skip the AF2 download; `bindcraft fetch-weights` later |

CUDA 13 dropped every GPU below compute capability 7.5, so a V100 or P100 takes
the CUDA 12 wheels however new its driver is — the installer detects this, but
name it yourself if you install on a login node whose card differs from the
compute nodes. Match the choice to the **compute** nodes.

Note `rocm` is offered. That does not rescue RFdiffusion, which remains
CUDA-only, but it means an AMD facility is not automatically out for BC2.

**Verify:**

```bash
python -c "import jax; print(jax.devices()); assert jax.default_backend() == 'gpu'"
python -m bindcraft.selfcheck cuda13     # name your accelerator
bindcraft --help
```

`[CpuDevice(id=0)]` means a campaign would run on the CPU. The installer flags
this itself when the machine has a GPU that JAX did not take.

### Running it

The CLI and input format are both different from v1:

```bash
bindcraft design my_design.json
```

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

Compare that with `configs/bindcraft/target_crown_loops.json` — same target,
same hotspots, different schema. There is no `--filters` / `--advanced` pair;
named presets supply those, overridden per-key with `--set`.

### GPU memory

Per worker, for a padded complex of N residues: 7.1 GB at 64, 8.0 at 128, 11.8
at 256, 18.0 at 384, 26.7 at 512. A campaign packs up to seven workers per card
automatically. Our targets are 42–180 residues, so one worker fits a 16 GB card
comfortably; cap packing with `BINDCRAFT_WORKERS_PER_GPU=1` if you hit limits.

---

## 1b. BindCraft v1 — what our configs currently target

The primary generator for this project as it stands. Install this one if you
want to run the campaigns exactly as specified in `PROTOCOL_MtrE.md`.

```bash
git clone https://github.com/martinpacesa/BindCraft ~/software/BindCraft
cd ~/software/BindCraft
bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'mamba'
```

Set `--cuda` to a version **at or below** what `nvidia-smi` reported. If you
leave it blank the installer guesses, and it often guesses wrong — a wrong
guess surfaces later as JAX silently falling back to CPU, which looks like
"BindCraft is just very slow" rather than like an error.

The script creates a `BindCraft` conda environment on Python 3.10, installs
JAX with CUDA support, ColabDesign and PyRosetta, downloads the AlphaFold2
weights (5.3 GB) into `params/`, and makes `functions/dssp` and
`functions/DAlphaBall.gcc` executable.

> **PyRosetta licensing.** The installer pulls in PyRosetta, which is free for
> academic use but requires a paid licence for commercial use. Fine for us;
> worth knowing before anyone from industry asks to reuse the pipeline.

**Verify:**

```bash
conda activate BindCraft
python -c "import jax; print(jax.devices())"     # must list a CudaDevice, not CpuDevice
python -c "import colabdesign, pyrosetta; print('ok')"
ls params/params_model_5_ptm.npz
```

If `jax.devices()` shows only `CpuDevice`, the CUDA build did not take. Re-run
the installer with the correct `--cuda`; do not try to patch JAX by hand.

**Run it on our crown target:**

```bash
conda activate BindCraft
cd ~/software/BindCraft
python -u ./bindcraft.py \
  --settings /path/to/mtre_binder_design/configs/bindcraft/target_crown_loops.json \
  --filters  /path/to/mtre_binder_design/configs/bindcraft/filters_mtre_strict.json \
  --advanced ./settings_advanced/default_4stage_multimer.json
```

Paths inside the target JSON (`starting_pdb`, `design_path`) are resolved
relative to your working directory, so either run from the directory that makes
them valid or edit them to absolute paths.

Expect to run **several hundred trajectories** before designs start passing
filters, and a few thousand on a hard target. This is normal and is why the GPU
matters.

---

## 2. RFdiffusion

The fold-diversity arm, and the tool for the macrocycle work. Harder to install,
because the pinned dependency stack is from 2021.

### Read this before choosing a route

The shipped `env/SE3nv.yml` pins **PyTorch 1.9 and cudatoolkit 11.1**. That
combination has no compiled kernels for anything newer than Ampere, so on an
**RTX 4090 / L40S (sm_89) or H100 (sm_90)** it fails with
`CUDA error: no kernel image is available` or
`nvrtc: invalid value for --gpu-architecture`. This is a known, long-standing
issue ([RFdiffusion #289](https://github.com/RosettaCommons/RFdiffusion/issues/289)),
and the README says as much: the maintainers ship one CUDA 11.1 file and leave
customisation to the user.

So pick by your hardware:

| Your GPU | Route |
| --- | --- |
| V100, RTX 2080/3090, A100, A5000/A6000 (sm_70–sm_86) | **2a**, the stock conda install |
| RTX 4090, L40S, H100, anything newer (sm_89+) | **2b**, Docker |

### 2a. Stock conda install (Ampere and older)

```bash
git clone https://github.com/RosettaCommons/RFdiffusion ~/software/RFdiffusion
cd ~/software/RFdiffusion

# model weights, ~5 GB
bash scripts/download_models.sh ~/software/RFdiffusion/models

conda env create -f env/SE3nv.yml
conda activate SE3nv
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../..
pip install -e .
```

Note the weights are served over plain **HTTP**, not HTTPS. If you are behind a
proxy that rewrites or blocks HTTP, that download is where it will fail.

### 2b. Docker (modern GPUs, and the more reliable route generally)

Needs Docker plus the NVIDIA Container Toolkit. This sidesteps the whole
dependency problem, because the image carries a stack that works.

```bash
git clone https://github.com/RosettaCommons/RFdiffusion ~/software/RFdiffusion
cd ~/software/RFdiffusion
docker build -f docker/Dockerfile -t rfdiffusion .

mkdir -p ~/rfd/{models,inputs,outputs}
bash scripts/download_models.sh ~/rfd/models
```

Then run with the directories mounted:

```bash
docker run -it --rm --gpus all \
  -v ~/rfd/models:/models -v ~/rfd/inputs:/inputs -v ~/rfd/outputs:/outputs \
  rfdiffusion \
  inference.model_directory_path=/models \
  inference.output_prefix=/outputs/test \
  inference.input_pdb=/inputs/target_loop2.pdb \
  inference.num_designs=2 \
  'contigmap.contigs=[A315-335/0 B315-335/0 20-30]'
```

There is also a prebuilt image at `rosettacommons/rfdiffusion` on Docker Hub if
you would rather not build.

The image is still CUDA 11.6 based. If even that is too old for your card, the
smallest working change is to bump the base image and the `torch`/`dgl` pins in
`docker/Dockerfile` to a matching CUDA — for example `torch==2.1.0+cu121` with
`dgl` built for cu121 — and rebuild. Change both together; mismatched torch and
dgl CUDA builds fail at import with an unhelpful message.

**Verify (either route):**

```bash
conda activate SE3nv            # route 2a only
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
cd ~/software/RFdiffusion/examples && bash design_ppi.sh
```

`design_ppi.sh` is the upstream smoke test and takes a few minutes. If it
produces PDBs in `example_outputs/`, the install is good.

---

## 3. Macrocycles: nothing extra to install

**RFpeptides is part of mainline RFdiffusion.** There is no separate repository
and no second environment. Cyclic design is turned on with two flags:

```
inference.cyclic=True
inference.cyc_chains='a'
```

Two traps, both of which will fail quietly rather than loudly:

- **The peptide contig goes first.** `cyc_chains='a'` cyclises output chain A,
  so the diffused peptide must be the first contig segment and the target comes
  after it. This is the reverse of ordinary binder design. Get it backwards and
  you cyclise part of your target.
- Keep `--config-name base` and `diffuser.T=50`. Both are part of the published
  protocol.

Our `configs/rfpeptides/run_macrocycles.sh` already has this right.

---

## 4. Smoke test against our targets

Once both are in, confirm the pipeline end to end on the smallest target
(`target_loop2.pdb`, 42 residues — it will run on almost anything):

```bash
export RFDIFFUSION=~/software/RFdiffusion
cd /path/to/mtre_binder_design

# RFdiffusion: 2 designs, a minute or two
python $RFDIFFUSION/scripts/run_inference.py \
  inference.output_prefix=./designs/smoke/test \
  inference.input_pdb=./targets/target_loop2.pdb \
  inference.num_designs=2 \
  'contigmap.contigs=[A315-335/0 B315-335/0 55-80]' \
  'ppi.hotspot_res=[A326,A322,A327]'

# then the geometry filter
python3 scripts/filter_designs.py --mode loop2 ./designs/smoke/*.pdb
```

The contig strings for all three targets are in
`targets/mtre_sites.json` under `contigs`, regenerated whenever
`prepare_mtre_targets.py` runs, so they cannot drift out of step with the PDBs.

---

## 5. When it goes wrong

| Symptom | Cause | Fix |
| --- | --- | --- |
| **BC2:** `bindcraft: this environment has no ...` | running under the wrong interpreter | `source .venv/bin/activate` |
| **BC2:** `module not installed: <name>` | wheel resolved but unpacked nothing | `bash install.sh` again in the same environment |
| **BC2:** works on login node, fails on compute nodes | installer read the login node's driver | `bash install.sh cuda12` |
| **BC2:** `checkpoints missing or unfinished` | interrupted download | delete `~/.cache/bindcraft/alphafold`, `bindcraft fetch-weights` |
| **BC2:** out of GPU memory | too many packed workers | `BINDCRAFT_WORKERS_PER_GPU=1` |
| `jax.devices()` shows only CpuDevice | wrong `--cuda` at install | re-run `install_bindcraft.sh` with the version `nvidia-smi` reports |
| `no kernel image is available` | GPU newer than the pinned CUDA | route 2b, Docker |
| `nvrtc: invalid value for --gpu-architecture` | same | route 2b, Docker |
| dgl import error after upgrading torch | torch and dgl built for different CUDA | reinstall both for the same CUDA |
| CUDA out of memory in BindCraft | target too large | trim further — ours are already trimmed; or drop to a two-protomer wedge |
| Weights download fails | `files.ipd.uw.edu` is plain HTTP | check proxy settings |
| BindCraft runs but accepts nothing | normal | several hundred trajectories is expected; thousands on hard targets |

---

## 6. What this does not cover

Neither tool can be installed in a Claude Code cloud session: those containers
have no GPU, and the network policy blocks both GitHub and the weights host.
Install on a lab workstation, a university HPC allocation, or a rented cloud GPU.

GitHub is not an option either. GitHub Actions does offer a GPU runner (one
Tesla T4, 16 GB), but jobs are capped at 6 hours with no extension, repository
cache is 10 GB — less than the two weight sets together — and artifact storage
is 500 MB on Free. It is fine for a smoke test on every push, and unusable for
a real campaign.

One hardware note if you are choosing between national facilities: **both tools
are CUDA-only.** NCI Gadi has NVIDIA GPUs and is fine. Pawsey's Setonix GPU
partition is AMD MI250X running ROCm, so neither pipeline will run there
without a porting effort that is well beyond the scope of this project.
