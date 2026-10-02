# Installing RFdiffusion and BindCraft locally

Verified against the repositories as of 2 October 2026:
RosettaCommons/RFdiffusion `86507b6`, martinpacesa/BindCraft `7713aa0`.

Budget about 2 hours, most of it downloads.

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

## 1. BindCraft

The primary generator. Install this one first — it is the better behaved of the
two and it is where most of our design runs will happen.

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
Install on your own machine, Spartan, or a rented cloud GPU.
