# Running this analysis locally

This pipeline was built in a remote Claude Code session, which has no access
to your filesystem. Running it on your own machine removes the upload problem
entirely: the maps stay where they are and the scripts read them directly.

## 1. Get the code

```bash
git clone https://github.com/MexBBoy/MexB_Analysis.git
cd MexB_Analysis
git checkout claude/mexb-substrate-tunnel-analysis-fdkv8o
```

## 2. Dependencies

```bash
pip install numpy scipy pyKVFinder openpyxl matplotlib mrcfile
```

Optional, all fetched automatically on first use if missing:

| tool | needed for | notes |
| --- | --- | --- |
| CAVER 3.0.3 | tunnel cross-check | Java required; downloaded by `run_caver.py` |
| CaverDock 1.2 | transport energetics | Linux x86-64 only; `caverdock_run.py` stages its Ubuntu 18.04 runtime |
| fpocket | pocket druggability | build from source into `work/fpocket`; skipped if absent |
| OpenBabel | PDBQT preparation | `apt install openbabel` |

## 3. Run it

```bash
bash run.sh                 # everything
SKIP_CAVER=1 bash run.sh    # skip the ~25 min CAVER cross-check
```

Regenerates every table, figure, ChimeraX session, the 3D viewer, the
manuscript draft and the workbook, from a clean `work/`.

## 4. Add the density maps

Put them anywhere and point at them — no cropping, no uploading:

```bash
python3 scripts/map_validation.py \
    --map /path/to/Amp_sharpened.mrc \
    --structure Amp_MexB_20260826 \
    --resolution 3.0
```

Writes `results/tables/map_validation_<structure>.csv`, scoring per-residue
and per-ligand RSCC plus per-atom map z-scores. It targets the flags the
coordinate analysis could not settle:

- **R971** — is the modelled rotamer supported, or is it a modelling error?
- **the three DDM molecules** — are all three real, or is the 98.2% occlusion
  figure resting on an over-modelled one?
- **the ampicillin pose** — how much confidence does the contact list deserve?
- **N676 / N718 / L827** — how well resolved is the tunnel constriction?

Half-maps are better than a single sharpened map here: sharpening can flatter
a weak side chain, which is exactly the R971 question.

`scripts/prepare_maps.py` is only needed for sending maps to a *remote*
session. Locally it is unnecessary.

## What is already established

- 51/53 checks against `PROTOCOL.md` section 6 pass (`results/tables/validation.csv`).
- The two failures are the protocol's own tunnel constants. This pipeline and
  CAVER 3.0.3 independently give 2.21 / 2.2175 Å with the same constriction
  residues, against the protocol's expected 2.01 Å at F615. Re-derive that
  constant before trusting it.
- CaverDock's lower-bound barrier peaks ~14 Å before the geometric
  bottleneck, so radius alone picks the wrong rate-limiting residues.
