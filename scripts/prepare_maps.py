#!/usr/bin/env python3
"""Cut cryo-EM maps down to something you can drag into a chat window.

ONE FILE, ONE COMMAND, needs only numpy. Copy this script to wherever the
full-size maps live and run:

    python3 prepare_maps.py

It finds every .mrc/.ccp4/.map next to it (or in --maps), matches each to a
PDB model by filename, cuts small full-resolution boxes around the regions
this analysis actually needs, converts them to 16-bit, and packs the lot into
a single .zip you can upload.

Typical result: a few hundred MB of maps becomes a handful of MB, with no
loss of detail where it matters.

    python3 prepare_maps.py --budget-mb 20      # aim for a smaller bundle
    python3 prepare_maps.py --radius 20         # bigger boxes
    python3 prepare_maps.py --context           # also add a binned whole map
"""
from __future__ import annotations

import argparse
import glob
import os
import struct
import sys
import zipfile

import numpy as np

# --------------------------------------------------------------- constants
# MexB (UniProt P52002) numbering, zero offset.
DBP = [136, 139, 178, 277, 279, 327, 573, 610, 612, 615, 617, 626, 628, 630]
PBP = [79, 128, 151, 152, 176, 180, 273, 274, 276, 668, 672, 674, 676,
       717, 819, 825, 828]
RELAY = [407, 408, 939, 971, 976]
SWITCH = list(range(613, 623))
SOLVENT = {"HOH", "WAT", "DOD"}
AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
      "MSE"}


# ------------------------------------------------------------------ MRC I/O
# A minimal reader/writer so no pip install is needed. The header layout is
# the CCP4/MRC2014 standard; the only subtle part is mapc/mapr/maps, which
# say which crystallographic axis each array axis is - assuming (z,y,x)
# silently shifts every coordinate for maps that are not stored that way.

class MRC:
    def __init__(self, path):
        with open(path, "rb") as fh:
            hdr = fh.read(1024)
            if len(hdr) < 1024:
                raise ValueError(f"{path}: truncated header")
            nx, ny, nz, mode = struct.unpack_from("<4i", hdr, 0)
            nxs, nys, nzs = struct.unpack_from("<3i", hdr, 16)
            mx, my, mz = struct.unpack_from("<3i", hdr, 28)
            cella = struct.unpack_from("<3f", hdr, 40)
            mapc, mapr, maps = struct.unpack_from("<3i", hdr, 64)
            origin = struct.unpack_from("<3f", hdr, 196)
            nsymbt = struct.unpack_from("<i", hdr, 92)[0]
            dt = {0: np.int8, 1: np.int16, 2: np.float32,
                  6: np.uint16, 12: np.float16}.get(mode)
            if dt is None:
                raise ValueError(f"{path}: unsupported MRC mode {mode}")
            fh.seek(1024 + max(nsymbt, 0))
            data = np.frombuffer(fh.read(nx * ny * nz * np.dtype(dt).itemsize),
                                 dtype=dt).astype(np.float32)
        data = data.reshape(nz, ny, nx)            # (slow, medium, fast)
        self.mode = mode
        self.voxel = np.array([cella[0] / max(mx, 1), cella[1] / max(my, 1),
                               cella[2] / max(mz, 1)], dtype=float)
        axes = [maps, mapr, mapc]                  # crystal axis per array axis
        order = [axes.index(i) for i in (1, 2, 3)]
        self.data = np.transpose(data, order)      # now (x, y, z)
        nstart = np.array([nxs, nys, nzs], dtype=float)
        org = np.array(origin, dtype=float)
        self.origin = org if np.any(org) else nstart * self.voxel
        self.shape = np.array(self.data.shape)

    def world_to_grid(self, xyz):
        return (np.asarray(xyz, float) - self.origin) / self.voxel

    def grid_to_world(self, ijk):
        return np.asarray(ijk, float) * self.voxel + self.origin


def write_mrc(path, data_xyz, voxel, origin, as_int16=True,
              parent_stats=None):
    """data_xyz is (x, y, z); written back in MRC's (z, y, x) order."""
    d = np.asarray(data_xyz, dtype=np.float32)
    arr = np.transpose(d, (2, 1, 0))
    if as_int16:
        lo, hi = float(arr.min()), float(arr.max())
        rng = (hi - lo) or 1.0
        # 16 bits over the local range is far finer than map noise
        scaled = np.round((arr - lo) / rng * 60000.0 - 30000.0)
        out = scaled.astype("<i2")
        mode = 1
        amin, amax, amean = lo, hi, float(arr.mean())
        rms = float(arr.std())
    else:
        out = arr.astype("<f4")
        mode = 2
        amin, amax, amean = float(arr.min()), float(arr.max()), float(arr.mean())
        rms = float(arr.std())
    nz, ny, nx = arr.shape
    hdr = bytearray(1024)
    struct.pack_into("<4i", hdr, 0, nx, ny, nz, mode)
    struct.pack_into("<3i", hdr, 16, 0, 0, 0)                  # nstart
    struct.pack_into("<3i", hdr, 28, nx, ny, nz)               # mx,my,mz
    struct.pack_into("<3f", hdr, 40, nx * voxel[0], ny * voxel[1],
                     nz * voxel[2])                            # cella
    struct.pack_into("<3f", hdr, 52, 90.0, 90.0, 90.0)         # cellb
    struct.pack_into("<3i", hdr, 64, 1, 2, 3)                   # mapc/r/s
    struct.pack_into("<3f", hdr, 76, amin, amax, amean)
    struct.pack_into("<i", hdr, 92, 0)                          # nsymbt
    struct.pack_into("<3f", hdr, 196, *[float(v) for v in origin])
    hdr[208:212] = b"MAP "
    struct.pack_into("<i", hdr, 212, 0x00004144)                # little endian
    struct.pack_into("<f", hdr, 216, rms)
    labels = []
    if as_int16:
        labels.append(f"CROP int16 lo={lo:.6g} hi={hi:.6g}")
    if parent_stats is not None:
        # a crop is mostly protein, so its own mean is not the map's mean;
        # keep the parent statistics or every z-score comes out too low
        labels.append(f"PARENT mean={parent_stats[0]:.6g} "
                      f"sigma={parent_stats[1]:.6g}")
    struct.pack_into("<i", hdr, 220, len(labels))
    for i, note in enumerate(labels[:10]):
        off = 224 + 80 * i
        hdr[off:off + min(len(note), 80)] = note[:80].encode()
    with open(path, "wb") as fh:
        fh.write(bytes(hdr))
        fh.write(out.tobytes())


# ------------------------------------------------------------------ PDB
def read_pdb(path):
    atoms = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        el = line[76:78].strip().upper()
        name = line[12:16].strip()
        if not el:
            el = name[0] if name[:1].isalpha() else name[1:2]
        if el == "H":
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        atoms.append({"name": name, "resname": line[17:20].strip(),
                      "chain": line[21], "resseq": int(line[22:26]),
                      "xyz": xyz, "het": line.startswith("HETATM")})
    return atoms


def sel(atoms, chain=None, resids=None, het=None, protein=None):
    out = atoms
    if chain is not None:
        out = [a for a in out if a["chain"] == chain]
    if resids is not None:
        rs = set(resids)
        out = [a for a in out if a["resseq"] in rs]
    if het is True:
        out = [a for a in out if a["het"] and a["resname"] not in SOLVENT]
    if protein is True:
        out = [a for a in out if a["resname"] in AA]
    return out


def centre(atoms):
    return np.array([a["xyz"] for a in atoms], dtype=float).mean(axis=0)


# ------------------------------------------------------------------ crops
def crop(m, cen, radius):
    c = m.world_to_grid(cen)
    r = np.ceil(radius / m.voxel).astype(int)
    lo = np.maximum(np.floor(c - r).astype(int), 0)
    hi = np.minimum(np.ceil(c + r).astype(int) + 1, m.shape)
    if np.any(hi <= lo):
        return None, None
    return m.data[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]], m.grid_to_world(lo)


def bin_map(m, factor):
    s = (m.shape // factor) * factor
    d = m.data[:s[0], :s[1], :s[2]]
    d = d.reshape(s[0] // factor, factor, s[1] // factor, factor,
                  s[2] // factor, factor).mean(axis=(1, 3, 5))
    return d, m.voxel * factor


def targets_for(atoms):
    """(label, centre) for every region worth keeping, per chain."""
    out = []
    chains = sorted({a["chain"] for a in atoms if a["resname"] in AA})
    ligs = {}
    for a in sel(atoms, het=True):
        ligs.setdefault((a["chain"], a["resseq"], a["resname"]), []).append(a)
    for (ch, rs, rn), ats in sorted(ligs.items()):
        out.append((f"ligand_{rn}{ch}{rs}", centre(ats)))
    for ch in chains:
        for label, resids in (("relay", RELAY), ("pockets", DBP + PBP),
                              ("switch", SWITCH)):
            a = sel(atoms, chain=ch, resids=resids, protein=True)
            if a:
                out.append((f"{label}_{ch}", centre(a)))
    return out


def match_structure(map_path, pdbs):
    """Pair a map with a model by longest shared filename prefix."""
    stem = os.path.basename(map_path).lower()
    best, score = None, 0
    for p in pdbs:
        ps = os.path.basename(p).lower().rsplit(".", 1)[0]
        n = 0
        for a, b in zip(stem, ps):
            if a != b:
                break
            n += 1
        # a shared token also counts, e.g. "amp" or "ddm"
        for tok in ("amp", "ddm", "lmt", "zz7"):
            if tok in stem and tok in ps:
                n += 10
        if n > score:
            best, score = p, n
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maps", nargs="*", default=None,
                    help="map files or a directory (default: alongside this "
                         "script and ./maps)")
    ap.add_argument("--pdb", nargs="*", default=None,
                    help="model files (default: any .pdb found nearby)")
    ap.add_argument("--radius", type=float, default=15.0)
    ap.add_argument("--budget-mb", type=float, default=45.0)
    ap.add_argument("--context", action="store_true",
                    help="also include a heavily binned whole map")
    ap.add_argument("--float32", action="store_true",
                    help="keep full float precision (roughly doubles size)")
    ap.add_argument("--out", default="mexb_map_crops.zip")
    a = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    cands = []
    if a.maps:
        for m in a.maps:
            cands += (sorted(glob.glob(os.path.join(m, "*")))
                      if os.path.isdir(m) else [m])
    else:
        for d in (here, os.getcwd(), os.path.join(here, "maps"),
                  os.path.join(os.getcwd(), "maps")):
            cands += glob.glob(os.path.join(d, "*"))
    maps = sorted({p for p in cands
                   if p.lower().endswith((".mrc", ".ccp4", ".map"))})
    pdbs = a.pdb or sorted({p for d in (here, os.getcwd(),
                                        os.path.join(here, "structures"),
                                        os.path.join(os.getcwd(), "structures"))
                            for p in glob.glob(os.path.join(d, "*.pdb"))})
    if not maps:
        raise SystemExit("No .mrc/.ccp4/.map files found. Pass --maps.")
    if not pdbs:
        raise SystemExit("No .pdb models found. Pass --pdb.")
    print(f"maps:   {len(maps)} found")
    print(f"models: {len(pdbs)} found\n")

    work = os.path.join(os.getcwd(), "_map_crops")
    os.makedirs(work, exist_ok=True)
    written = []
    for mp in maps:
        pdb = match_structure(mp, pdbs)
        if pdb is None:
            print(f"  {os.path.basename(mp)}: no matching model, skipped")
            continue
        try:
            m = MRC(mp)
        except Exception as e:
            print(f"  {os.path.basename(mp)}: unreadable ({e})")
            continue
        full_mb = m.data.nbytes / 1e6
        print(f"  {os.path.basename(mp)}  {tuple(int(x) for x in m.shape)} "
              f"voxels, {m.voxel[0]:.2f} A/vox, {full_mb:.0f} MB "
              f"-> model {os.path.basename(pdb)}")
        pstats = (float(m.data.mean()), float(m.data.std()))
        atoms = read_pdb(pdb)
        stem = os.path.basename(mp).rsplit(".", 1)[0]
        for label, cen in targets_for(atoms):
            sub, org = crop(m, cen, a.radius)
            if sub is None or sub.size == 0:
                continue
            out = os.path.join(work, f"{stem}__{label}.mrc")
            write_mrc(out, sub, m.voxel, org, as_int16=not a.float32,
                      parent_stats=pstats)
            written.append(out)
        if a.context:
            factor = max(2, int(np.ceil((m.shape.prod() / 4e6) ** (1 / 3))))
            d, vox = bin_map(m, factor)
            out = os.path.join(work, f"{stem}__context_bin{factor}.mrc")
            write_mrc(out, d, vox, m.origin, as_int16=not a.float32,
                      parent_stats=pstats)
            written.append(out)

    if not written:
        raise SystemExit("Nothing was written - check the map/model pairing.")

    with zipfile.ZipFile(a.out, "w", zipfile.ZIP_DEFLATED,
                         compresslevel=6) as z:
        for f in written:
            z.write(f, os.path.basename(f))
    mb = os.path.getsize(a.out) / 1e6
    print(f"\n  {len(written)} crops -> {a.out}  ({mb:.1f} MB)")
    if mb > a.budget_mb:
        print(f"  Larger than the {a.budget_mb:.0f} MB budget. Re-run with "
              f"--radius {max(8, a.radius - 5):.0f}, or drop --context.")
    else:
        print("  Small enough to upload. Attach this single file.")
    print(f"  (intermediate crops left in {work}; safe to delete)")


if __name__ == "__main__":
    main()
