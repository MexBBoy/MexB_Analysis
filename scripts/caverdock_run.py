#!/usr/bin/env python3
"""Pull a ligand through a MexB tunnel with CaverDock 1.2.

CAVER and this pipeline both describe the tunnel geometrically. CaverDock
adds the energetics: it docks the ligand at successive discs along the
tunnel and returns lower- and upper-bound energy profiles, so a geometric
constriction can be checked against an actual energy barrier.

CaverDock is not redistributed here. It is downloaded from
loschmidt.chemi.muni.cz on first use, together with the Ubuntu 18.04
runtime libraries its binary needs.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import CXDIR, TABLES, WORK_DIR, Structure, coords, fmt, \
    load_structures, write_csv

CD = os.path.join(WORK_DIR, "caverdock-1.2")
LIBS = os.path.join(WORK_DIR, "cdlibs", "ex", "usr")



CD_URL = ("https://loschmidt.chemi.muni.cz/static/releases/caverdock/1.2/"
          "caverdock-1.2-ubuntu18.04.tar.gz")
# CaverDock 1.2 ships a single Ubuntu 18.04 binary. On any newer distribution
# its Boost/OpenMPI/hwloc sonames are missing, so those runtime libraries are
# fetched from the Ubuntu archive and used via LD_LIBRARY_PATH rather than
# installed system-wide.
RUNTIME_DEBS = [
    ("http://archive.ubuntu.com/ubuntu/pool/main/b/boost1.65.1/"
     "libboost-system1.65.1_1.65.1+dfsg-0ubuntu5_amd64.deb"),
    ("http://archive.ubuntu.com/ubuntu/pool/main/b/boost1.65.1/"
     "libboost-thread1.65.1_1.65.1+dfsg-0ubuntu5_amd64.deb"),
    ("http://archive.ubuntu.com/ubuntu/pool/main/b/boost1.65.1/"
     "libboost-serialization1.65.1_1.65.1+dfsg-0ubuntu5_amd64.deb"),
    ("http://archive.ubuntu.com/ubuntu/pool/main/b/boost1.65.1/"
     "libboost-program-options1.65.1_1.65.1+dfsg-0ubuntu5_amd64.deb"),
    ("http://archive.ubuntu.com/ubuntu/pool/universe/o/openmpi/"
     "libopenmpi2_2.1.1-8_amd64.deb"),
    ("http://archive.ubuntu.com/ubuntu/pool/universe/o/openmpi/"
     "openmpi-common_2.1.1-8_all.deb"),
    ("http://archive.ubuntu.com/ubuntu/pool/universe/o/openmpi/"
     "openmpi-bin_2.1.1-8_amd64.deb"),
    ("https://launchpad.net/ubuntu/+archive/primary/+files/"
     "libhwloc5_1.11.9-1_amd64.deb"),
]


def ensure_caverdock():
    """Download CaverDock and the runtime it needs, if not already present."""
    import platform
    import shutil as _sh
    import tarfile
    import urllib.request

    binary = os.path.join(CD, "bin", "caverdock")
    if not os.path.exists(binary):
        if platform.system() != "Linux" or platform.machine() != "x86_64":
            raise SystemExit(
                "CaverDock 1.2 is distributed only as a Linux x86-64 binary. "
                "On another platform use the Singularity image from "
                "caver.cz (caverdock-1.2.sif) and point CD at it.")
        os.makedirs(WORK_DIR, exist_ok=True)
        tgz = os.path.join(WORK_DIR, "caverdock-1.2.tar.gz")
        if not os.path.exists(tgz):
            print(f"  downloading CaverDock from {CD_URL}")
            urllib.request.urlretrieve(CD_URL, tgz)
        with tarfile.open(tgz) as t:
            t.extractall(WORK_DIR)
        os.chmod(binary, 0o755)
        os.chmod(os.path.join(CD, "bin", "caverdock_split"), 0o755)

    if not os.path.exists(os.path.join(LIBS, "bin", "mpirun.openmpi")):
        if _sh.which("dpkg-deb") is None:
            raise SystemExit("dpkg-deb is needed to unpack CaverDock's "
                             "runtime libraries")
        d = os.path.join(WORK_DIR, "cdlibs")
        os.makedirs(d, exist_ok=True)
        for url in RUNTIME_DEBS:
            deb = os.path.join(d, os.path.basename(url))
            if not os.path.exists(deb):
                print(f"  fetching {os.path.basename(url)}")
                urllib.request.urlretrieve(url, deb)
            subprocess.run(["dpkg-deb", "-x", deb, os.path.join(d, "ex")],
                           check=True)
        print("  CaverDock runtime staged in work/cdlibs")
    return binary


def env():
    e = dict(os.environ)
    lib = os.path.join(LIBS, "lib", "x86_64-linux-gnu")
    e["LD_LIBRARY_PATH"] = (f"{lib}:{os.path.join(lib,'openmpi','lib')}:"
                            + e.get("LD_LIBRARY_PATH", ""))
    e["OPAL_PREFIX"] = LIBS
    e["PATH"] = os.path.join(LIBS, "bin") + ":" + e.get("PATH", "")
    # the container has no IB fabric; keep OpenMPI on shared memory
    e["OMPI_MCA_btl"] = "self,vader,tcp"
    return e


def read_trace(path):
    pts, rad = [], []
    for line in open(path):
        if line.startswith("HETATM"):
            pts.append((float(line[30:38]), float(line[38:46]),
                        float(line[46:54])))
            rad.append(float(line[60:66]))
    return np.array(pts), np.array(rad)


def write_dsd(path, pts, rad, spacing=0.5, max_len=None):
    """CaverDock tunnel: x y z dx dy dz radius, one disc per line."""
    s = np.concatenate([[0], np.cumsum(
        np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    if max_len:
        keep = s <= max_len
        pts, rad, s = pts[keep], rad[keep], s[keep]
    want = np.arange(0, s[-1], spacing)
    P = np.stack([np.interp(want, s, pts[:, i]) for i in range(3)], axis=1)
    R = np.interp(want, s, rad)
    D = np.gradient(P, axis=0)
    D /= (np.linalg.norm(D, axis=1)[:, None] + 1e-9)
    with open(path, "w") as fh:
        for p, d, r in zip(P, D, R):
            fh.write(f"{p[0]:.4f} {p[1]:.4f} {p[2]:.4f} "
                     f"{d[0]:.6f} {d[1]:.6f} {d[2]:.6f} {r:.4f}\n")
    return P, R


def write_receptor(path, atoms):
    tmp = path.replace(".pdbqt", "_tmp.pdb")
    with open(tmp, "w") as fh:
        for i, a in enumerate(atoms, start=1):
            nm = a.name if len(a.name) >= 4 else f" {a.name:<3.3s}"
            fh.write(f"ATOM  {i:5d} {nm:<4.4s} {a.resname:>3.3s} "
                     f"{a.chain}{a.resseq:4d}    "
                     f"{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}"
                     f"  1.00  0.00          {a.element:>2.2s}\n")
        fh.write("END\n")
    subprocess.run(["obabel", tmp, "-O", path, "-xr", "--partialcharge",
                    "gasteiger"], check=True, capture_output=True)
    return path


def write_ligand(path, atoms):
    tmp = path.replace(".pdbqt", "_tmp.pdb")
    with open(tmp, "w") as fh:
        for i, a in enumerate(atoms, start=1):
            nm = a.name if len(a.name) >= 4 else f" {a.name:<3.3s}"
            fh.write(f"HETATM{i:5d} {nm:<4.4s} {a.resname:>3.3s} "
                     f"{a.chain}{a.resseq:4d}    "
                     f"{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}"
                     f"  1.00  0.00          {a.element:>2.2s}\n")
        fh.write("END\n")
    subprocess.run(["obabel", tmp, "-O", path, "-h", "--partialcharge",
                    "gasteiger"], check=True, capture_output=True)
    return path


def parse_profile(path, spacing):
    """CaverDock trajectory PDBQT.

    Each disc carries a line
        REMARK CAVERDOCK TUNNEL: <disc> <energy> <radius> <internal>
    The fourth field is CaverDock's own coordinate, not arc length in A, so
    position along the tunnel is recomputed from the disc index and the disc
    spacing we generated the .dsd with.
    """
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("REMARK CAVERDOCK TUNNEL:"):
            p = line.split()
            try:
                disc = int(float(p[3]))
                out.append({"disc": disc, "energy": float(p[4]),
                            "radius": float(p[5]),
                            "arc": disc * spacing})
            except (ValueError, IndexError):
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", default="Amp_MexB_20260826")
    ap.add_argument("--chain", default="E")
    ap.add_argument("--seed", default="ZZ72000")
    ap.add_argument("--rank", default="t1")
    ap.add_argument("--max-len", type=float, default=32.0,
                    help="trim the tunnel to this arc length (A)")
    ap.add_argument("--spacing", type=float, default=0.5)
    ap.add_argument("--exhaustiveness", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=1)
    ap.add_argument("--parse-only", action="store_true",
                    help="re-parse an existing run without redocking")
    ap.add_argument("--ranks", type=int, default=4,
                    help="MPI ranks: 1 master + workers, minimum 2")
    a = ap.parse_args()

    binary = ensure_caverdock()

    s = next(x for x in load_structures() if x.name == a.structure)
    d = os.path.join(WORK_DIR, "caverdock_runs",
                     f"{a.structure}_{a.chain}_{a.seed}_{a.rank}")
    os.makedirs(d, exist_ok=True)

    trace = os.path.join(
        CXDIR, f"{a.structure}_protein_{a.chain}_{a.seed}_{a.rank}"
               f"_tunnel.pdb")
    pts, rad = read_trace(trace)
    P, R = write_dsd(os.path.join(d, "tunnel.dsd"), pts, rad,
                     spacing=a.spacing, max_len=a.max_len)
    print(f"  tunnel: {len(P)} discs over {a.max_len:.0f} A, "
          f"radius {R.min():.2f}-{R.max():.2f} A")

    prot = [x for x in s.protein_atoms if not x.is_hydrogen]
    lig = [x for x in s.het_atoms if not x.is_hydrogen]
    if not a.parse_only:
        print("  preparing receptor (this takes a minute)")
        write_receptor(os.path.join(d, "receptor.pdbqt"), prot)
        write_ligand(os.path.join(d, "ligand.pdbqt"), lig)

    lo, hi = P.min(axis=0), P.max(axis=0)
    cen = (lo + hi) / 2
    size = (hi - lo) + 12.0
    conf = os.path.join(d, "caverdock.conf")
    with open(conf, "w") as fh:
        fh.write(f"receptor = receptor.pdbqt\nligand = ligand.pdbqt\n"
                 f"tunnel = tunnel.dsd\n\n"
                 f"center_x = {cen[0]:.3f}\ncenter_y = {cen[1]:.3f}\n"
                 f"center_z = {cen[2]:.3f}\n"
                 f"size_x = {size[0]:.1f}\nsize_y = {size[1]:.1f}\n"
                 f"size_z = {size[2]:.1f}\n\n"
                 f"exhaustiveness = {a.exhaustiveness}\ncpu = {a.cpu}\n"
                 f"multiple_search = 3\n")
    print(f"  box {size[0]:.0f} x {size[1]:.0f} x {size[2]:.0f} A "
          f"centred on the tunnel")

    if a.parse_only:
        print("  --parse-only: reusing the existing trajectory")
    else:
        print("  running CaverDock ...")
    # CaverDock refuses to run a tunnel job on a single rank: it needs a
    # master plus at least one worker, so it must be launched under mpirun.
    mpirun = os.path.join(LIBS, "bin", "mpirun.openmpi")
    r = None
    cmd = [mpirun, "-n", str(a.ranks), "--allow-run-as-root",
           "--oversubscribe", binary, "--config", "caverdock.conf",
           "--out", "traj"]
    if not a.parse_only:
        r = subprocess.run(cmd, cwd=d, env=env(), capture_output=True,
                           text=True, timeout=21600)
    log = os.path.join(d, "caverdock.log")
    if r is not None:
        open(log, "w").write(r.stdout + "\n---stderr---\n" + r.stderr)
        if r.returncode != 0:
            # the upper-bound stage often fails where the lower bound
            # succeeded; a partial result is still worth parsing
            print(f"  CaverDock exited {r.returncode} (see {log}) - parsing "
                  f"whatever trajectory was produced")

    rows = []
    for kind in ("lb", "ub"):
        prof = parse_profile(os.path.join(d, f"traj-{kind}.pdbqt"),
                             a.spacing)
        if not prof:
            continue
        E = [p["energy"] for p in prof]
        imax = int(np.argmax(E))
        imin = int(np.argmin([p["radius"] for p in prof]))
        print(f"  {kind}: {len(prof)} discs, E {min(E):.2f} to {max(E):.2f} "
              f"kcal/mol; barrier from the site {max(E) - E[0]:.2f}")
        print(f"       energy max at {prof[imax]['arc']:.1f} A "
              f"(radius {prof[imax]['radius']:.2f} A)")
        print(f"       narrowest point at {prof[imin]['arc']:.1f} A "
              f"(radius {prof[imin]['radius']:.2f} A, "
              f"E {prof[imin]['energy']:.2f})")
        for p in prof:
            rows.append([a.structure, a.chain, a.seed, kind, p["disc"],
                         fmt(p["arc"], 2), fmt(p["radius"], 2),
                         fmt(p["energy"], 2)])
    if not rows:
        print("  no profile parsed")
        return
    write_csv(os.path.join(TABLES, "caverdock_profile.csv"),
              ["structure", "chain", "seed", "bound", "disc",
               "position_along_tunnel_A", "tunnel_radius_A",
               "energy_kcal_mol"], rows)
    print(f"  wrote results/tables/caverdock_profile.csv ({len(rows)} rows)")


if __name__ == "__main__":
    main()
