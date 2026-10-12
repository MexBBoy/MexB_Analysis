#!/usr/bin/env python3
"""Trim a target structure for BindCraft2, and say whether the result fits.

Memory in BC2 grows with the SQUARE of the padded complex, so trimming is not
tidying - it is the difference between a campaign that runs and one that is
refused. This reports the cost before you commit to a file.

Two ways to cut, combinable:

  --drop 21-101,345-349     remove residue ranges (the "loose ends")
  --keep-above -5           keep residues whose CA sits above this z

--keep-above is the useful one for a membrane protein already placed in a
membrane frame (OPM/PPM output, or anything oriented with the membrane normal
on z). OPM puts the bilayer centre at z=0 with the hydrophobic slab at roughly
+-15 A, so for an outer-membrane beta-barrel, z>0 is broadly the extracellular
half and z<0 the periplasmic half. Check the printed profile rather than
assuming - frames differ.

Standard library only: runs under a system python3 with no environment active.

    python3 trim_target.py in.pdb --keep-above -5 -o out.pdb
    python3 trim_target.py in.pdb --drop 21-101,345-349 --binder 96
"""

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict

# BC2's documented per-worker budget, verified against the five figures in
# docs/source/installation.md (7.1 GB at 64 residues through 26.7 at 512).
MODEL_RESIDENT_GB = 3.4
PER_RESIDUE_SQ_GB = 38e-6
HEADROOM_GB = 4.0
LENGTH_BUCKET = 32


def per_worker_gb(n: int) -> float:
    return 2.0 * (MODEL_RESIDENT_GB + PER_RESIDUE_SQ_GB * n * n)


def padded(n: int, bucket: int = LENGTH_BUCKET) -> int:
    return -(-n // bucket) * bucket


def parse_ranges(text: str) -> set[int]:
    out: set[int] = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part.lstrip("-"):
            lo, hi = part.rsplit("-", 1)
            out.update(range(int(lo), int(hi) + 1))
        else:
            out.add(int(part))
    return out


def read_atoms(path: str) -> list[str]:
    with open(path) as fh:
        return fh.readlines()


def residue_key(line: str) -> tuple[str, int]:
    return line[21], int(line[22:26])


def summarise(lines: list[str], label: str) -> int:
    chains: "OrderedDict[str, list[tuple[int, float]]]" = OrderedDict()
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            ch, num = residue_key(line)
            chains.setdefault(ch, []).append((num, float(line[46:54])))
    total = 0
    print(f"\n{label}")
    for ch, residues in chains.items():
        nums = [n for n, _ in residues]
        zs = [z for _, z in residues]
        runs, start, prev = [], nums[0], nums[0]
        for n in nums[1:]:
            if n == prev + 1:
                prev = n
            else:
                runs.append((start, prev))
                start = prev = n
        runs.append((start, prev))
        total += len(residues)
        print(f"  chain {ch}: {len(residues):4d} residues  "
              f"{'/'.join(f'{a}-{b}' for a, b in runs)}  "
              f"z {min(zs):+.0f}..{max(zs):+.0f}")
    print(f"  TOTAL {total} residues across {len(chains)} chain(s)")
    return total


def z_profile(lines: list[str], step: int = 10) -> None:
    zs = [float(l[46:54]) for l in lines
          if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "CA"]
    if not zs:
        return
    lo, hi = int(min(zs) // step) * step, int(max(zs) // step + 1) * step
    print("\nresidues by depth (all chains), 10 A bands:")
    widest = 0
    bands = []
    for z in range(hi, lo - step, -step):
        n = sum(1 for v in zs if z <= v < z + step)
        bands.append((z, n))
        widest = max(widest, n)
    for z, n in bands:
        bar = "#" * int(40 * n / widest) if widest else ""
        print(f"  z {z:+5d}..{z + step:+5d} {n:5d} {bar}")


def verdict(target_residues: int, binder: int) -> None:
    n = target_residues + padded(binder)
    gb = per_worker_gb(n)
    usable = 24.0 - HEADROOM_GB
    workers = int(usable // gb)
    print(f"\nBindCraft2 on a 24 GB card, with a {binder} aa binder "
          f"(pads to {padded(binder)}):")
    print(f"  padded complex N = {n}")
    print(f"  {gb:.1f} GB per worker  ->  ", end="")
    if workers >= 1:
        print(f"FITS, {workers} worker(s)")
    else:
        print(f"DOES NOT FIT ({gb:.0f} GB needed, {usable:.0f} GB available)")
        ceiling = 0
        while per_worker_gb(ceiling + 1 + padded(binder)) <= usable:
            ceiling += 1
        print(f"  trim the target to about {ceiling} residues to fit on one worker")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("-o", "--out", help="write the trimmed structure here")
    ap.add_argument("--drop", default="", help="residue ranges to remove, e.g. 21-101,345-349")
    ap.add_argument("--keep-above", type=float, help="keep residues whose CA z is at or above this")
    ap.add_argument("--keep-below", type=float, help="keep residues whose CA z is at or below this")
    ap.add_argument("--chains", default="", help="keep only these chains, e.g. A,B")
    ap.add_argument("--binder", type=int, default=96, help="binder length for the memory estimate")
    args = ap.parse_args()

    lines = read_atoms(args.pdb)
    if not any(l.startswith("ATOM") for l in lines):
        sys.exit(f"{args.pdb} holds no ATOM records - is it mmCIF, or an error page?")

    summarise(lines, f"BEFORE  {args.pdb}")
    z_profile(lines)

    drop = parse_ranges(args.drop) if args.drop else set()
    keep_chains = {c.strip() for c in args.chains.split(",") if c.strip()}

    # Decide per residue from its CA, then apply that to all of its atoms, so a
    # z cut never keeps half a side chain.
    ca_z = {residue_key(l): float(l[46:54]) for l in lines
            if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "CA"}

    def keep(line: str) -> bool:
        if not line.startswith(("ATOM", "HETATM")):
            return not line.startswith(("ANISOU", "TER"))
        key = residue_key(line)
        if keep_chains and key[0] not in keep_chains:
            return False
        if key[1] in drop:
            return False
        z = ca_z.get(key)
        if z is None:
            return False
        if args.keep_above is not None and z < args.keep_above:
            return False
        if args.keep_below is not None and z > args.keep_below:
            return False
        return True

    kept = [l for l in lines if keep(l)]
    total = summarise(kept, "AFTER")
    verdict(total, args.binder)

    if args.out:
        with open(args.out, "w") as fh:
            fh.writelines(l for l in kept if l.startswith(("ATOM", "HETATM")))
            fh.write("END\n")
        print(f"\nwrote {args.out}")
    else:
        print("\n(no -o given, nothing written)")


if __name__ == "__main__":
    main()
