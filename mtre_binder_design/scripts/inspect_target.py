#!/usr/bin/env python3
"""Print what a PDB/mmCIF target actually contains: chains, residue ranges,
counts and sequence. Standard library only - runs on any python3."""
import sys, re
from collections import OrderedDict

THREE_TO_ONE = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

path = sys.argv[1]
chains = OrderedDict()
with open(path) as fh:
    for line in fh:
        if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA':
            ch, num, res = line[21], int(line[22:26]), line[17:20].strip()
            chains.setdefault(ch, []).append((num, THREE_TO_ONE.get(res, 'X')))

if not chains:
    sys.exit(f"no CA atoms found in {path} - is it mmCIF? convert to PDB, or tell me")

total = 0
for ch, residues in chains.items():
    nums = [n for n, _ in residues]
    runs, start, prev = [], nums[0], nums[0]
    for n in nums[1:]:
        if n == prev + 1: prev = n
        else: runs.append((start, prev)); start = prev = n
    runs.append((start, prev))
    seq = ''.join(a for _, a in residues)
    total += len(residues)
    print(f"chain {ch}: {len(residues)} residues, "
          f"numbered {'/'.join(f'{a}-{b}' for a, b in runs)}"
          f"{'  <-- GAPS' if len(runs) > 1 else ''}")
    for i in range(0, len(seq), 60):
        print(f"   {nums[i]:>5} {seq[i:i+60]}")
print(f"\nTOTAL {total} residues across {len(chains)} chain(s): {','.join(chains)}")
