"""Shared constants, PDB I/O and geometry for the MexB analysis pipeline.

Constants come from PROTOCOL.md section 0 and must not be re-derived.
"""
from __future__ import annotations

import os
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT_DIR = os.path.join(REPO, "structures")
WORK_DIR = os.path.join(REPO, "work")
TABLES = os.path.join(REPO, "results", "tables")
FIGURES = os.path.join(REPO, "results", "figures")
CXDIR = os.path.join(REPO, "results", "chimerax")
DATA = os.path.join(REPO, "data")

for _d in (WORK_DIR, TABLES, FIGURES, CXDIR):
    os.makedirs(_d, exist_ok=True)


# ---------------------------------------------------------------- constants

def _rng(*pairs):
    out = []
    for a, b in pairs:
        out.extend(range(a, b + 1))
    return out


DOMAINS = {
    "porter": _rng((40, 130), (136, 180), (278, 281), (284, 328),
                   (571, 666), (679, 717), (814, 818), (821, 858)),
    "docking": _rng((181, 277), (718, 813)),
    "TM": _rng((10, 35), (337, 495), (516, 565), (876, 1030)),          # trimmed
    "TM_untrimmed": _rng((10, 35), (337, 359), (361, 565), (862, 1030)),
}

SUBDOMAINS = {
    "PN1": _rng((40, 130)),
    "PN2": _rng((136, 180), (278, 281), (284, 328)),
    "PC1": _rng((571, 666)),
    "PC2": _rng((679, 717), (814, 818), (821, 858)),
}

DBP = [136, 139, 178, 277, 279, 327, 573, 610, 612, 615, 617, 626, 628, 630]
PBP = [79, 128, 151, 152, 176, 180, 273, 274, 276, 668, 672, 674, 676,
       717, 819, 825, 828]
SWITCH_LOOP = list(range(613, 623))

RELAY = {407: "ASP", 408: "ASP", 939: "LYS", 971: "ARG", 976: "THR"}
RELAY_ATOMS = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "LYS": ("NZ",),
    "ARG": ("NE", "NH1", "NH2"),
    "THR": ("OG1",),
}

REFERENCE_STATES = {          # PROTOCOL section 0, angstrom
    "Access":    (26.5, 26.2),
    "Binding":   (28.3, 29.1),
    "Extrusion": (29.8, 24.5),
}

# Regions for stage-3 per-region summary.
REGIONS = {
    "TM1": _rng((10, 35)),
    "PN1": SUBDOMAINS["PN1"],
    "PN2": SUBDOMAINS["PN2"],
    "DN": _rng((181, 277)),
    "TM2": _rng((337, 359)),
    "Ialpha": _rng((360, 380)),
    "TM3-6": _rng((381, 495)),
    "loop496-515": _rng((496, 515)),
    "TM6b": _rng((516, 565)),
    "PC1": SUBDOMAINS["PC1"],
    "PC2": SUBDOMAINS["PC2"],
    "DC": _rng((718, 813)),
    "junction859-875": _rng((859, 875)),
    "TM7-12": _rng((876, 1030)),
}

AA3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "ASH": "D", "GLH": "E", "LYN": "K",
}

AROMATIC = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "CG", "CD1", "NE1"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
}

# Kyte-Doolittle hydropathy
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
    "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
    "Y": -1.3, "V": 4.2,
}

# van der Waals radii (angstrom) for grid/tunnel work
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
       "H": 1.20, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}

# Protomers whose own density does not support the model, established by
# map validation on 2026-08-30 (PROTOCOL known issues 2 and 3). Chain F of
# the ampicillin reconstruction has whole-protomer median RSCC 0.182, 21st
# percentile of its map, with relay residues scoring negative. Anything
# derived from it - relay distances, state assignment - is not determinable.
UNSUPPORTED_PROTOMERS = {("Amp_MexB_20260826", "F"): (
    "median RSCC 0.182 (21st pct of its map); relay residues negative "
    "(THR976F -0.164, LYS939F -0.081, ARG971F -0.004)")}


def density_warning(structure, chain):
    return UNSUPPORTED_PROTOMERS.get((structure, chain))


SOLVENT_HET = {"HOH", "WAT", "DOD"}
DETERGENTS = {"LMT", "LDA", "DDM", "BOG", "OLC", "PLM", "MC3", "PEE",
              "PGT", "CLR", "C8E", "DMU", "TRD", "UND", "D10", "D12"}


def vdw(element: str) -> float:
    return VDW.get(element.strip().upper(), 1.70)


# ------------------------------------------------------------------- PDB I/O

class Atom:
    __slots__ = ("serial", "name", "resname", "chain", "resseq", "icode",
                 "xyz", "occ", "bfac", "element", "hetatm")

    def __init__(self, serial, name, resname, chain, resseq, icode, xyz,
                 occ, bfac, element, hetatm):
        self.serial = serial
        self.name = name
        self.resname = resname
        self.chain = chain
        self.resseq = resseq
        self.icode = icode
        self.xyz = xyz
        self.occ = occ
        self.bfac = bfac
        self.element = element
        self.hetatm = hetatm

    @property
    def is_hydrogen(self):
        return self.element == "H"

    def __repr__(self):
        return (f"<Atom {self.chain}/{self.resname}{self.resseq}/{self.name}>")


class Structure:
    """A parsed PDB: flat atom list plus convenience indexes."""

    def __init__(self, path):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.atoms: list[Atom] = []
        self._parse()

    def _parse(self):
        with open(self.path) as fh:
            for line in fh:
                rec = line[:6]
                if rec not in ("ATOM  ", "HETATM"):
                    continue
                element = line[76:78].strip().upper()
                name = line[12:16].strip()
                if not element:            # fall back to the atom name
                    element = name[0] if name[0].isalpha() else name[1]
                try:
                    xyz = (float(line[30:38]), float(line[38:46]),
                           float(line[46:54]))
                except ValueError:
                    continue
                self.atoms.append(Atom(
                    serial=int(line[6:11]),
                    name=name,
                    resname=line[17:20].strip(),
                    chain=line[21].strip() or " ",
                    resseq=int(line[22:26]),
                    icode=line[26].strip(),
                    xyz=xyz,
                    occ=float(line[54:60] or 1.0),
                    bfac=float(line[60:66] or 0.0),
                    element=element,
                    hetatm=(rec == "HETATM"),
                ))

    # --- selections -------------------------------------------------------
    @property
    def protein_atoms(self):
        return [a for a in self.atoms
                if not a.hetatm and a.resname in AA3to1]

    @property
    def het_atoms(self):
        return [a for a in self.atoms
                if a.hetatm and a.resname not in SOLVENT_HET]

    @property
    def chains(self):
        seen = []
        for a in self.protein_atoms:
            if a.chain not in seen:
                seen.append(a.chain)
        return seen

    def ligands(self):
        """[(chain, resseq, resname, [atoms]), ...] for non-solvent HETATM."""
        groups = {}
        for a in self.het_atoms:
            groups.setdefault((a.chain, a.resseq, a.resname), []).append(a)
        return [(c, r, n, ats) for (c, r, n), ats in sorted(groups.items())]

    def ca(self, chain):
        """{resseq: xyz} for CA atoms of one protein chain."""
        out = {}
        for a in self.protein_atoms:
            if a.chain == chain and a.name == "CA":
                out[a.resseq] = np.array(a.xyz)
        return out

    def sequence(self, chain):
        """{resseq: one-letter} from the model (CA presence defines a residue)."""
        out = {}
        for a in self.protein_atoms:
            if a.chain == chain and a.name == "CA":
                out[a.resseq] = AA3to1.get(a.resname, "X")
        return out

    def residue_atoms(self, chain, resseq, heavy_only=True):
        return [a for a in self.atoms
                if a.chain == chain and a.resseq == resseq
                and not a.hetatm and (not heavy_only or not a.is_hydrogen)]

    def chain_atoms(self, chain, heavy_only=True, protein_only=True):
        src = self.protein_atoms if protein_only else self.atoms
        return [a for a in src if a.chain == chain
                and (not heavy_only or not a.is_hydrogen)]


def coords(atoms):
    return np.array([a.xyz for a in atoms], dtype=float)


# ---------------------------------------------------------------- geometry

def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Return (R, t) minimising |R @ mobile.T + t - target.T|."""
    mc = mobile.mean(axis=0)
    tc = target.mean(axis=0)
    P = mobile - mc
    Q = target - tc
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = tc - R @ mc
    return R, t


def apply_rt(R, t, X):
    return (R @ X.T).T + t


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def rotation_angle_axis(R):
    """Rotation angle (deg) and unit axis of a proper rotation matrix."""
    tr = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    ang = float(np.degrees(np.arccos(tr)))
    w, v = np.linalg.eig(R)
    idx = int(np.argmin(np.abs(w - 1.0)))
    axis = np.real(v[:, idx])
    n = np.linalg.norm(axis)
    axis = axis / n if n else np.array([0.0, 0.0, 1.0])
    return ang, axis


def superpose_on(mob: Structure, mob_chain, ref: Structure, ref_chain,
                 resids):
    """Kabsch fit of mob chain onto ref chain over CA of `resids` present in
    both. Returns (R, t, n_used, fit_rmsd)."""
    mca, rca = mob.ca(mob_chain), ref.ca(ref_chain)
    common = [r for r in resids if r in mca and r in rca]
    if len(common) < 3:
        raise ValueError("too few common residues for superposition")
    M = np.array([mca[r] for r in common])
    T = np.array([rca[r] for r in common])
    R, t = kabsch(M, T)
    return R, t, len(common), rmsd(apply_rt(R, t, M), T)


def centroid(ca_map, resids):
    pts = [ca_map[r] for r in resids if r in ca_map]
    if not pts:
        return None
    return np.array(pts).mean(axis=0)


def load_reference_sequence():
    path = os.path.join(DATA, "P52002.fasta")
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


def load_structures():
    """All PDBs under structures/, sorted by name."""
    paths = sorted(os.path.join(STRUCT_DIR, f)
                   for f in os.listdir(STRUCT_DIR) if f.endswith(".pdb"))
    return [Structure(p) for p in paths]


def write_csv(path, header, rows):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    return path


def fmt(x, nd=2):
    if x is None:
        return ""
    if isinstance(x, float):
        if np.isnan(x):
            return ""
        return f"{x:.{nd}f}"
    return x
