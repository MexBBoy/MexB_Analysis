#!/usr/bin/env python3
"""Shared cryo-EM map handling: reading, world<->grid transforms, cropping.

Kept separate from map_validation.py so the crop utility can be run on a
machine that only has numpy and mrcfile installed.
"""
from __future__ import annotations

import numpy as np

try:
    import mrcfile
except ImportError:  # pragma: no cover
    mrcfile = None


class Map:
    """An MRC/CCP4 map with a correct grid <-> world transform.

    MRC stores the array in (slow, medium, fast) axis order, which is not
    necessarily (z, y, x): mapc/mapr/maps say which crystallographic axis
    each array axis corresponds to. Getting this wrong silently shifts every
    density value, so it is handled explicitly here rather than assumed.
    """

    def __init__(self, path):
        if mrcfile is None:
            raise SystemExit("pip install mrcfile")
        self.path = path
        with mrcfile.open(path, permissive=True, mode="r") as m:
            data = np.asarray(m.data, dtype=np.float32)
            h = m.header
            self.voxel = np.array([float(m.voxel_size.x),
                                   float(m.voxel_size.y),
                                   float(m.voxel_size.z)])
            mapc, mapr, maps = int(h.mapc), int(h.mapr), int(h.maps)
            nstart = np.array([int(h.nxstart), int(h.nystart),
                               int(h.nzstart)])
            origin = np.array([float(h.origin.x), float(h.origin.y),
                               float(h.origin.z)])
            self.labels = []
            try:
                for i in range(int(h.nlabl)):
                    t = h.label[i].tobytes().decode(errors="replace").strip()
                    if t:
                        self.labels.append(t)
            except Exception:
                pass
        # array axes are (slow, medium, fast) = (maps, mapr, mapc)
        # transpose so that axis 0 -> x, 1 -> y, 2 -> z
        axes = [maps, mapr, mapc]              # crystallographic axis per array axis
        order = [axes.index(i) for i in (1, 2, 3)]
        self.data = np.transpose(data, order)
        # nstart is given in crystallographic axis order already
        self.origin = origin + nstart * self.voxel
        if np.allclose(origin, 0) and np.allclose(nstart, 0):
            self.origin = np.zeros(3)
        # crops written by prepare_maps.py are stored as int16 over a
        # recorded range; undo that here or every density value - and every
        # z-score derived from it - is in the wrong units
        for t in getattr(self, "labels", []):
            if t.startswith("CROP int16"):
                try:
                    lo = float(t.split("lo=")[1].split()[0])
                    hi = float(t.split("hi=")[1].split()[0])
                except (IndexError, ValueError):
                    break
                self.data = ((self.data + 30000.0) / 60000.0
                             * (hi - lo) + lo).astype(np.float32)
                break
        self.shape = np.array(self.data.shape)

    def world_to_grid(self, xyz):
        return (np.asarray(xyz, dtype=float) - self.origin) / self.voxel

    def grid_to_world(self, ijk):
        return np.asarray(ijk, dtype=float) * self.voxel + self.origin

    def sample(self, xyz):
        """Trilinear interpolation at world coordinates."""
        g = self.world_to_grid(np.atleast_2d(xyz))
        lo = np.floor(g).astype(int)
        frac = g - lo
        out = np.full(len(g), np.nan, dtype=np.float64)
        ok = np.all((lo >= 0) & (lo < self.shape - 1), axis=1)
        if not ok.any():
            return out
        l, f = lo[ok], frac[ok]
        acc = np.zeros(len(l))
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    w = (((1 - f[:, 0]) if dx == 0 else f[:, 0])
                         * ((1 - f[:, 1]) if dy == 0 else f[:, 1])
                         * ((1 - f[:, 2]) if dz == 0 else f[:, 2]))
                    acc += w * self.data[l[:, 0] + dx, l[:, 1] + dy,
                                         l[:, 2] + dz]
        out[ok] = acc
        return out

    @property
    def parent_stats(self):
        """Mean and sigma of the map this was cropped from, if recorded.

        A crop is mostly protein, so its own mean sits well above the parent
        map's; using it would push every z-score down by a constant.
        """
        for t in getattr(self, "labels", []):
            if t.startswith("PARENT"):
                try:
                    mu = float(t.split("mean=")[1].split()[0])
                    sd = float(t.split("sigma=")[1].split()[0])
                    return mu, sd
                except (IndexError, ValueError):
                    return None
        return None

    def stats(self):
        p = self.parent_stats
        if p:
            return p
        d = self.data
        return float(d.mean()), float(d.std())


def crop(map_obj, centre, radius, out_path):
    """Write a cubic sub-map of `radius` A about `centre`."""
    if mrcfile is None:
        raise SystemExit("pip install mrcfile")
    c = map_obj.world_to_grid(centre)
    r = np.ceil(radius / map_obj.voxel).astype(int)
    lo = np.maximum(np.floor(c - r).astype(int), 0)
    hi = np.minimum(np.ceil(c + r).astype(int) + 1, map_obj.shape)
    sub = map_obj.data[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    new_origin = map_obj.grid_to_world(lo)
    with mrcfile.new(out_path, overwrite=True) as m:
        # write back in (z, y, x) array order, the MRC default
        m.set_data(np.transpose(sub, (2, 1, 0)).astype(np.float32))
        m.voxel_size = tuple(float(v) for v in map_obj.voxel)
        m.header.origin.x = float(new_origin[0])
        m.header.origin.y = float(new_origin[1])
        m.header.origin.z = float(new_origin[2])
        m.update_header_from_data()
    return sub.shape, out_path
