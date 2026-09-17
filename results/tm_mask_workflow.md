# Making a TM-domain mask for 3DVA and 3D Flex

MexB, amphipol-exchanged. Residue ranges are this project's own TM definition
(`REGIONS` in `scripts/mexb_common.py`), so they can be checked against the
rest of the analysis.

## 0. Two masks, not one

These jobs use a mask for different purposes, and the same file will not do
for both.

| | what the mask is for | what it should cover |
|---|---|---|
| **3DVA** | a focus mask: variability is solved only inside it | the TM barrel only |
| **3D Flex** | the domain the tetrahedral mesh is built over | the whole protein |

For 3D Flex, masking down to the TM alone is a mistake. The mesh deforms
whatever it spans and the reconstruction uses the whole particle, so a
TM-only mesh leaves the porter and docking domains unmodelled. If the
question is "how does the TM move", build the mesh over the whole protein
and read the TM motion out of the result - optionally declaring the docking
domain rigid in Mesh Prep, rather than deleting it from the mask.

What both masks share is that they should **exclude the amphipol belt**. That
density is disordered and low-resolution; inside a 3DVA focus mask it
dominates the leading components with belt wobble, and inside a Flex mesh it
absorbs deformation that belongs to the protein. Building the mask from the
fitted atomic model rather than by thresholding the map excludes it for free.

## 1. TM residue ranges

Per protomer, 407 residues:

```
10-35, 337-495, 516-565, 859-1030
```

made up of TM1 (10-35), TM2 (337-359), Ialpha (360-380), TM3-6 (381-495),
TM6b (516-565), the junction (859-875) and TM7-12 (876-1030).

Everything else - the porter and docking domains - is:

```
1-9, 36-336, 496-515, 566-858, 1031-1046
```

Mask all three chains. The three protomers sit in different states, so the
map was refined in C1; a mask around one protomer of a pseudo-symmetric
trimer invites alignment ambiguity in 3DVA.

## 2. Export the consensus volume from CryoSPARC

From the refinement the model was built into (Homogeneous or NU-Refine),
download `*_volume_map.mrc`. Note its **box size and pixel size** - the mask
must match both exactly.

## 3. ChimeraX: make the shape

ChimeraX makes the *shape*; CryoSPARC makes it into a *soft mask*. Doing the
softening in CryoSPARC avoids hand-tuning thresholds here.

```
open consensus.mrc                          # becomes #1
open fitted_model.pdb                       # becomes #2

# synthetic density from the TM atoms only, on the map's own grid
molmap #2/A,B,C:10-35,337-495,516-565,859-1030 6 onGrid #1
                                            # becomes #3

save tm_shape.mrc model #3
```

`onGrid #1` is the part that matters: without it molmap picks its own box and
origin and the mask will not align to the map in CryoSPARC.

The `6` is the resolution in A at which the atoms are blurred. 5-8 A gives a
shape that follows the helices without being lumpy; lower makes a mask that
hugs each helix and can clip density at the edges.

Check before saving: `volume #3 level 0.1` and confirm it envelopes the TM
density of #1 with a little room, and that no amphipol density pokes out.

For the **3D Flex** mask, repeat with the whole chain instead:

```
molmap #2/A,B,C 6 onGrid #1
save protein_shape.mrc model #3
```

## 4. CryoSPARC: shape to soft mask

1. **Import Volumes** - bring in `tm_shape.mrc`. Import as a volume, not a
   mask, since it still needs thresholding.
2. **Volume Tools** - inputs: the imported volume.
   - Type of input volume: `volume`
   - Type of output volume: `mask`
   - Threshold: pick from the ChimeraX level that looked right (with molmap
     output, something around 0.05-0.2 usually)
   - Dilation radius (px): 3-6
   - Soft padding width (px): 6-10

   Dilation puts the hard edge outside the density; soft padding is the
   cosine falloff. Too tight and you cut signal and create edge artefacts;
   too loose and the amphipol comes back in. Inspect the output.

3. Confirm the output box and pixel size match the particles.

## 5. 3DVA

**3D Variability Analysis** - inputs: particles from the consensus
refinement, plus the TM mask.

- Number of modes to solve: 3 to start
- Filter resolution: 5-8 A. This is the main knob. Too high and the modes
  chase noise; for TM helix motion 6 A is a reasonable start.
- Leave symmetry at C1.

Then **3D Variability Display**:

- `simple` for the raw component volumes
- `intermediates` for a movie along each mode - best for seeing helix motion
- `cluster` to split particles into discrete classes you can refine
  separately

Sanity checks: a mode that is mostly mask-edge ripple means the mask is too
tight; a mode that is mostly belt motion means the amphipol got in.

## 6. 3D Flex

1. **3D Flex Data Prep** - particles + consensus volume. Set the training box
   size down (a crop to ~128-192 px is usual) or training is very slow.
2. **3D Flex Mesh Prep** - consensus volume + the **whole-protein** mask from
   step 3, not the TM mask.
   - Base number of tetra cells: start ~20, raise for finer motion
   - Use the segmentation options to mark the docking domain rigid if the
     interest is TM motion specifically
3. **3D Flex Training** - 2 latent dimensions to start. The rigidity prior is
   the parameter to tune: too low and the map tears, too high and nothing
   moves.
4. **3D Flex Reconstruction** - produces the improved map.
5. **3D Flex Generator** - volumes sampled along the latent space, for
   figures.

## 7. Checking the result is real

The same trap as any variability analysis: both methods will always return
something. Worth having:

- Run 3DVA with the mask on **half the particles**, then the other half, and
  confirm the leading modes describe the same motion.
- Compare the amplitude of the motion against the same-state protomer
  scatter measured here: per-helix centroid shifts between two protomers in
  the same assigned state run 0.33-1.09 A, median 0.66 A. A 3DVA mode moving
  helices by less than about 1 A is within that floor.
- The TM7 displacement between assigned states is 2.8 A, so a genuine
  state-related mode should be of that order, not tenfold larger.
