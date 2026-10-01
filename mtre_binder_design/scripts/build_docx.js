// Builds the MtrE binder-design write-up as a Word document.
//   node scripts/build_docx.js [outfile]
const fs = require('fs');
const path = require('path');
const d = require('docx');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  ImageRun, Footer, PageNumber, LevelFormat, ExternalHyperlink,
} = d;

const ROOT = path.resolve(__dirname, '..');
const OUT = process.argv[2] || path.join(ROOT, 'MtrE_binder_design_plan.docx');

// A4 portrait, default 1" margins -> usable width in DXA.
const CONTENT = 9026;
const GREY = 'F2F3F5';
const RULE = 'BFC4CB';

const p = (text, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 140, line: 276 },
  alignment: opts.align,
  children: [new TextRun({ text, size: opts.size ?? 21, italics: opts.italics, bold: opts.bold, color: opts.color })],
});

const rich = (children, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 140, line: 276 },
  children,
});

const t = (text, o = {}) => new TextRun({ text, size: o.size ?? 21, bold: o.bold, italics: o.italics, color: o.color });

const link = (text, url) => new ExternalHyperlink({
  children: [new TextRun({ text, size: 21, style: 'Hyperlink' })],
  link: url,
});

const h1 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_1, spacing: { before: 320, after: 160 },
  children: [new TextRun({ text, size: 28, bold: true, color: '1A1D21' })],
});

const h2 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_2, spacing: { before: 240, after: 120 },
  children: [new TextRun({ text, size: 23, bold: true, color: '2E3338' })],
});

const numberedReasons = (x) => numbered(x, 'reasons');
const numberedPipeline = (x) => numbered(x, 'pipeline');
const numberedCascade = (x) => numbered(x, 'cascade');

const bullet = (text) => new Paragraph({
  numbering: { reference: 'dots', level: 0 },
  spacing: { after: 90, line: 276 },
  children: [new TextRun({ text, size: 21 })],
});

const numbered = (text, ref = 'reasons') => new Paragraph({
  numbering: { reference: ref, level: 0 },
  spacing: { after: 90, line: 276 },
  children: [new TextRun({ text, size: 21 })],
});

const rule = () => new Paragraph({
  spacing: { before: 60, after: 160 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: RULE } },
  children: [],
});

// --- table helper -----------------------------------------------------------
function table(headers, rows, widths) {
  const cell = (text, { head = false, w } = {}) => new TableCell({
    width: { size: w, type: WidthType.DXA },
    margins: { top: 80, bottom: 80, left: 110, right: 110 },
    shading: head ? { type: ShadingType.CLEAR, fill: GREY, color: 'auto' } : undefined,
    children: String(text).split(' ').map((line) => new Paragraph({
      spacing: { after: 0, line: 252 },
      children: [new TextRun({ text: line, size: 19, bold: head })],
    })),
  });
  return new Table({
    width: { size: CONTENT, type: WidthType.DXA },
    columnWidths: widths,
    borders: {
      top: { style: BorderStyle.SINGLE, size: 4, color: RULE },
      bottom: { style: BorderStyle.SINGLE, size: 4, color: RULE },
      left: { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' },
      right: { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 2, color: RULE },
      insideVertical: { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' },
    },
    rows: [
      new TableRow({
        tableHeader: true,
        children: headers.map((hd, i) => cell(hd, { head: true, w: widths[i] })),
      }),
      ...rows.map((r) => new TableRow({ children: r.map((c, i) => cell(c, { w: widths[i] })) })),
    ],
  });
}

// --- figure -----------------------------------------------------------------
const figPath = path.join(ROOT, 'figures', 'fig1_targets.png');
const figure = new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { before: 120, after: 80 },
  children: [new ImageRun({
    type: 'png',
    data: fs.readFileSync(figPath),
    transformation: { width: 600, height: 616 },
  })],
});

const caption = new Paragraph({
  spacing: { after: 220, line: 252 },
  children: [
    t('Figure 1. ', { bold: true, size: 19 }),
    t('What each campaign is trying to achieve. ', { bold: true, size: 19 }),
    t('(a) MtrE in cross-section, drawn from the coordinates of PDB 4MT0 rather than schematically: '
      + 'the wall profile, the depths and the open lumen are measured. The barrel lumen is an open funnel '
      + 'from the extracellular mouth down to z = +34, and the channel’s only constriction — the '
      + 'aspartate ring D422/D425 — lies at z = −46, on the periplasmic side and out of reach. '
      + '(b–d) Mock-ups of the three binder concepts on the same scale. Grey arrow, direction of drug '
      + 'efflux; red cross, efflux blocked. Only (b) and (c) occlude the conduit; (d) engages the surface '
      + 'and would not be expected to stop efflux on its own.', { size: 19 }),
  ],
});

// --- document ---------------------------------------------------------------
const doc = new Document({
  creator: 'MtrE binder design',
  title: 'De novo binder design against MtrE',
  styles: {
    default: { document: { run: { font: 'Calibri', size: 21 } } },
    characterStyles: [{
      id: 'Hyperlink', name: 'Hyperlink',
      run: { color: '1F6FB4', underline: { type: 'single' } },
    }],
  },
  numbering: {
    config: [
      {
        reference: 'dots',
        levels: [{
          level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 360, hanging: 200 } } },
        }],
      },
      ...['reasons', 'pipeline', 'cascade'].map((reference) => ({
        reference,
        levels: [{
          level: 0, format: LevelFormat.DECIMAL, text: '%1.', alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 360, hanging: 240 } } },
        }],
      })),
    ],
  },
  sections: [{
    properties: {},
    footers: {
      default: new Footer({
        children: [new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [new TextRun({ children: [PageNumber.CURRENT], size: 18, color: '7A8089' })],
        })],
      }),
    },
    children: [
      new Paragraph({
        spacing: { after: 60 },
        children: [new TextRun({
          text: 'De novo binder design against MtrE',
          size: 36, bold: true, color: '1A1D21',
        })],
      }),
      new Paragraph({
        spacing: { after: 40 },
        children: [new TextRun({
          text: 'An alternative to Boltz-led screening, and what the structure says we should target',
          size: 24, color: '2E3338',
        })],
      }),
      new Paragraph({
        spacing: { after: 220 },
        children: [new TextRun({ text: '1 October 2026', size: 19, color: '7A8089' })],
      }),
      rule(),

      h1('Summary'),
      p('Boltz is producing too many false positives to screen with, and the reason is structural rather than '
        + 'a matter of thresholds. It is a structure predictor, not a binding discriminator, and it has no lipid '
        + 'bilayer — so on an outer-membrane protein it confidently docks binders onto surfaces that are '
        + 'buried in membrane in a real cell. Tightening confidence cut-offs does not remove those, because the '
        + 'model is confident and self-consistent about them.'),
      p('The proposal is to stop using a co-folding model as the primary screen and move to a generate-then-filter '
        + 'pipeline that has been benchmarked experimentally, with Boltz demoted to one vote among three. Three '
        + 'things follow from measuring the MtrE structure directly:'),
      bullet('The β-barrel lumen is an open funnel. Pore radius never drops below 7.4 Å between the '
        + 'extracellular mouth and z = +34 — about 32 Å of channel with no constriction. A binder applied '
        + 'to intact cells can insert that far without crossing the membrane. This is a far better site than the '
        + 'surface crown, which protrudes only 7–10 Å above the lipid.'),
      bullet('The only gate in the channel is unreachable. The aspartate ring sits 110 Å below the surface on '
        + 'the periplasmic side. It is the obvious “cork” site and we should not spend any effort on it.'),
      bullet('BindCraft should be the primary generator, not RFdiffusion and not Boltz. On a single target with one '
        + 'experimental readout, BindCraft returned a 50% hit rate from 48 designs against 24% from 96 RFdiffusion '
        + 'designs.'),
      p('Target files, generator configs and the filtering code are written and in the repository. What is needed '
        + 'is GPU time.'),

      h1('Why the current approach keeps producing false positives'),
      p('This is a known and documented failure mode, not something we have done wrong. The BindCraft authors '
        + 'tested AlphaFold3 as a filter and still found a large proportion of false positives, and Boltz is the '
        + 'same class of model. Three problems compound on this particular target:'),
      numberedReasons('No bilayer. The model will dock a binder onto the lipid-facing belt of the β-barrel — a '
        + 'hydrophobic, highly designable-looking patch that is physically occluded by the outer membrane. A large '
        + 'share of our false positives are probably exactly this.'),
      numberedReasons('Flexible loops. The extracellular loops are short and mobile, so per-residue confidence is '
        + 'uninformative there and pLDDT/PAE cannot be read the way they are read on a folded domain.'),
      numberedReasons('No negative controls. Confidence scores only mean something relative to a decoy distribution, and '
        + 'we have not been generating one.'),
      p('The fix is not a better threshold. It is to put a deterministic geometric filter in front of the scoring '
        + 'models, and to score every surviving design against decoys as well as against MtrE.'),

      h1('What the structure says'),
      p('Everything below was measured from PDB 4MT0 (MtrE open state, 3.29 Å, trimeric) by a script that '
        + 're-checks each assumption on every run and refuses to write target files if any check fails.'),
      figure,
      caption,

      h2('A numbering trap worth knowing about'),
      p('4MT0 uses precursor numbering, which includes the 20-residue signal peptide, so author numbering = mature '
        + 'numbering + 20. The aspartate ring is called D402/D405 in the papers and D422/D425 in the coordinates '
        + '— these are the same two residues. Every residue identifier in our configs is author numbering, '
        + 'matching the PDB. A design run launched in the wrong frame targets a site 20 residues away and returns '
        + 'confident, worthless designs.'),

      h1('The three campaigns, and what each is trying to achieve'),
      p('These are not three ways of doing the same thing. They differ in what a success would actually be, and '
        + 'that distinction should drive which one we prioritise.'),

      h2('A — Lumen plug: block efflux by filling the conduit'),
      p('Objective: occlude the MtrE channel so the MtrCDE pump cannot expel antibiotics, turning the binder into '
        + 'an adjuvant that restores sensitivity to drugs the strain is currently resistant to.'),
      p('Why this site: the lumen is open for 32 Å from the outside with no gate, and the mouth is about 15 Å '
        + 'across — wide enough for a single α-helix with side chains, but not for a helical hairpin. That '
        + 'constraint independently reproduces the colicin E1 / TolC cryo-EM structure, in which colicin plugs TolC '
        + 'as a single-pass folded helix. Isolated colicin E1 fragments potentiate aminoglycosides, fluoroquinolones '
        + 'and macrolides by this mechanism, and klebicin C blocks drug efflux through TolC the same way.'),
      p('Success looks like: reduced Nile red or ethidium efflux, and a 2–64-fold MIC shift in a checkerboard '
        + 'against azithromycin, ciprofloxacin and ceftriaxone, lost in a ΔmtrE strain.'),
      p('Main risk: demanding geometry, and a mini-protein at a mucosal surface is exposed to proteolysis. This is '
        + 'the highest-payoff and highest-risk arm.'),

      h2('B — Macrocycle: block efflux with something small and synthetic'),
      p('Objective: the same mechanism as A, but in a chemistry that is cheap to make, protease-hardenable and '
        + 'avoids a biologic development path. Two sites are worth running — the lumen mouth (which can block '
        + 'efflux) and the Loop 2 groove (which probably cannot, but is a smaller, better-defined pocket).'),
      p('Why this is tractable: the sites are extracellular, so membrane permeability — normally the constraint '
        + 'that kills macrocycle programmes — does not apply. There is also direct precedent in this exact '
        + 'system: rationally designed self-inhibitory peptides against MtrC, MtrD and MtrE gave 2- to 64-fold '
        + 'increases in susceptibility across FA1090 and WHO K/P/X with no human-cell toxicity.'),
      p('Success looks like: the same functional readouts as A, from an 8–16mer that can be made by solid-phase '
        + 'synthesis.'),
      p('Main risk: macrocycle design is less mature than minibinder design, and the published tool is built for '
        + '7–20mers.'),

      h2('C — Crown groove: bind the surface of live cells'),
      p('Objective: a high-affinity, selective handle on the gonococcal surface. This is a targeting and detection '
        + 'reagent, not necessarily an efflux inhibitor.'),
      p('Why this site is the safest place to start: accessibility is already proven. Loop 2 is a 13-residue '
        + 'surface-exposed epitope conserved in more than 98% of gonococcal isolates; antibodies against it bind '
        + 'live gonococci and are bactericidal, and both a Loop 2 peptide vaccine and an anti-Loop 2 monoclonal '
        + 'confer complement-dependent protection. We also found that Loop 2 of one protomer packs directly against '
        + 'Loop 1 of its neighbour (W333 to L114, 6.4 Å), giving a composite inter-protomer groove that buries '
        + 'far more surface than either loop alone — which is the structural answer to the flexible-loop '
        + 'problem.'),
      p('Success looks like: binding to live FA1090 by flow cytometry, lost in ΔmtrE.'),
      p('What it will probably not do: stop efflux. The loops are not the conduit. If we want both, a C-type binder '
        + 'can be fused to an A-type plug.'),

      table(
        ['', 'A — lumen plug', 'B — macrocycle', 'C — crown groove'],
        [
          ['Site', 'Barrel lumen, z +34 to +66', 'Lumen mouth and Loop 2 groove', 'Loop1/Loop2 seam, two protomers'],
          ['Hotspots', 'A123, A312, A311, A132', 'A111, A120, A123; A326, A322, A327', 'A326, A333, C114, C119'],
          ['Size', '60–110 aa', '8–16 aa, cyclic', '55–100 aa'],
          ['Generator', 'BindCraft + RFdiffusion', 'RFpeptides', 'BindCraft + RFdiffusion'],
          ['Blocks efflux?', 'Yes, by design', 'Mouth arm only', 'Unlikely'],
          ['Risk', 'High', 'Medium', 'Low'],
        ],
        [1400, 2542, 2542, 2542],
      ),
      p('Recommended order: C first, because proven accessibility makes it the fastest way to validate the whole '
        + 'assay cascade; B in parallel, because it is cheap; A once the cascade works, because it has the highest '
        + 'payoff and the hardest geometry.', { after: 200 }),

      h1('The proposed pipeline'),
      p('Generation uses BindCraft as the primary tool with RFdiffusion as a fold-diversity arm, following the '
        + 'protocol of Clement et al. (2026), who benchmarked both on one target with one experimental readout. '
        + 'Filtering runs cheapest-first, so we only spend scoring compute on designs that are geometrically '
        + 'possible:'),
      numberedPipeline('Geometry. Deterministic and free. Rejects any design placing atoms below the outer-leaflet '
        + 'boundary, and enforces real insertion depth for plug designs. Written and tested.'),
      numberedPipeline('AlphaFold2 initial-guess, interface pAE < 10.'),
      numberedPipeline('Rosetta: ΔΔG, shape complementarity ≥ 0.62, buried unsatisfied polars ≤ 2.'),
      numberedPipeline('Boltz-2 as a third orthogonal vote — never alone, and never as the primary screen.'),
      numberedPipeline('Decoy panel: rescore survivors against E. coli TolC and P. aeruginosa OprM, and discard anything '
        + 'that scores as well on a decoy. This also generates the selectivity data we need anyway.'),
      numberedPipeline('Scramble control: shuffle interface residues and rescore. If scrambles score nearly as well, the '
        + 'metric is reading fold quality rather than interface quality.'),
      p('Our filter set also tightens BindCraft’s defaults — i_pTM 0.50 to 0.65, i_pAE 0.35 to 0.28, shape '
        + 'complementarity 0.55 to 0.62, interface residues 7 to 12. This cuts accepted-design yield substantially, '
        + 'which is the intention: compute is far cheaper than 96 wells of protein expression.'),

      h1('Experimental cascade'),
      p('Target production is the hard part. MtrE should be expressed with a signal peptide for outer-membrane '
        + 'targeting, purified in DDM or LDAO, and reconstituted into nanodiscs — this is what the colicin E1 / '
        + 'TolC structural work used, and detergent micelles give artefact-prone BLI. Biotinylate the scaffold, not '
        + 'MtrE.'),
      numberedCascade('BLI on MtrE-nanodisc, 200 s association and 500 s dissociation, 96-well.'),
      numberedCascade('Empty-nanodisc counter-screen on every binder. This is the single most important false-positive '
        + 'control in the project.'),
      numberedCascade('TolC-nanodisc counter-screen for selectivity.'),
      numberedCascade('Whole-cell flow cytometry on live FA1090 and WHO K/P/X, with an isogenic ΔmtrE control.'),
      numberedCascade('Nile red or ethidium bromide efflux accumulation assay.'),
      numberedCascade('Antibiotic checkerboard MIC, matched to the published self-inhibitory-peptide benchmark so the '
        + 'numbers are directly comparable.'),
      numberedCascade('Thermal stability and aggregation; then cryo-EM of the lead complex.'),
      p('Two orthogonal primary assays, not one. In the Clement study, BLI response and ligand competition '
        + 'correlated at r = −0.55 to −0.72, and it was that cross-check that separated real binders from '
        + 'artefacts.'),

      h1('What we need'),
      p('Ready now, in the repository: the target preparation script with its self-checks, three trimmed target '
        + 'structures (42, 90 and 147 residues), BindCraft target and filter configs for all three campaigns, '
        + 'RFdiffusion and RFpeptides run scripts, the geometry filter, and a protocol document recording the '
        + 'measured ground truth.'),
      p('Needed: GPU time. BindCraft and RFdiffusion both require CUDA and cannot run on CPU. The authors recommend '
        + '32 GB of GPU memory for large targets; our trimmed targets sit comfortably below that. Spartan would do, as would a short cloud-GPU rental for a first pass.'),

      h1('Corrections and caveats'),
      p('Two elements of the first draft of this plan were wrong, and both were caught by measuring the structure '
        + 'rather than reasoning from the literature. First, the aspartate ring was proposed as the macrocycle '
        + '“cork” target; it is periplasmic and unreachable, and that arm now points at the lumen mouth. '
        + 'Second, a C3-symmetric binder was proposed via RFdiffusion’s symmetric mode; that mode builds '
        + 'symmetric oligomers de novo and does not design a C3 binder onto a C3 target. For three-fold avidity we '
        + 'will design a monomeric groove binder and trimerise it experimentally with a foldon fusion.'),
      bullet('LOS shielding. The crown protrudes only 7–10 Å above the outer-leaflet aromatic girdle, and '
        + 'the lipooligosaccharide inner core extends roughly 10–15 Å above the headgroups, so part of that '
        + 'shell is sterically contested on a live cell. The anti-Loop 2 antibody data says the region is reachable '
        + 'in practice. The signature of a problem here would be designs that bind purified protein and fail on '
        + 'whole cells.'),
      bullet('Gating state. 4MT0 is the open state, which MtrC binding stabilises. A plug that only engages the open '
        + 'channel would be conditionally active on actively-effluxing cells — acceptable for an adjuvant, but '
        + 'it changes the assay design.'),
      bullet('No experimental MtrCDE assembly structure exists. We can model it by homology to AcrAB–TolC, but '
        + 'only to reason about gating, never as a design target.'),
      bullet('Framing. These are efflux-pump adjuvants rather than standalone antimicrobials — that is what the '
        + 'precedents actually deliver. N. gonorrhoeae work is BSL-2.'),

      h1('Key references'),
      rich([t('Clement J, Lkhagvajav T, Hoare BL, et al. ', {}), t('A complete RXFP1–relaxin interaction model '
        + 'unlocks the design of potent mini-protein modulators. ', { italics: true }),
        t('bioRxiv 2026.06.19.733483. '), link('doi.org/10.64898/2026.06.19.733483', 'https://doi.org/10.64898/2026.06.19.733483')], { after: 110 }),
      rich([t('Pacesa M, et al. '), t('One-shot design of functional protein binders with BindCraft. ', { italics: true }),
        t('Nature, 2025. '), link('nature.com/articles/s41586-025-09429-6', 'https://www.nature.com/articles/s41586-025-09429-6')], { after: 110 }),
      rich([t('Rettie SA, et al. '), t('Accurate de novo design of high-affinity protein-binding macrocycles using '
        + 'deep learning (RFpeptides). ', { italics: true }), t('Nature Chemical Biology, 2025. '),
        link('nature.com/articles/s41589-025-01929-w', 'https://www.nature.com/articles/s41589-025-01929-w')], { after: 110 }),
      rich([t('Lei H-T, et al. '), t('Crystal structure of the open state of the Neisseria gonorrhoeae MtrE outer '
        + 'membrane channel. ', { italics: true }), t('PLOS ONE, 2014 (PDB 4MT0). '),
        link('PMC4046963', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4046963/')], { after: 110 }),
      rich([t('Housden NG, Grinter R, Kleanthous C, et al. '), t('Toxin import through the antibiotic efflux channel '
        + 'TolC (klebicin C). ', { italics: true }), t('Nature Communications, 2021. '),
        link('nature.com/articles/s41467-021-24930-y', 'https://www.nature.com/articles/s41467-021-24930-y')], { after: 110 }),
      rich([t('Budiardjo SJ, et al. '), t('Colicin E1 opens its hinge to plug TolC. ', { italics: true }),
        t('eLife, 2022 (PDB 6WXI). '), link('elifesciences.org/articles/73297', 'https://elifesciences.org/articles/73297')], { after: 110 }),
      rich([t('Abdali N, et al. '), t('Self-inhibitory peptides targeting the Neisseria gonorrhoeae MtrCDE efflux '
        + 'pump increase antibiotic susceptibility. ', { italics: true }), t('Antimicrob Agents Chemother, 2022. '),
        link('PMC8765275', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8765275/')], { after: 110 }),
      rich([t('Wang S, et al. '), t('Gonococcal MtrE and its surface-expressed Loop 2 are immunogenic and elicit '
        + 'bactericidal antibodies. ', { italics: true }), t('Journal of Infection, 2018. '),
        link('PubMed 29902495', 'https://pubmed.ncbi.nlm.nih.gov/29902495/')], { after: 110 }),
      rich([t('MtrE Loop2-specific multiple antigenic peptide vaccine and monoclonal antibody confer '
        + 'complement-dependent protection against Neisseria gonorrhoeae. ', { italics: true }),
        t('npj Vaccines, 2026. '), link('nature.com/articles/s41541-026-01412-0', 'https://www.nature.com/articles/s41541-026-01412-0')], { after: 110 }),

      rule(),
      p('Figure 1a, the measured geometry and the target files are reproduced by '
        + 'scripts/prepare_mtre_targets.py and scripts/figure_targets.py in the mtre_binder_design directory.',
        { size: 18, italics: true, color: '7A8089' }),
    ],
  }],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(OUT, buf);
  console.log(`wrote ${OUT} (${(buf.length / 1024).toFixed(0)} kB)`);
});
