// Build the ANPDB exploratory novelty + multi-method target-prediction summary deck.
// Replaces the earlier T2D-only deck — this version reflects the full exploratory
// study (1,012 truly-novel compounds, PIDGIN + ChEMBL-NN consensus, IDF clustering,
// ophthalmology cluster, scaffold/ADME/pathway/UMAP/ethno layers).

const pptxgen = require("pptxgenjs");
const path = require("path");

const ROOT  = path.resolve(__dirname, "..");
const FIG   = path.join(ROOT, "output", "figures");
const SCAF  = path.join(ROOT, "output", "scaffolds");
const ADME  = path.join(ROOT, "output", "adme");
const PATH_ = path.join(ROOT, "output", "pathways");
const CHEM  = path.join(ROOT, "output", "chemspace");
const NETW  = path.join(ROOT, "output", "network", "figures");
const OPH   = path.join(ROOT, "output", "network", "ophthalmology", "figures");
const DOCK  = path.join(ROOT, "output", "docking");
const OUT   = path.join(ROOT, "output", "anpdb_exploratory_summary.pptx");

const C = {
  navy:    "0F1B3D",
  blue:    "065A82",
  teal:    "1C7293",
  amber:   "EE854A",
  cream:   "F4F7FB",
  white:   "FFFFFF",
  ink:     "21295C",
  muted:   "5C6B73",
  divider: "C9D6E2",
  green:   "2C8C5A",
  red:     "B33A3A",
};

const F = { head: "Cambria", body: "Calibri" };

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.author = "ANPDB exploratory target-prediction study";
pres.title  = "ANPDB Novelty + Multi-method Target Prediction";
const W = 13.3, H = 7.5;

const TOTAL_PAGES = 24;

function addFooter(slide, page) {
  slide.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: H - 0.32, w: W, h: 0.32, fill: { color: C.navy }, line: { color: C.navy },
  });
  slide.addText("ANPDB exploratory novelty + multi-method target prediction", {
    x: 0.5, y: H - 0.32, w: 8, h: 0.32, fontFace: F.body, fontSize: 9,
    color: C.cream, valign: "middle", margin: 0,
  });
  slide.addText(`${page} / ${TOTAL_PAGES}`, {
    x: W - 1.5, y: H - 0.32, w: 1, h: 0.32, fontFace: F.body, fontSize: 9,
    color: C.cream, valign: "middle", align: "right", margin: 0,
  });
}

function addContentTitle(slide, kicker, title) {
  slide.addShape(pres.shapes.RECTANGLE, {
    x: 0.5, y: 0.45, w: 0.12, h: 0.7, fill: { color: C.amber }, line: { color: C.amber },
  });
  slide.addText(kicker, {
    x: 0.75, y: 0.42, w: 11, h: 0.32, fontFace: F.body, fontSize: 11, bold: true,
    color: C.amber, charSpacing: 6, margin: 0,
  });
  slide.addText(title, {
    x: 0.75, y: 0.74, w: 12, h: 0.55, fontFace: F.head, fontSize: 26, bold: true,
    color: C.ink, margin: 0,
  });
}

// Card with image + commentary on the right
function addImageWithCommentary(slide, imgPath, kicker, kickerColor, commentary) {
  slide.addImage({
    path: imgPath,
    x: 0.6, y: 1.55, w: 8.7, h: 5.3, sizing: { type: "contain", w: 8.7, h: 5.3 },
  });
  slide.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.55, w: 3.2, h: 5.3, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
    shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
  });
  slide.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.55, w: 3.2, h: 0.45, fill: { color: kickerColor }, line: { color: kickerColor },
  });
  slide.addText(kicker, {
    x: 9.7, y: 1.6, w: 3.0, h: 0.32, fontFace: F.body, fontSize: 10, bold: true,
    color: C.white, charSpacing: 6, margin: 0,
  });
  slide.addText(commentary, {
    x: 9.7, y: 2.1, w: 3.0, h: 4.6, fontFace: F.body, fontSize: 11, color: C.ink,
    lineSpacingMultiple: 1.3, margin: 0, valign: "top",
  });
}

// ---------- Slide 1: Title ----------
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: H - 0.6, w: W, h: 0.6, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("NOVELTY  /  MULTI-METHOD TARGET PREDICTION  /  POLYPHARMACOLOGY", {
    x: 1, y: 1.2, w: 11.5, h: 0.4, fontFace: F.body, fontSize: 12, bold: true,
    color: C.amber, charSpacing: 8,
  });
  s.addText("Novel compounds in the African Natural Products Database\nand their predicted disease-relevant targets", {
    x: 1, y: 1.75, w: 11.3, h: 2.4, fontFace: F.head, fontSize: 34, bold: true,
    color: C.white, lineSpacingMultiple: 1.1,
  });
  s.addText([
    { text: "ANPDB v.2026  ·  COCONUT 05-2026  ·  ChEMBL 36  ·  PIDGINv4  ·  ChEMBL k-NN", options: { color: "B0C4DE" } },
  ], {
    x: 1, y: 4.3, w: 11, h: 0.5, fontFace: F.body, fontSize: 14,
  });
  s.addShape(pres.shapes.LINE, {
    x: 1, y: 5.0, w: 4, h: 0, line: { color: C.amber, width: 1.5 },
  });
  s.addText("Three actionable leads:  CYP1B1 / glaucoma  ·  SGLT2 / type 2 diabetes  ·  MAO-B / Parkinson's", {
    x: 1, y: 5.1, w: 11.5, h: 0.4, fontFace: F.body, fontSize: 13,
    color: C.cream, italic: true,
  });
  s.addText("Generated " + new Date().toISOString().slice(0,10), {
    x: 1, y: 5.55, w: 11, h: 0.4, fontFace: F.body, fontSize: 11,
    color: "9AAFCC", italic: true,
  });
}

// ---------- Slide 2: The question ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "THE QUESTION", "What can ANPDB tell us beyond what is already in COCONUT and ChEMBL?");

  s.addText([
    { text: "Gap in the literature", options: { bold: true, color: C.ink, fontSize: 16, breakLine: true } },
    { text: "ANPDB (11,448 compounds, 2026 release) is the largest African natural-product collection. No published study has applied modern machine-learning target prediction at scale to any African NP database — ANPDB's own paper relies on ePharmaLib pharmacophore screens. ", options: { color: C.ink } },
    { text: "We address this gap.", options: { bold: true, color: C.amber } },
  ], { x: 0.75, y: 1.7, w: 6.0, h: 4.5, fontFace: F.body, fontSize: 13, color: C.ink, lineSpacingMultiple: 1.3, valign: "top" });

  const cardX = 7.2, cardW = 5.4;
  const items = [
    { tag: "STAGE 1", text: "Define structurally novel ANPDB compounds vs COCONUT 2.0 + ChEMBL 36 (two-pass: InChIKey skeleton + Tanimoto T<0.85)." },
    { tag: "STAGE 2", text: "Predict targets with two orthogonal methods: PIDGINv4 random-forest classifiers + local ChEMBL k-NN on PIDGIN's bioactivity dataset." },
    { tag: "STAGE 3", text: "Triangulate consensus + IDF-weighted disease clustering + ethnopharmacological evidence to surface drug leads." },
  ];
  let y = 1.7;
  for (const it of items) {
    s.addShape(pres.shapes.RECTANGLE, {
      x: cardX, y, w: cardW, h: 1.4, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: cardX, y, w: 0.1, h: 1.4, fill: { color: C.teal }, line: { color: C.teal },
    });
    s.addText(it.tag, {
      x: cardX + 0.25, y: y + 0.15, w: cardW - 0.4, h: 0.3,
      fontFace: F.body, fontSize: 10, bold: true, color: C.teal, charSpacing: 6, margin: 0,
    });
    s.addText(it.text, {
      x: cardX + 0.25, y: y + 0.45, w: cardW - 0.4, h: 0.85,
      fontFace: F.body, fontSize: 12, color: C.ink, margin: 0, valign: "top",
    });
    y += 1.55;
  }
  addFooter(s, 2);
}

// ---------- Slide 3: Datasets ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "INPUTS", "Three databases anchor the novelty filter");
  const col = [
    { name: "ANPDB",     size: "11,448",     unit: "natural products",
      sub:  "Northern, Eastern + Southern Africa, U Freiburg curation (NAR 2026).",
      role: "QUERY SET" },
    { name: "COCONUT",   size: "738,827",    unit: "structures",
      sub:  "Open natural-product collection (May 2026), 479,721 unique InChIKey skeletons.",
      role: "NP REFERENCE" },
    { name: "ChEMBL 36", size: "2,715,471",  unit: "InChIKey skeletons",
      sub:  "Bioactive small-molecule reference (chembl_chemreps).",
      role: "BIOACTIVE REFERENCE" },
  ];
  const y0 = 1.85, h = 4.4;
  const colW = 3.95, gap = 0.32;
  const x0 = (W - 3*colW - 2*gap) / 2;
  col.forEach((c, i) => {
    const x = x0 + i * (colW + gap);
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: y0, w: colW, h, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 10, offset: 3, angle: 90, opacity: 0.08 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: y0, w: colW, h: 0.5, fill: { color: i===0 ? C.amber : C.teal }, line: { color: i===0 ? C.amber : C.teal },
    });
    s.addText(c.role, {
      x: x + 0.3, y: y0 + 0.05, w: colW - 0.6, h: 0.4,
      fontFace: F.body, fontSize: 10, bold: true, color: C.white, charSpacing: 6, margin: 0,
    });
    s.addText(c.name, {
      x: x + 0.3, y: y0 + 0.7, w: colW - 0.6, h: 0.6,
      fontFace: F.head, fontSize: 22, bold: true, color: C.ink, margin: 0,
    });
    s.addText(c.size, {
      x: x + 0.3, y: y0 + 1.55, w: colW - 0.6, h: 1.0,
      fontFace: F.head, fontSize: 36, bold: true, color: i===0 ? C.amber : C.blue, margin: 0,
    });
    s.addText(c.unit, {
      x: x + 0.3, y: y0 + 2.55, w: colW - 0.6, h: 0.4,
      fontFace: F.body, fontSize: 13, color: C.muted, margin: 0,
    });
    s.addText(c.sub, {
      x: x + 0.3, y: y0 + 3.05, w: colW - 0.6, h: 1.2,
      fontFace: F.body, fontSize: 12, color: C.ink, margin: 0, valign: "top",
    });
  });
  addFooter(s, 3);
}

// ---------- Slide 4: Pipeline ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "PIPELINE", "Standardise → filter → predict (×2) → consensus");

  const steps = [
    { n: "1", title: "Standardise + filter",
      detail: "Salt-strip + uncharge + normalise (RDKit MolStandardize). InChIKey skeleton match vs COCONUT∪ChEMBL, then Morgan-ECFP4 Tanimoto < 0.85 vs COCONUT.",
      out:    "1,012 truly novel" },
    { n: "2", title: "PIDGINv4",
      detail: "Random-forest classifiers on ECFP4 (PubChem + ChEMBL_28) at 100 nM, 1,718 targets. ad60 + ad90 reliability cutoffs.",
      out:    "949 cpds with hits" },
    { n: "3", title: "ChEMBL k-NN",
      detail: "Local k-NN target inference, K=20, sim_floor=0.30, support≥2, score≥0.40, against PIDGIN's bioactivity_dataset (263k actives).",
      out:    "344 predictions" },
    { n: "4", title: "Consensus + downstream",
      detail: "Intersect PIDGIN ∩ ChEMBL-NN. IDF-weighted disease clustering. Pathway enrichment, ADME, scaffolds, UMAP, ethnopharmacology.",
      out:    "83 consensus pairs" },
  ];
  const sx = 0.6, sy = 1.7, sw = 12.1, sh = 1.18;
  steps.forEach((st, i) => {
    const y = sy + i * (sh + 0.13);
    s.addShape(pres.shapes.RECTANGLE, {
      x: sx, y, w: sw, h: sh,
      fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.OVAL, {
      x: sx + 0.25, y: y + 0.22, w: 0.75, h: 0.75, fill: { color: C.navy }, line: { color: C.navy },
    });
    s.addText(st.n, {
      x: sx + 0.25, y: y + 0.22, w: 0.75, h: 0.75,
      fontFace: F.head, fontSize: 22, bold: true, color: C.amber,
      align: "center", valign: "middle", margin: 0,
    });
    s.addText(st.title, {
      x: sx + 1.2, y: y + 0.13, w: 4.6, h: 0.4,
      fontFace: F.head, fontSize: 17, bold: true, color: C.ink, margin: 0,
    });
    s.addText(st.detail, {
      x: sx + 1.2, y: y + 0.55, w: 7.6, h: 0.6,
      fontFace: F.body, fontSize: 11, color: C.muted, margin: 0, valign: "top",
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: sx + sw - 2.9, y: y + 0.18, w: 2.7, h: sh - 0.36, fill: { color: C.cream }, line: { color: C.divider, width: 0.5 },
    });
    s.addText("OUTPUT", {
      x: sx + sw - 2.9, y: y + 0.22, w: 2.7, h: 0.25,
      fontFace: F.body, fontSize: 8, bold: true, color: C.muted, charSpacing: 6, align: "center", margin: 0,
    });
    s.addText(st.out, {
      x: sx + sw - 2.9, y: y + 0.45, w: 2.7, h: 0.55,
      fontFace: F.head, fontSize: 16, bold: true, color: C.amber, align: "center", margin: 0,
    });
  });
  addFooter(s, 4);
}

// ---------- Slide 5: Novelty results ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "NOVELTY", "1,012 ANPDB compounds (8.8 %) survive both filters");

  // Big number callout
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 1.7, w: 5.5, h: 4.7, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 1.7, w: 5.5, h: 0.18, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("TRULY NOVEL ANPDB", {
    x: 0.95, y: 1.95, w: 5.1, h: 0.4,
    fontFace: F.body, fontSize: 11, bold: true, color: C.amber, charSpacing: 8, margin: 0,
  });
  s.addText("1,012", {
    x: 0.95, y: 2.45, w: 5.1, h: 1.6,
    fontFace: F.head, fontSize: 88, bold: true, color: C.white, margin: 0,
  });
  s.addText("of 11,448 ANPDB compounds", {
    x: 0.95, y: 4.05, w: 5.1, h: 0.4,
    fontFace: F.body, fontSize: 14, color: "B0C4DE", margin: 0,
  });
  s.addText("8.8 %", {
    x: 0.95, y: 4.5, w: 5.1, h: 0.7,
    fontFace: F.head, fontSize: 36, bold: true, color: C.amber, margin: 0,
  });
  s.addText("Standardised parent structure absent from both COCONUT (skeleton + Tanimoto<0.85) and ChEMBL (skeleton).", {
    x: 0.95, y: 5.3, w: 5.1, h: 0.9,
    fontFace: F.body, fontSize: 11, color: "CADCFC", italic: true, margin: 0, valign: "top",
  });

  const tx = 6.7, ty = 1.7;
  const rows = [
    [{ text: "Stage", options: { bold: true, color: C.white, fill: { color: C.blue } } },
     { text: "Compounds", options: { bold: true, color: C.white, fill: { color: C.blue }, align: "right" } },
     { text: "%", options: { bold: true, color: C.white, fill: { color: C.blue }, align: "right" } }],
    [ "ANPDB total",                       { text: "11,448", options: { align: "right" } }, { text: "100.0", options: { align: "right" } }],
    [ "After standardisation (parent)",    { text: "11,448", options: { align: "right" } }, { text: "100.0", options: { align: "right" } }],
    [ "Pass 1 — InChIKey skeleton novel",  { text: "1,364",  options: { align: "right" } }, { text: "11.9",  options: { align: "right" } }],
    [ "Pass 2 — Tanimoto < 0.85 vs COCONUT", { text: "1,012",  options: { align: "right" } }, { text: "8.8",   options: { align: "right" } }],
    [ { text: "  truly novel set", options: { bold: true, color: C.amber } },
      { text: "1,012", options: { align: "right", bold: true, color: C.amber } },
      { text: "8.8",   options: { align: "right", bold: true, color: C.amber } }],
  ];
  s.addTable(rows, {
    x: tx, y: ty, w: 6.0, colW: [3.4, 1.4, 1.2],
    rowH: 0.45,
    fontFace: F.body, fontSize: 12, color: C.ink,
    border: { type: "solid", pt: 0.5, color: C.divider },
    fill: { color: C.white },
  });
  s.addText("Notes", {
    x: tx, y: 5.1, w: 6.0, h: 0.35,
    fontFace: F.body, fontSize: 11, bold: true, color: C.muted, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "•  RDKit MolStandardize (LargestFragment + Normalize + Uncharge).", options: { breakLine: true } },
    { text: "•  Tanimoto on Morgan ECFP4 (radius 2, 2048 bits).", options: { breakLine: true } },
    { text: "•  COCONUT fingerprints cached locally (738,823 valid).", options: {} },
  ], {
    x: tx, y: 5.5, w: 6.0, h: 1.5,
    fontFace: F.body, fontSize: 11, color: C.ink, margin: 0, lineSpacingMultiple: 1.3,
  });
  addFooter(s, 5);
}

// ---------- Slide 6: Chemical space (UMAP) ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "CHEMICAL SPACE", "Novel ANPDB compounds occupy distinct regions of UMAP");
  addImageWithCommentary(s,
    path.join(CHEM, "fig_chemical_space_umap.png"),
    "READING", C.teal,
    [
      { text: "Each dot is a Morgan ECFP4 fingerprint, projected by UMAP (Jaccard, k=30).\n\n", options: {} },
      { text: "Grey: 5,000-cpd random COCONUT background.\n", options: { breakLine: true } },
      { text: "Red: 1,012 novel ANPDB compounds.\n\n", options: { color: C.red, bold: true, breakLine: true } },
      { text: "The novel set is not concentrated in one region — pockets of distinct chemistry exist throughout the natural-product space.", options: {} },
    ]);
  addFooter(s, 6);
}

// ---------- Slide 7: Scaffolds ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "SCAFFOLD DIVERSITY", "700 unique Bemis-Murcko scaffolds across 986 valid compounds");
  s.addImage({
    path: path.join(SCAF, "fig_scaffold_diversity.png"),
    x: 0.6, y: 1.6, w: 6.5, h: 4.0, sizing: { type: "contain", w: 6.5, h: 4.0 },
  });
  s.addImage({
    path: path.join(SCAF, "fig_scaffold_freq.png"),
    x: 7.3, y: 1.6, w: 5.5, h: 5.4, sizing: { type: "contain", w: 5.5, h: 5.4 },
  });

  // bottom callout strip
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.6, y: 5.85, w: 6.5, h: 1.2, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText([
    { text: "Diversity ratio  ", options: { color: C.amber, bold: true, fontSize: 12, charSpacing: 5 } },
    { text: "0.71 ", options: { color: C.white, fontSize: 18, bold: true, fontFace: F.head } },
    { text: "scaffolds/compound\n", options: { color: C.cream, fontSize: 12, breakLine: true } },
    { text: "Singletons        ", options: { color: C.amber, bold: true, fontSize: 12, charSpacing: 5 } },
    { text: "80 % ", options: { color: C.white, fontSize: 18, bold: true, fontFace: F.head } },
    { text: "of scaffolds appear in only one compound\n", options: { color: C.cream, fontSize: 12, breakLine: true } },
    { text: "Top scaffold      ", options: { color: C.amber, bold: true, fontSize: 12, charSpacing: 5 } },
    { text: "benzene (3.4 %), flavone (2.2 %), anthraquinone (1.0 %)", options: { color: C.cream, fontSize: 12 } },
  ], { x: 0.8, y: 5.95, w: 6.3, h: 1.05, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.3 });
  addFooter(s, 7);
}

// ---------- Slide 8: PIDGIN target prediction ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "TARGET PREDICTION (1/2)", "PIDGINv4 — random-forest classifiers across 1,718 targets");

  s.addText([
    { text: "Method", options: { bold: true, color: C.ink, fontSize: 15, breakLine: true } },
    { text: "Random-forest binary classifiers on ECFP4, trained on PubChem + ChEMBL_28 actives at 100 nM. Reliability-density applicability domain (AD) per prediction.\n\n", options: { color: C.ink } },
    { text: "We did not pre-filter for Lipinski/PAINS/ADMET. ", options: { bold: true, color: C.amber } },
    { text: "An earlier pre-filtered run excluded 92 % of the truly-novel set.", options: { color: C.ink } },
  ], { x: 0.75, y: 1.65, w: 6.0, h: 4.5, fontFace: F.body, fontSize: 12, color: C.ink, lineSpacingMultiple: 1.3, valign: "top" });

  const tiles = [
    { label: "Compounds scored", val: "1,012", note: "all truly-novel" },
    { label: "Compounds with ad60 hit", val: "949",   note: "94 % of input" },
    { label: "Distinct targets hit (ad60)", val: "1,378", note: "out of 1,718" },
  ];
  const tx0 = 7.0, ty0 = 1.7, tw = 1.85, th = 2.2, tg = 0.18;
  tiles.forEach((t, i) => {
    const x = tx0 + i * (tw + tg);
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: ty0, w: tw, h: th, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: ty0, w: tw, h: 0.1, fill: { color: C.teal }, line: { color: C.teal },
    });
    s.addText(t.label, {
      x: x + 0.15, y: ty0 + 0.2, w: tw - 0.3, h: 0.6,
      fontFace: F.body, fontSize: 9, bold: true, color: C.muted, charSpacing: 5, margin: 0,
    });
    s.addText(t.val, {
      x: x + 0.15, y: ty0 + 0.8, w: tw - 0.3, h: 0.9,
      fontFace: F.head, fontSize: 22, bold: true, color: C.blue, margin: 0,
    });
    s.addText(t.note, {
      x: x + 0.15, y: ty0 + 1.7, w: tw - 0.3, h: 0.45,
      fontFace: F.body, fontSize: 9, color: C.muted, margin: 0, italic: true, valign: "top",
    });
  });

  // bottom strip with caveat
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 4.2, w: 11.85, h: 2.7, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
    shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
  });
  s.addText("Promiscuity caveat", {
    x: 0.95, y: 4.35, w: 11.5, h: 0.4,
    fontFace: F.body, fontSize: 12, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  s.addText("PIDGIN's training set is biased towards human therapeutic targets. The raw hit-list is dominated by promiscuous CYPs (P05177/CYP1A2 alone returns ~ 40 cpds), and disease annotations are uneven (e.g. \"Schizophrenia\" appears on 121 hit compounds simply because many enzymes have a schizophrenia line in their disease record). We address this with (a) an orthogonal k-NN method on the same active universe, and (b) IDF re-weighting at the disease level — see Slides 9–11.",
    {
      x: 0.95, y: 4.7, w: 11.5, h: 2.1,
      fontFace: F.body, fontSize: 12, color: C.ink, margin: 0, valign: "top", lineSpacingMultiple: 1.3,
    });
  addFooter(s, 8);
}

// ---------- Slide 9: ChEMBL k-NN ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "TARGET PREDICTION (2/2)", "Local ChEMBL k-NN — orthogonal validator");

  s.addText([
    { text: "Method", options: { bold: true, color: C.ink, fontSize: 15, breakLine: true } },
    { text: "For each query: top-K nearest ChEMBL actives by Morgan-ECFP4 Tanimoto. Each neighbour contributes ", options: { color: C.ink } },
    { text: "score(target) += sim(q,n) ", options: { fontFace: "Consolas", color: C.ink } },
    { text: "for every UniProt it is annotated as ACTIVE against. Aggregate → normalise by Σ-similarity.\n\n", options: { color: C.ink } },
    { text: "Parameters", options: { bold: true, color: C.ink, breakLine: true } },
    { text: "K=20  ·  sim_floor=0.30  ·  min support=2  ·  min normalised score=0.40", options: { color: C.ink, fontFace: "Consolas" } },
  ], { x: 0.75, y: 1.65, w: 6.0, h: 4.5, fontFace: F.body, fontSize: 12, color: C.ink, lineSpacingMultiple: 1.3, valign: "top" });

  const tiles = [
    { label: "ChEMBL active universe",  val: "263,224", note: "compounds @ 100 nM" },
    { label: "Distinct targets",        val: "2,634",   note: "active in any cpd" },
    { label: "Predictions retained",    val: "344",     note: "293 cpds × 104 tgts" },
  ];
  const tx0 = 7.0, ty0 = 1.7, tw = 1.85, th = 2.2, tg = 0.18;
  tiles.forEach((t, i) => {
    const x = tx0 + i * (tw + tg);
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: ty0, w: tw, h: th, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: ty0, w: tw, h: 0.1, fill: { color: i===2 ? C.amber : C.teal }, line: { color: i===2 ? C.amber : C.teal },
    });
    s.addText(t.label, {
      x: x + 0.15, y: ty0 + 0.2, w: tw - 0.3, h: 0.6,
      fontFace: F.body, fontSize: 9, bold: true, color: C.muted, charSpacing: 5, margin: 0,
    });
    s.addText(t.val, {
      x: x + 0.15, y: ty0 + 0.8, w: tw - 0.3, h: 0.9,
      fontFace: F.head, fontSize: 20, bold: true, color: i===2 ? C.amber : C.blue, margin: 0,
    });
    s.addText(t.note, {
      x: x + 0.15, y: ty0 + 1.7, w: tw - 0.3, h: 0.45,
      fontFace: F.body, fontSize: 9, color: C.muted, margin: 0, italic: true, valign: "top",
    });
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 4.2, w: 11.85, h: 2.7, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText("Why this method matters", {
    x: 0.95, y: 4.32, w: 11.5, h: 0.4,
    fontFace: F.body, fontSize: 12, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "Model-free.", options: { bold: true, color: C.amber } },
    { text: " A target inference is just \"my query is similar to several known actives at this target.\" No classifier failure mode, no AD complications.\n\n", options: { color: C.cream } },
    { text: "Same universe.", options: { bold: true, color: C.amber } },
    { text: " Uses PIDGIN's own bioactivity_dataset, so an agreement reflects two independent signal-extraction methods on the same evidence — not a difference in training data.\n\n", options: { color: C.cream } },
    { text: "Recovered SGLT2.", options: { bold: true, color: C.amber } },
    { text: " The canonical T2D target was missing from PIDGIN's annotation-driven T2D analysis but emerged through ChEMBL-NN consensus — see Slide 11.", options: { color: C.cream } },
  ], { x: 0.95, y: 4.7, w: 11.5, h: 2.1, fontFace: F.body, fontSize: 11, valign: "top", margin: 0, lineSpacingMultiple: 1.3 });
  addFooter(s, 9);
}

// ---------- Slide 10: Consensus ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "MULTI-METHOD CONSENSUS", "PIDGIN ∩ ChEMBL-NN — 83 high-confidence pairs");

  // Venn-style stats area (using two large overlapping cards)
  const cy = 1.7, ch = 2.4;
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: cy, w: 5.7, h: ch, fill: { color: C.blue }, line: { color: C.blue },
  });
  s.addText("PIDGINv4 (ad60)", {
    x: 0.95, y: cy + 0.18, w: 5.3, h: 0.4,
    fontFace: F.body, fontSize: 12, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  s.addText("16,427", {
    x: 0.95, y: cy + 0.5, w: 5.3, h: 0.9,
    fontFace: F.head, fontSize: 38, bold: true, color: C.white, margin: 0,
  });
  s.addText("(compound, target) pairs    ·    949 cpds    ·    1,378 targets",
    { x: 0.95, y: cy + 1.5, w: 5.3, h: 0.7, fontFace: F.body, fontSize: 11, color: "B0C4DE", margin: 0, valign: "top" });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 6.85, y: cy, w: 5.7, h: ch, fill: { color: C.teal }, line: { color: C.teal },
  });
  s.addText("ChEMBL k-NN", {
    x: 7.05, y: cy + 0.18, w: 5.3, h: 0.4,
    fontFace: F.body, fontSize: 12, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  s.addText("344", {
    x: 7.05, y: cy + 0.5, w: 5.3, h: 0.9,
    fontFace: F.head, fontSize: 38, bold: true, color: C.white, margin: 0,
  });
  s.addText("predictions    ·    293 cpds    ·    104 targets",
    { x: 7.05, y: cy + 1.5, w: 5.3, h: 0.7, fontFace: F.body, fontSize: 11, color: "B0C4DE", margin: 0, valign: "top" });

  // Big consensus block
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 4.4, w: 11.85, h: 2.5, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 4.4, w: 11.85, h: 0.18, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText("CONSENSUS  =  PIDGIN ∩ ChEMBL-NN", {
    x: 0.95, y: 4.65, w: 11.5, h: 0.4,
    fontFace: F.body, fontSize: 12, bold: true, color: C.navy, charSpacing: 8, margin: 0,
  });
  s.addText("83 pairs", {
    x: 0.95, y: 5.0, w: 5.5, h: 1.0,
    fontFace: F.head, fontSize: 44, bold: true, color: C.navy, margin: 0,
  });
  s.addText("80 compounds  ×  21 targets", {
    x: 0.95, y: 5.95, w: 6.0, h: 0.45,
    fontFace: F.body, fontSize: 14, bold: true, color: C.navy, margin: 0,
  });
  s.addText([
    { text: "Top consensus targets\n", options: { bold: true, color: C.navy, fontSize: 12, charSpacing: 4, breakLine: true } },
    { text: "CYP1B1 (Q16678) — 20 cpds — glaucoma\n", options: { color: C.navy, fontSize: 11, breakLine: true } },
    { text: "FimH-like (P08191) — 19 cpds\n", options: { color: C.navy, fontSize: 11, breakLine: true } },
    { text: "SGLT2 (P31639) — 14 cpds — type 2 diabetes\n", options: { color: C.navy, fontSize: 11, breakLine: true } },
    { text: "HSD11B1 (P80365) — 5 cpds — metabolic\n", options: { color: C.navy, fontSize: 11, breakLine: true } },
    { text: "MAO-B (P27338) — 1 cpd — Parkinson's", options: { color: C.navy, fontSize: 11 } },
  ], { x: 7.5, y: 4.8, w: 5.0, h: 2.0, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.2 });
  addFooter(s, 10);
}

// ---------- Slide 11: SGLT2 — the win ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "WHY MULTI-METHOD MATTERS", "SGLT2: invisible to PIDGIN, recovered through consensus");

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 1.7, w: 11.85, h: 5.3, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
    shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
  });

  s.addText("The vignette", {
    x: 0.95, y: 1.85, w: 11.5, h: 0.4,
    fontFace: F.body, fontSize: 11, bold: true, color: C.amber, charSpacing: 8, margin: 0,
  });
  s.addText("SGLT2 (P31639) is the canonical type-2-diabetes drug target. Dapagliflozin, empagliflozin, canagliflozin all target it.",
    { x: 0.95, y: 2.2, w: 11.5, h: 0.6, fontFace: F.body, fontSize: 14, color: C.ink, margin: 0, valign: "top", italic: true });

  // Two columns: PIDGIN view vs Consensus view
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.95, y: 3.0, w: 5.7, h: 3.7, fill: { color: C.cream }, line: { color: C.divider, width: 0.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.95, y: 3.0, w: 5.7, h: 0.4, fill: { color: C.red }, line: { color: C.red },
  });
  s.addText("PIDGIN-only T2D analysis", {
    x: 1.1, y: 3.05, w: 5.4, h: 0.32,
    fontFace: F.body, fontSize: 11, bold: true, color: C.white, charSpacing: 5, margin: 0,
  });
  s.addText([
    { text: "SGLT2 not flagged.\n\n", options: { bold: true, color: C.red, fontSize: 14, breakLine: true } },
    { text: "PIDGIN's per-target disease lookup did not carry the canonical \"Diabetes Mellitus, Non-Insulin-Dependent\" string for P31639. The strict-NIDDM strength filter (≥4) therefore excluded SGLT2 from the T2D leaderboard.\n\n", options: { color: C.ink, fontSize: 11 } },
    { text: "0 ", options: { fontSize: 24, bold: true, color: C.red, fontFace: F.head } },
    { text: "novel ANPDB compounds reach an annotated T2D target via PIDGIN alone after strength filtering.", options: { color: C.ink, fontSize: 11 } },
  ], { x: 1.1, y: 3.5, w: 5.4, h: 3.1, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.3 });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 6.85, y: 3.0, w: 5.7, h: 3.7, fill: { color: C.cream }, line: { color: C.divider, width: 0.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 6.85, y: 3.0, w: 5.7, h: 0.4, fill: { color: C.green }, line: { color: C.green },
  });
  s.addText("PIDGIN ∩ ChEMBL-NN consensus", {
    x: 7.0, y: 3.05, w: 5.4, h: 0.32,
    fontFace: F.body, fontSize: 11, bold: true, color: C.white, charSpacing: 5, margin: 0,
  });
  s.addText([
    { text: "SGLT2 emerges.\n\n", options: { bold: true, color: C.green, fontSize: 14, breakLine: true } },
    { text: "ChEMBL-NN finds 14 ANPDB-novel compounds whose nearest active neighbours include known SGLT2 inhibitors (T > 0.45). PIDGIN agrees on the same 14 at ad60. The consensus does not depend on the disease lookup string.\n\n", options: { color: C.ink, fontSize: 11 } },
    { text: "14 ", options: { fontSize: 24, bold: true, color: C.green, fontFace: F.head } },
    { text: "novel compounds with cross-method support for SGLT2 — direct T2D leads.", options: { color: C.ink, fontSize: 11 } },
  ], { x: 7.0, y: 3.5, w: 5.4, h: 3.1, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.3 });
  addFooter(s, 11);
}

// ---------- Slide 12: Docking validation ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "STRUCTURAL VALIDATION", "AutoDock Vina docks 34 consensus leads with strong affinity");

  // Left: boxplot
  s.addImage({
    path: path.join(DOCK, "fig_dock_boxplot.png"),
    x: 0.55, y: 1.55, w: 7.6, h: 5.4, sizing: { type: "contain", w: 7.6, h: 5.4 },
  });

  // Right: per-target stat cards
  const cards = [
    { tname: "CYP1B1", colour: C.teal,  pdb: "4I8V", n: 19, best: -13.5, top: "Mol_11315" },
    { tname: "SGLT2",  colour: C.amber, pdb: "7VSI", n: 14, best: -11.8, top: "Mol_13144" },
    { tname: "MAO-B",  colour: C.green, pdb: "2V5Z", n: 1,  best: -10.3, top: "Mol_14056" },
  ];
  let cy = 1.6;
  cards.forEach((c) => {
    s.addShape(pres.shapes.RECTANGLE, {
      x: 8.4, y: cy, w: 4.4, h: 1.55, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: 8.4, y: cy, w: 0.13, h: 1.55, fill: { color: c.colour }, line: { color: c.colour },
    });
    s.addText(c.tname, {
      x: 8.6, y: cy + 0.1, w: 4.0, h: 0.35,
      fontFace: F.head, fontSize: 17, bold: true, color: c.colour, margin: 0,
    });
    s.addText(`PDB ${c.pdb}  ·  n=${c.n}`, {
      x: 8.6, y: cy + 0.45, w: 4.0, h: 0.3,
      fontFace: F.body, fontSize: 10, color: C.muted, margin: 0,
    });
    s.addText(`${c.best.toFixed(1)} kcal/mol`, {
      x: 8.6, y: cy + 0.75, w: 4.0, h: 0.5,
      fontFace: F.head, fontSize: 22, bold: true, color: C.ink, margin: 0,
    });
    s.addText(`top hit: ${c.top}`, {
      x: 8.6, y: cy + 1.22, w: 4.0, h: 0.3,
      fontFace: F.body, fontSize: 11, italic: true, color: C.muted, margin: 0,
    });
    cy += 1.7;
  });

  // bottom strip
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.55, y: 6.4, w: 12.25, h: 0.6, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText([
    { text: "All 34 valid poses ≤ −7 kcal/mol  ·  ", options: { color: C.amber, bold: true, fontSize: 12 } },
    { text: "82 % ≤ −9 kcal/mol  ·  ", options: { color: C.cream, fontSize: 12 } },
    { text: "Vina exhaustiveness=16, 9 modes/ligand", options: { color: "B0C4DE", fontSize: 11, italic: true } },
  ], { x: 0.75, y: 6.45, w: 12.0, h: 0.55, fontFace: F.body, valign: "middle", margin: 0 });
  addFooter(s, 12);
}

// ---------- Slide 13: Docking validation controls ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "VALIDATION CONTROLS", "All 3 receptors recover crystal pose to < 2 Å (top-9 modes)");

  // Three-panel validation figure (decoy + redock + CYP1B1 mode scatter)
  s.addImage({
    path: path.join(DOCK, "validation", "fig_dock_validation_v2.png"),
    x: 0.4, y: 1.45, w: 8.7, h: 5.7, sizing: { type: "contain", w: 8.7, h: 5.7 },
  });

  // Right: cards
  // Card A: redock best-of-9
  let cy = 1.5;
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.25, y: cy, w: 3.7, h: 2.5, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
    shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.25, y: cy, w: 3.7, h: 0.4, fill: { color: C.green }, line: { color: C.green },
  });
  s.addText("REDOCK  ·  best-of-9 RMSD", {
    x: 9.4, y: cy + 0.05, w: 3.5, h: 0.32,
    fontFace: F.body, fontSize: 9, bold: true, color: C.white, charSpacing: 5, margin: 0,
  });
  s.addText([
    { text: "CYP1B1   ", options: { color: C.ink, fontSize: 10, fontFace: "Consolas" } },
    { text: "0.39 Å ✓\n", options: { bold: true, color: C.green, fontSize: 11, fontFace: "Consolas", breakLine: true } },
    { text: "SGLT2    ", options: { color: C.ink, fontSize: 10, fontFace: "Consolas" } },
    { text: "1.28 Å ✓\n", options: { bold: true, color: C.green, fontSize: 11, fontFace: "Consolas", breakLine: true } },
    { text: "MAO-B    ", options: { color: C.ink, fontSize: 10, fontFace: "Consolas" } },
    { text: "0.92 Å ✓\n\n", options: { bold: true, color: C.green, fontSize: 11, fontFace: "Consolas", breakLine: true } },
    { text: "Mode 1 RMSD\n", options: { color: C.muted, fontSize: 9, breakLine: true } },
    { text: "CYP1B1: 3.84 Å (flipped, scoring artefact within 1 kcal/mol). Mode 1 sub-2 Å for SGLT2 (1.65) + MAO-B (1.61).",
      options: { color: C.muted, fontSize: 9, italic: true } },
  ], { x: 9.4, y: cy + 0.45, w: 3.5, h: 2.05, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.15 });

  // Card B: decoy
  cy += 2.65;
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.25, y: cy, w: 3.7, h: 2.5, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
    shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.25, y: cy, w: 3.7, h: 0.4, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("DECOY  ·  Mann–Whitney U", {
    x: 9.4, y: cy + 0.05, w: 3.5, h: 0.32,
    fontFace: F.body, fontSize: 9, bold: true, color: C.white, charSpacing: 5, margin: 0,
  });
  s.addText([
    { text: "Consensus < decoys (stronger binding):\n\n",
      options: { bold: true, color: C.amber, fontSize: 10, breakLine: true } },
    { text: "CYP1B1  ", options: { color: C.ink, fontSize: 10, fontFace: "Consolas" } },
    { text: "p=1.8×10⁻⁶ ✓\n", options: { bold: true, color: C.green, fontSize: 10, fontFace: "Consolas", breakLine: true } },
    { text: "SGLT2   ", options: { color: C.ink, fontSize: 10, fontFace: "Consolas" } },
    { text: "p=5.2×10⁻⁴ ✓\n", options: { bold: true, color: C.green, fontSize: 10, fontFace: "Consolas", breakLine: true } },
    { text: "MAO-B   ", options: { color: C.ink, fontSize: 10, fontFace: "Consolas" } },
    { text: "p=0.07 (n=1)\n\n", options: { color: C.muted, fontSize: 10, fontFace: "Consolas", breakLine: true } },
    { text: "n_consensus=34  ·  n_decoy=78  ·  decoys = no PIDGIN nor ChEMBL-NN hit.",
      options: { color: C.muted, fontSize: 9, italic: true } },
  ], { x: 9.4, y: cy + 0.45, w: 3.5, h: 2.05, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.15 });

  addFooter(s, 13);
}

// ---------- Slide 14: IDF disease leaderboard ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "DISEASE PROFILE", "IDF re-weighting deflates promiscuous terms");
  addImageWithCommentary(s,
    path.join(NETW, "fig_disease_leaderboard_ad60_idf.png"),
    "READING", C.teal,
    [
      { text: "IDF(d) = log( N_targets / df(d) )\n\n", options: { fontFace: "Consolas", fontSize: 10 } },
      { text: "N=1,718 targets, 2,178 disease terms.\n\n", options: { breakLine: false } },
      { text: "Promiscuous terms (Schizophrenia, Autosomal recessive predisposition) shrink — the leaderboard is now dominated by glaucoma + congenital eye disease, leukaemia, hearing loss, and a small set of specific neurological / metabolic terms.", options: {} },
    ]);
  addFooter(s, 14);
}

// ---------- Slide 15: IDF cluster summary ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "DISEASE CLUSTERING", "Six IDF-weighted clusters; clusters 4 + 6 = ophthalmology");
  addImageWithCommentary(s,
    path.join(NETW, "fig_cluster_summary_ad60_idf.png"),
    "READING", C.teal,
    [
      { text: "Hierarchical clustering, k=6, cosine + average linkage on the IDF-weighted compound × disease matrix (top 30 diseases).\n\n", options: {} },
      { text: "Clusters 4 + 6 ", options: { bold: true, color: C.amber } },
      { text: "are dominated by ophthalmology terms — 98 compounds (33 % of the predicted-hit set).\n\n", options: { breakLine: false } },
      { text: "Other clusters split along leukaemia, neurological, metabolic axes.", options: {} },
    ]);
  addFooter(s, 15);
}

// ---------- Slide 16: Ophthalmology — buckets ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "OPHTHALMOLOGY CLUSTER", "98 compounds spanning glaucoma, cataract, hearing & cornea");
  addImageWithCommentary(s,
    path.join(OPH, "fig_disease_buckets.png"),
    "OBSERVATION", C.amber,
    [
      { text: "Six near-synonym buckets:\n", options: { bold: true, breakLine: true } },
      { text: "•  Glaucoma — 59 cpds\n", options: { breakLine: true } },
      { text: "•  Hearing / Otology — 50\n", options: { breakLine: true } },
      { text: "•  Cataract / Lens — 48\n", options: { breakLine: true } },
      { text: "•  Cornea / Sclera — 38\n", options: { breakLine: true } },
      { text: "•  Retina / Optic — 33\n\n", options: { breakLine: true } },
      { text: "Major compounds reach ≥3 buckets — true polypharmacology.", options: { italic: true } },
    ]);
  addFooter(s, 16);
}

// ---------- Slide 17: Ophthalmology — top targets ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "OPHTHALMOLOGY TARGETS", "CYP1B1 + CA4 + FA2H drive the eye-disease signal");
  addImageWithCommentary(s,
    path.join(OPH, "fig_top_targets.png"),
    "TARGETS", C.teal,
    [
      { text: "CYP1B1 (Q16678) — 35 cpds\n", options: { bold: true, breakLine: true } },
      { text: "Mutations cause primary congenital glaucoma (PCG). Validated drug target.\n\n", options: { fontSize: 10, breakLine: true } },
      { text: "CA4 (P22748) — 28 cpds\n", options: { bold: true, breakLine: true } },
      { text: "Carbonic anhydrase IV. Existing class: dorzolamide (glaucoma eye drops).\n\n", options: { fontSize: 10, breakLine: true } },
      { text: "FA2H (P55789) — 20 cpds\n", options: { bold: true, breakLine: true } },
      { text: "Fatty-acid 2-hydroxylase. Mutations → spastic paraplegia + leukodystrophy with ophthalmic features.", options: { fontSize: 10 } },
    ]);
  addFooter(s, 17);
}

// ---------- Slide 18: Source species ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "SOURCE SPECIES", "Top species cluster around traditional eye-disease use");
  addImageWithCommentary(s,
    path.join(OPH, "fig_top_species.png"),
    "BIOLOGICAL ORIGIN", C.amber,
    [
      { text: "Top contributors:\n", options: { bold: true, breakLine: true } },
      { text: "Tephrosia purpurea\n", options: { italic: true, breakLine: true } },
      { text: "Solanum nigrum\n", options: { italic: true, breakLine: true } },
      { text: "Erythrina abyssinica\n", options: { italic: true, breakLine: true } },
      { text: "Combretum molle\n", options: { italic: true, breakLine: true } },
      { text: "Croton macrostachyus\n\n", options: { italic: true, breakLine: true } },
      { text: "Slide 19 cross-checks these against documented African eye-medicine ethnobotany.", options: {} },
    ]);
  addFooter(s, 18);
}

// ---------- Slide 19: Ethnopharmacology table ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "ETHNOPHARMACOLOGICAL VALIDATION", "4/10 leading species have documented traditional eye-disease use");

  const headerStyle = { bold: true, color: C.white, fill: { color: C.navy }, align: "left" };
  const rows = [
    [
      { text: "Species", options: headerStyle },
      { text: "Traditional eye-related use", options: headerStyle },
      { text: "Compound classes", options: headerStyle },
      { text: "Ref", options: headerStyle },
    ],
    [
      { text: "Tephrosia purpurea", options: { italic: true, color: C.amber, bold: true } },
      "Leaf juice — eye inflammation, night blindness",
      "Rotenoids, flavanones",
      "Soni 2013",
    ],
    [
      { text: "Solanum nigrum", options: { italic: true, color: C.amber, bold: true } },
      "Leaf decoction — conjunctivitis, eye wash (sub-Saharan)",
      "Steroidal glycoalkaloids, polyphenols",
      "Jain 2011",
    ],
    [
      { text: "Securidaca longipedunculata", options: { italic: true, color: C.amber, bold: true } },
      "Root-bark wash — eye infections (Mozambique, Tanzania)",
      "Xanthones, methyl salicylate, saponins",
      "Neuwinger 2000",
    ],
    [
      { text: "Croton macrostachyus", options: { italic: true, color: C.amber, bold: true } },
      "Leaf/latex — trachoma, conjunctivitis (Ethiopian TM)",
      "Clerodane diterpenoids, crotepoxide",
      "Mesfin 2009",
    ],
    [
      { text: "Warburgia ugandensis", options: { italic: true } },
      "Leaf steam — eye irritation (Kenyan, minor)",
      "Drimane sesquiterpenoids",
      "Maroyi 2014",
    ],
    [
      { text: "Erythrina abyssinica", options: { italic: true, color: C.muted } },
      { text: "no specific eye use documented", options: { color: C.muted, italic: true } },
      { text: "Pterocarpans, isoflavonoids", options: { color: C.muted } },
      { text: "Yenesew 2004", options: { color: C.muted } },
    ],
    [
      { text: "Combretum molle", options: { italic: true, color: C.muted } },
      { text: "no specific eye use documented", options: { color: C.muted, italic: true } },
      { text: "Tannins, stilbenoids", options: { color: C.muted } },
      { text: "Eloff 2008", options: { color: C.muted } },
    ],
  ];
  s.addTable(rows, {
    x: 0.55, y: 1.55, w: 12.2, colW: [2.3, 4.4, 3.4, 2.1],
    rowH: 0.55,
    fontFace: F.body, fontSize: 11, color: C.ink,
    border: { type: "solid", pt: 0.5, color: C.divider },
    fill: { color: C.white },
    valign: "middle",
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 6.1, w: 11.85, h: 0.85, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText([
    { text: "Take-away   ", options: { color: C.amber, bold: true, fontSize: 12, charSpacing: 5 } },
    { text: "An unsupervised IDF cluster (built only on PIDGIN target predictions) recovers a group of compounds whose source species are independently documented to treat eye disease. ", options: { color: C.cream, fontSize: 12 } },
    { text: "External validation, no annotation peeking.", options: { color: C.amber, italic: true, bold: true, fontSize: 12 } },
  ], { x: 0.95, y: 6.2, w: 11.5, h: 0.7, fontFace: F.body, valign: "middle", margin: 0, lineSpacingMultiple: 1.3 });
  addFooter(s, 19);
}

// ---------- Slide 20: Pathway enrichment ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "PATHWAY ENRICHMENT", "Xenobiotic + steroid metabolism dominate the consensus targets");
  addImageWithCommentary(s,
    path.join(PATH_, "fig_pathway_enrichment.png"),
    "g:Profiler", C.teal,
    [
      { text: "21 consensus UniProt → g:Profiler (g:SCS, p<0.05)\n\n", options: { fontSize: 10 } },
      { text: "11 enriched terms\n", options: { bold: true, breakLine: true } },
      { text: "•  KEGG Tryptophan metabolism (p=4.2e-4)\n", options: { fontSize: 10, breakLine: true } },
      { text: "•  KEGG Steroid hormone biosynthesis (p=1.9e-3)\n", options: { fontSize: 10, breakLine: true } },
      { text: "•  Reactome Phase I — Functionalisation\n", options: { fontSize: 10, breakLine: true } },
      { text: "•  GO Response to xenobiotic (5/396)\n\n", options: { fontSize: 10, breakLine: true } },
      { text: "Coherent with the CYP/MAO/HSD11B1 backbone.", options: { italic: true, fontSize: 10 } },
    ]);
  addFooter(s, 20);
}

// ---------- Slide 21: ADME / BOILED-Egg ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "ADME — TOP-50 SHORTLIST", "BOILED-Egg + Lipinski/Veber/Egan compliance");
  s.addImage({
    path: path.join(ADME, "fig_boiled_egg.png"),
    x: 0.6, y: 1.55, w: 6.5, h: 5.4, sizing: { type: "contain", w: 6.5, h: 5.4 },
  });
  s.addImage({
    path: path.join(ADME, "fig_adme_rules.png"),
    x: 7.3, y: 1.85, w: 5.6, h: 3.1, sizing: { type: "contain", w: 5.6, h: 3.1 },
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 7.3, y: 5.1, w: 5.6, h: 1.85, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText([
    { text: "Headline\n", options: { color: C.amber, bold: true, fontSize: 11, charSpacing: 6, breakLine: true } },
    { text: "53 % Lipinski-compliant (≤1 violation)\n", options: { color: C.cream, fontSize: 11, breakLine: true } },
    { text: "29 % predicted GI-absorbed\n", options: { color: C.cream, fontSize: 11, breakLine: true } },
    { text: "4 % predicted BBB-permeant\n\n", options: { color: C.cream, fontSize: 11, breakLine: true } },
    { text: "Median MW=404 · LogP=1.2 · TPSA=140\n", options: { color: C.cream, fontSize: 10, italic: true, breakLine: true } },
    { text: "Typical NP profile — Ro5 filtering would have removed many top consensus leads.", options: { color: C.amber, fontSize: 10, italic: true } },
  ], { x: 7.5, y: 5.2, w: 5.4, h: 1.7, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.3 });
  addFooter(s, 21);
}

// ---------- Slide 22: Three actionable leads ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "ACTIONABLE LEADS", "Three target hypotheses that survive triangulation");

  const leads = [
    { color: C.amber,
      tag: "GLAUCOMA",
      target: "CYP1B1 (Q16678)",
      n: "20",
      caption: "Twenty consensus compounds, predominantly from T. purpurea and S. nigrum. CYP1B1 mutations cause primary congenital glaucoma; the ophthalmology cluster + traditional eye-disease use of source species converge on the same hypothesis." },
    { color: C.green,
      tag: "TYPE 2 DIABETES",
      target: "SGLT2 (P31639)",
      n: "14",
      caption: "Fourteen consensus compounds — gliflozin-class target. Recovered only because ChEMBL-NN bypassed PIDGIN's incomplete disease annotation. Direct path to T2D follow-up assays." },
    { color: C.teal,
      tag: "PARKINSON'S",
      target: "MAO-B (P27338)",
      n: "1",
      caption: "Single high-confidence consensus hit. Top-neighbour Tanimoto > 0.6 against a known MAO-B-active ChEMBL ligand. Single-target lead, ready for orthogonal docking + biochemical screen." },
  ];

  let y = 1.7;
  leads.forEach((l, i) => {
    s.addShape(pres.shapes.RECTANGLE, {
      x: 0.75, y, w: 11.85, h: 1.65, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: 0.75, y, w: 0.16, h: 1.65, fill: { color: l.color }, line: { color: l.color },
    });
    s.addText(l.tag, {
      x: 1.05, y: y + 0.15, w: 3.3, h: 0.32,
      fontFace: F.body, fontSize: 10, bold: true, color: l.color, charSpacing: 8, margin: 0,
    });
    s.addText(l.target, {
      x: 1.05, y: y + 0.45, w: 3.3, h: 0.5,
      fontFace: F.head, fontSize: 19, bold: true, color: C.ink, margin: 0,
    });
    s.addText(l.n, {
      x: 1.05, y: y + 0.95, w: 1.4, h: 0.65,
      fontFace: F.head, fontSize: 36, bold: true, color: l.color, margin: 0,
    });
    s.addText("compounds", {
      x: 2.4, y: y + 1.18, w: 1.95, h: 0.4,
      fontFace: F.body, fontSize: 11, color: C.muted, italic: true, margin: 0,
    });
    s.addText(l.caption, {
      x: 4.5, y: y + 0.18, w: 8.0, h: 1.3,
      fontFace: F.body, fontSize: 12, color: C.ink, margin: 0, valign: "top", lineSpacingMultiple: 1.3,
    });
    y += 1.78;
  });
  addFooter(s, 22);
}

// ---------- Slide 23: Differentiation vs published work ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "DIFFERENTIATION", "Where this work sits in the African NP literature");

  const headerStyle = { bold: true, color: C.white, fill: { color: C.navy }, align: "left" };
  const rows = [
    [
      { text: "Database / Paper", options: headerStyle },
      { text: "Year", options: headerStyle },
      { text: "What was done", options: headerStyle },
      { text: "ML target prediction?", options: headerStyle },
    ],
    [ "AfroDb (Ntie-Kang)",        "2013", "Druglikeness + chemical-space PCA", { text: "no", options: { color: C.red, bold: true } } ],
    [ "SANCDB v1 (Hatherley)",     "2015", "DB construction + scaffolds",       { text: "no", options: { color: C.red, bold: true } } ],
    [ "NANPDB (Onguéné)",          "2017", "Physicochemical descriptors",        { text: "no", options: { color: C.red, bold: true } } ],
    [ "EANPDB review (Simoben)",   "2020", "Comparative DB review",              { text: "n/a", options: { color: C.muted } } ],
    [ "SANCDB v2 (Hatherley)",     "2021", "DB update + ADMET",                 { text: "no", options: { color: C.red, bold: true } } ],
    [ "ANPDB NAR paper",           "2026", "DB release + ePharmaLib pharmacophore",  { text: "structure-only (not ML)", options: { color: C.red, italic: true } } ],
    [ "Value Addition review",     "2026", "Survey — explicitly flags AI/ML as gap", { text: "n/a", options: { color: C.muted } } ],
    [
      { text: "This work", options: { bold: true, color: C.amber } },
      { text: "2026", options: { bold: true, color: C.amber } },
      { text: "PIDGINv4 + ChEMBL-NN consensus + IDF clustering + ethno", options: { bold: true, color: C.amber } },
      { text: "yes (first)", options: { bold: true, color: C.green } },
    ],
  ];
  s.addTable(rows, {
    x: 0.55, y: 1.55, w: 12.2, colW: [3.7, 0.8, 5.4, 2.3],
    rowH: 0.45,
    fontFace: F.body, fontSize: 11, color: C.ink,
    border: { type: "solid", pt: 0.5, color: C.divider },
    fill: { color: C.white },
    valign: "middle",
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 6.0, w: 11.85, h: 0.95, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText([
    { text: "First   ", options: { color: C.amber, bold: true, fontSize: 12, charSpacing: 5 } },
    { text: "systematic ML target prediction on any African NP database. ", options: { color: C.cream, fontSize: 12 } },
    { text: "First   ", options: { color: C.amber, bold: true, fontSize: 12, charSpacing: 5 } },
    { text: "multi-method consensus on African NPs. ", options: { color: C.cream, fontSize: 12 } },
    { text: "First   ", options: { color: C.amber, bold: true, fontSize: 12, charSpacing: 5 } },
    { text: "two-pass novelty filter against COCONUT 2.0 + ChEMBL 36 on ANPDB.", options: { color: C.cream, fontSize: 12 } },
  ], { x: 0.95, y: 6.1, w: 11.5, h: 0.8, fontFace: F.body, valign: "middle", margin: 0, lineSpacingMultiple: 1.3 });
  addFooter(s, 23);
}

// ---------- Slide 24: Summary & next ----------
{
  const s = pres.addSlide();
  s.background = { color: C.navy };

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 0.65, w: 0.12, h: 0.7, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("SUMMARY", {
    x: 1.0, y: 0.62, w: 11, h: 0.32, fontFace: F.body, fontSize: 11, bold: true,
    color: C.amber, charSpacing: 6, margin: 0,
  });
  s.addText("Where this leaves us", {
    x: 1.0, y: 0.94, w: 12, h: 0.55, fontFace: F.head, fontSize: 28, bold: true, color: C.white, margin: 0,
  });

  const lx = 0.75, ly = 1.85, lw = 5.9;
  const rx = 7.0, ry = 1.85, rw = 5.6;

  s.addText("Headline numbers", {
    x: lx, y: ly, w: lw, h: 0.4, fontFace: F.body, fontSize: 13, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  const lines = [
    ["1,012",  "truly novel ANPDB compounds (post-standardisation, T<0.85 vs COCONUT)"],
    ["83",     "consensus (compound, target) pairs — PIDGIN ∩ ChEMBL-NN"],
    ["98",     "compounds in the IDF-weighted ophthalmology cluster"],
    ["−13.5",  "best Vina affinity (Mol_11315 / CYP1B1) — 82 % of poses ≤ −9 kcal/mol"],
    ["1.8e−6", "Mann–Whitney p — consensus binds stronger than property-matched decoys"],
    ["20 / 14 / 1", "consensus compounds for CYP1B1 / SGLT2 / MAO-B"],
  ];
  let y = ly + 0.45;
  lines.forEach(([num, label]) => {
    s.addText(num, {
      x: lx, y, w: 2.2, h: 0.55, fontFace: F.head, fontSize: 24, bold: true, color: C.amber, margin: 0, align: "left",
    });
    s.addText(label, {
      x: lx + 2.25, y: y + 0.08, w: lw - 2.25, h: 0.55,
      fontFace: F.body, fontSize: 11, color: C.cream, margin: 0, valign: "top",
    });
    y += 0.65;
  });

  s.addText("Recommended next steps", {
    x: rx, y: ry, w: rw, h: 0.4, fontFace: F.body, fontSize: 13, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  const next = [
    "MD simulations (10 ns + MM-PBSA) on top CYP1B1 + SGLT2 hits to confirm pose stability and ΔG.",
    "PASS Online activity-class confirmation for the top-50 consensus compounds (third orthogonal channel).",
    "Wet-lab validation — biochemical assay on top 3 leads per target (CYP1B1 EROD, SGLT2 glucose uptake, MAO-B kynuramine).",
    "Submit to Frontiers in Pharmacology (Ethnopharmacology or In Silico Methods section).",
  ];
  let y2 = ry + 0.45;
  next.forEach((t, i) => {
    s.addShape(pres.shapes.OVAL, {
      x: rx, y: y2 + 0.06, w: 0.32, h: 0.32, fill: { color: C.amber }, line: { color: C.amber },
    });
    s.addText(String(i + 1), {
      x: rx, y: y2 + 0.06, w: 0.32, h: 0.32,
      fontFace: F.head, fontSize: 14, bold: true, color: C.navy,
      align: "center", valign: "middle", margin: 0,
    });
    s.addText(t, {
      x: rx + 0.42, y: y2, w: rw - 0.42, h: 1.0,
      fontFace: F.body, fontSize: 11, color: C.cream, margin: 0, valign: "top", lineSpacingMultiple: 1.25,
    });
    y2 += 1.05;
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: H - 0.35, w: W, h: 0.35, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText(`ANPDB exploratory novelty + multi-method target prediction      ·      ${TOTAL_PAGES} / ${TOTAL_PAGES}`, {
    x: 0.5, y: H - 0.35, w: W - 1, h: 0.35,
    fontFace: F.body, fontSize: 10, bold: true, color: C.navy, charSpacing: 6, valign: "middle", margin: 0,
  });
}

pres.writeFile({ fileName: OUT }).then(() => {
  console.log("wrote", OUT);
}).catch((e) => {
  console.error(e);
  process.exit(1);
});
