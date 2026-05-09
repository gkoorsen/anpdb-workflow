// Build the ANPDB novelty + T2D target-prediction summary deck.
const pptxgen = require("pptxgenjs");
const path = require("path");

const FIG = path.resolve(__dirname, "../output/figures");
const OUT = path.resolve(__dirname, "../output/anpdb_t2dm_summary.pptx");

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
};

const F = { head: "Cambria", body: "Calibri" };

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.3 x 7.5
pres.author = "ANPDB Novelty + T2D Target Prediction";
pres.title  = "ANPDB Novelty + T2D Target Prediction";
const W = 13.3, H = 7.5;

function addFooter(slide, pageLabel) {
  slide.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: H - 0.32, w: W, h: 0.32, fill: { color: C.navy }, line: { color: C.navy },
  });
  slide.addText("ANPDB novelty + T2D target prediction", {
    x: 0.5, y: H - 0.32, w: 6, h: 0.32, fontFace: F.body, fontSize: 9,
    color: C.cream, valign: "middle", margin: 0,
  });
  slide.addText(pageLabel, {
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
    x: 0.75, y: 0.74, w: 12, h: 0.55, fontFace: F.head, fontSize: 28, bold: true,
    color: C.ink, margin: 0,
  });
}

// ---------- Slide 1: Title ----------
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: H - 0.6, w: W, h: 0.6, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("CHEMICAL NOVELTY  /  TARGET PREDICTION  /  TYPE 2 DIABETES", {
    x: 1, y: 1.3, w: 11, h: 0.4, fontFace: F.body, fontSize: 12, bold: true,
    color: C.amber, charSpacing: 8,
  });
  s.addText("Novel compounds in the African Natural Products Database\nand their predicted Type 2 diabetes targets", {
    x: 1, y: 1.85, w: 11.3, h: 2.4, fontFace: F.head, fontSize: 36, bold: true,
    color: C.white, lineSpacingMultiple: 1.05,
  });
  s.addText([
    { text: "Pipeline summary    ·    ", options: { color: C.cream } },
    { text: "ANPDB v.2026  ·  COCONUT 05-2026  ·  ChEMBL 36  ·  PIDGINv4", options: { color: "B0C4DE" } },
  ], {
    x: 1, y: 4.6, w: 11, h: 0.5, fontFace: F.body, fontSize: 14,
  });
  s.addShape(pres.shapes.LINE, {
    x: 1, y: 5.25, w: 4, h: 0, line: { color: C.amber, width: 1.5 },
  });
  s.addText("Generated " + new Date().toISOString().slice(0,10), {
    x: 1, y: 5.35, w: 11, h: 0.4, fontFace: F.body, fontSize: 11,
    color: "9AAFCC", italic: true,
  });
}

// ---------- Slide 2: The question ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "THE QUESTION", "Which ANPDB compounds are structurally novel?");

  // Left column — operational definition
  s.addText([
    { text: "Operational definition", options: { bold: true, color: C.ink, fontSize: 16, breakLine: true } },
    { text: "An ANPDB compound is treated as ", options: { color: C.ink } },
    { text: "novel ", options: { bold: true, color: C.amber } },
    { text: "if its 14-character InChIKey skeleton block does not appear in either reference natural-product or bioactivity database. A second pass uses Morgan-ECFP4 Tanimoto similarity to flag close analogues that the exact match misses.",
      options: { color: C.ink } },
  ], { x: 0.75, y: 1.7, w: 6.0, h: 3.5, fontFace: F.body, fontSize: 14, lineSpacingMultiple: 1.25, valign: "top" });

  // Right column — three contrasting framings card
  const cardX = 7.2, cardW = 5.4;
  const items = [
    { tag: "STRUCTURAL", text: "InChIKey skeleton match — fast, deterministic, ignores stereochemistry." },
    { tag: "FUZZY",      text: "Morgan-ECFP4 Tanimoto (radius 2, 2048 bits, threshold 0.85) — catches close analogues." },
    { tag: "STANDARDISED", text: "Salt-stripped + neutralised parent structures (planned cleanup step)." },
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
      fontFace: F.body, fontSize: 13, color: C.ink, margin: 0, valign: "top",
    });
    y += 1.55;
  }
  addFooter(s, "2");
}

// ---------- Slide 3: Datasets ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "INPUTS", "Three databases anchor the analysis");

  const col = [
    { name: "ANPDB",  size: "11,448",  unit: "natural products",
      sub:  "Northern + Eastern + Southern Africa, curated by U Freiburg (NAR 2026).",
      role: "QUERY SET" },
    { name: "COCONUT", size: "479,721",  unit: "InChIKey skeletons",
      sub:  "Open natural-product collection (May 2026 release, 738k structures, ChEMBL-curation pipeline).",
      role: "NP REFERENCE" },
    { name: "ChEMBL 36", size: "2,715,471", unit: "InChIKey skeletons",
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
      fontFace: F.head, fontSize: 38, bold: true, color: i===0 ? C.amber : C.blue, margin: 0,
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
  addFooter(s, "3");
}

// ---------- Slide 4: Pipeline ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "PIPELINE", "Three sequential passes, increasingly stringent");

  const steps = [
    { n: "1", title: "Exact match",
      detail: "InChIKey skeleton block of every ANPDB compound is looked up in COCONUT and ChEMBL.",
      out:    "1,370 unmatched (12.0%)" },
    { n: "2", title: "Tanimoto fuzzy",
      detail: "Morgan-ECFP4 (radius 2, 2048 bits) Tanimoto vs every COCONUT compound. Threshold 0.85.",
      out:    "1,014 truly novel (8.9%)" },
    { n: "3", title: "Target prediction",
      detail: "Pre-filtered ANPDB (Lipinski + PAINS + ADMET) scored against PIDGINv4, two AD cuts (60th and 90th percentile), with disease annotations merged in.",
      out:    "20 T2D-related targets" },
  ];
  const sx = 0.75, sy = 1.85, sw = 11.8, sh = 1.55;
  steps.forEach((st, i) => {
    const y = sy + i * (sh + 0.18);
    s.addShape(pres.shapes.RECTANGLE, {
      x: sx, y, w: sw, h: sh,
      fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.OVAL, {
      x: sx + 0.3, y: y + 0.32, w: 0.9, h: 0.9, fill: { color: C.navy }, line: { color: C.navy },
    });
    s.addText(st.n, {
      x: sx + 0.3, y: y + 0.32, w: 0.9, h: 0.9,
      fontFace: F.head, fontSize: 28, bold: true, color: C.amber,
      align: "center", valign: "middle", margin: 0,
    });
    s.addText(st.title, {
      x: sx + 1.4, y: y + 0.18, w: 5.0, h: 0.5,
      fontFace: F.head, fontSize: 20, bold: true, color: C.ink, margin: 0,
    });
    s.addText(st.detail, {
      x: sx + 1.4, y: y + 0.7, w: 6.5, h: 0.85,
      fontFace: F.body, fontSize: 12, color: C.muted, margin: 0, valign: "top",
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: sx + sw - 3.4, y: y + 0.3, w: 3.0, h: sh - 0.6, fill: { color: C.cream }, line: { color: C.divider, width: 0.5 },
    });
    s.addText("OUTPUT", {
      x: sx + sw - 3.4, y: y + 0.34, w: 3.0, h: 0.3,
      fontFace: F.body, fontSize: 9, bold: true, color: C.muted, charSpacing: 6, align: "center", margin: 0,
    });
    s.addText(st.out, {
      x: sx + sw - 3.4, y: y + 0.65, w: 3.0, h: 0.7,
      fontFace: F.head, fontSize: 17, bold: true, color: C.amber, align: "center", margin: 0,
    });
  });
  addFooter(s, "4");
}

// ---------- Slide 5: Pass 1 results ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "PASS 1 — EXACT INCHIKEY MATCH", "12.0 % of ANPDB has no skeleton match anywhere");

  // Big number callout
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 1.85, w: 5.5, h: 4.5, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 1.85, w: 5.5, h: 0.18, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("UNMATCHED COMPOUNDS", {
    x: 0.95, y: 2.15, w: 5.1, h: 0.4,
    fontFace: F.body, fontSize: 11, bold: true, color: C.amber, charSpacing: 8, margin: 0,
  });
  s.addText("1,370", {
    x: 0.95, y: 2.6, w: 5.1, h: 1.8,
    fontFace: F.head, fontSize: 96, bold: true, color: C.white, margin: 0,
  });
  s.addText("of 11,448 ANPDB compounds", {
    x: 0.95, y: 4.35, w: 5.1, h: 0.4,
    fontFace: F.body, fontSize: 14, color: "B0C4DE", margin: 0,
  });
  s.addText("12.0 %", {
    x: 0.95, y: 4.85, w: 5.1, h: 0.7,
    fontFace: F.head, fontSize: 40, bold: true, color: C.amber, margin: 0,
  });
  s.addText("Have no exact InChIKey-skeleton match in COCONUT or ChEMBL.", {
    x: 0.95, y: 5.6, w: 5.1, h: 0.6,
    fontFace: F.body, fontSize: 12, color: "CADCFC", italic: true, margin: 0, valign: "top",
  });

  // Right-side table
  const tx = 6.7, ty = 1.85;
  const rows = [
    [{ text: "Set", options: { bold: true, color: C.white, fill: { color: C.blue } } },
     { text: "Count", options: { bold: true, color: C.white, fill: { color: C.blue }, align: "right" } },
     { text: "%", options: { bold: true, color: C.white, fill: { color: C.blue }, align: "right" } }],
    [ "ANPDB (total)",        { text: "11,448", options: { align: "right" } }, { text: "100.0", options: { align: "right" } }],
    [ "  in ChEMBL",          { text: "4,808",  options: { align: "right" } }, { text: "42.0",  options: { align: "right" } }],
    [ "  in COCONUT",         { text: "10,020", options: { align: "right" } }, { text: "87.5",  options: { align: "right" } }],
    [ "  in either",          { text: "10,078", options: { align: "right" } }, { text: "88.0",  options: { align: "right" } }],
    [ "  in both",            { text: "4,750",  options: { align: "right" } }, { text: "41.5",  options: { align: "right" } }],
    [ { text: "  novel (in neither)", options: { bold: true, color: C.amber } },
      { text: "1,370",  options: { align: "right", bold: true, color: C.amber } },
      { text: "12.0",   options: { align: "right", bold: true, color: C.amber } }],
  ];
  s.addTable(rows, {
    x: tx, y: ty, w: 5.9, colW: [3.0, 1.6, 1.3],
    rowH: 0.42,
    fontFace: F.body, fontSize: 13, color: C.ink,
    border: { type: "solid", pt: 0.5, color: C.divider },
    fill: { color: C.white },
  });
  s.addText("Notes", {
    x: tx, y: 5.0, w: 5.9, h: 0.35,
    fontFace: F.body, fontSize: 11, bold: true, color: C.muted, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "•  Match key: 14-char InChIKey skeleton (stereochemistry-agnostic).", options: { breakLine: true } },
    { text: "•  COCONUT skeletons: 479,721 distinct.", options: { breakLine: true } },
    { text: "•  ChEMBL 36 skeletons: 2,715,471 distinct.", options: {} },
  ], {
    x: tx, y: 5.4, w: 5.9, h: 1.5,
    fontFace: F.body, fontSize: 12, color: C.ink, margin: 0, lineSpacingMultiple: 1.3,
  });
  addFooter(s, "5");
}

// ---------- Slide 6: Pass 1 caveats ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "PASS 1 — INTERPRETATION", "Three reasons the 1,370 number is conservative");

  const items = [
    { title: "Salt forms",
      body: "ANPDB stores some entries as sodium / chloride / sulphate adducts; COCONUT stores the parent. Their InChIKey blocks differ.",
      ex: "Mol_00032 vanilloside — has counter-ion '.C…'; parent likely already in COCONUT." },
    { title: "Synthetic derivatives",
      body: "Mosher-ester (α-methoxy-α-trifluoromethyl-phenylacetate) derivatives appear because they were used to assign stereochemistry, not because they are natural.",
      ex: "Mol_02941, Mol_02943, Mol_03208 — Mosher-ester derivatives of khayanone." },
    { title: "Tautomers / charge states",
      body: "Different tautomers of the same compound have distinct InChIKey blocks even when the underlying skeleton is the same.",
      ex: "Sulphates and free acids of the same flavonoid often appear separately." },
  ];
  let y = 1.85;
  items.forEach((it, i) => {
    s.addShape(pres.shapes.RECTANGLE, {
      x: 0.75, y, w: 11.8, h: 1.55,
      fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
      shadow: { type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.06 },
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: 0.75, y, w: 0.12, h: 1.55, fill: { color: C.amber }, line: { color: C.amber },
    });
    s.addText(it.title, {
      x: 1.0, y: y + 0.15, w: 4.5, h: 0.45,
      fontFace: F.head, fontSize: 18, bold: true, color: C.ink, margin: 0,
    });
    s.addText(it.body, {
      x: 1.0, y: y + 0.6, w: 7.0, h: 0.95,
      fontFace: F.body, fontSize: 12, color: C.ink, margin: 0, valign: "top",
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: 8.3, y: y + 0.18, w: 4.0, h: 1.2, fill: { color: C.cream }, line: { color: C.divider, width: 0.5 },
    });
    s.addText("EXAMPLE", {
      x: 8.45, y: y + 0.22, w: 3.7, h: 0.25,
      fontFace: F.body, fontSize: 9, bold: true, color: C.muted, charSpacing: 6, margin: 0,
    });
    s.addText(it.ex, {
      x: 8.45, y: y + 0.5, w: 3.75, h: 0.85,
      fontFace: F.body, fontSize: 11, color: C.ink, italic: true, margin: 0, valign: "top",
    });
    y += 1.7;
  });
  addFooter(s, "6");
}

// ---------- Slide 7: Pass 2 ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "PASS 2 — TANIMOTO FUZZY MATCH", "1,014 compounds remain novel under similarity-based search");

  s.addText([
    { text: "Method", options: { bold: true, color: C.ink, fontSize: 16, breakLine: true } },
    { text: "For each of the 1,370 Pass-1 residual compounds, computed Morgan circular fingerprint (radius 2, 2,048 bits) and compared against all 738,823 fingerprintable COCONUT structures using bulk Tanimoto.\n\n",
      options: { color: C.ink } },
    { text: "Threshold: ", options: { bold: true, color: C.ink } },
    { text: "Tanimoto ≥ 0.85 flags a close analogue. Below this, the compound is treated as ", options: { color: C.ink } },
    { text: "truly novel.", options: { bold: true, color: C.amber } },
  ], {
    x: 0.75, y: 1.85, w: 6.0, h: 4.5, fontFace: F.body, fontSize: 13, color: C.ink, lineSpacingMultiple: 1.3, valign: "top",
  });

  const rx = 7.0, ry = 1.85;
  const buckets = [
    ["<0.50",      91],
    ["0.50–0.70", 418],
    ["0.70–0.80", 339],
    ["0.80–0.85", 166],
    ["0.85–0.90", 147],
    ["0.90–0.95",  84],
    ["0.95–1.00", 125],
  ];
  s.addChart(pres.charts.BAR, [{
    name: "Compounds", labels: buckets.map(b => b[0]), values: buckets.map(b => b[1]),
  }], {
    x: rx, y: ry, w: 5.55, h: 4.4, barDir: "bar",
    chartColors: [C.blue],
    chartArea: { fill: { color: C.white }, roundedCorners: false },
    plotArea:  { fill: { color: C.white } },
    catAxisLabelColor: C.muted, valAxisLabelColor: C.muted,
    valGridLine: { color: C.divider, size: 0.5 },
    catGridLine: { style: "none" },
    showValue: true, dataLabelPosition: "outEnd", dataLabelColor: C.ink,
    showLegend: false, showTitle: true, title: "Tanimoto distribution vs COCONUT (n=1,370)", titleFontFace: F.head, titleFontSize: 12, titleColor: C.ink,
    catAxisLabelFontFace: F.body, valAxisLabelFontFace: F.body,
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 6.4, w: 11.8, h: 0.6, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText([
    { text: "356 ", options: { color: C.amber, fontFace: F.head, bold: true, fontSize: 18 } },
    { text: "have a close analogue (Tanimoto ≥ 0.85)    ·    ", options: { color: C.white, fontSize: 14 } },
    { text: "1,014 ", options: { color: C.amber, fontFace: F.head, bold: true, fontSize: 18 } },
    { text: "remain truly novel under fuzzy matching", options: { color: C.white, fontSize: 14 } },
  ], { x: 0.95, y: 6.42, w: 11.5, h: 0.55, valign: "middle", fontFace: F.body, margin: 0 });
  addFooter(s, "7");
}

// ---------- Slide 8: PIDGIN context ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "PIVOT — TARGET PREDICTION", "Drug-likeness pre-filter narrows the input pool");

  s.addText([
    { text: "Workflow inputs to PIDGINv4", options: { bold: true, color: C.ink, fontSize: 16, breakLine: true } },
    { text: "ANPDB compounds were filtered for ", options: { color: C.ink } },
    { text: "Lipinski's Rule of Five, PAINS reactivity, and ADMET", options: { bold: true, color: C.amber } },
    { text: " before being scored against PIDGINv4. Two applicability-domain (AD) cutoffs are reported throughout: ", options: { color: C.ink } },
    { text: "ad60", options: { bold: true } }, { text: " (≥ 60th percentile, looser) and ", options: {} },
    { text: "ad90", options: { bold: true } }, { text: " (≥ 90th percentile, strict).", options: {} },
  ], { x: 0.75, y: 1.85, w: 12, h: 1.5, fontFace: F.body, fontSize: 14, color: C.ink, lineSpacingMultiple: 1.3, valign: "top" });

  // Three stat tiles
  const tiles = [
    { label: "PIDGINv4 prediction rows", val: "1,546,201", note: "all (compound × target × threshold) tuples" },
    { label: "Hits at ad60 (looser)",    val: "1,779",     note: "391 compounds hit 226 targets" },
    { label: "Hits at ad90 (strict)",    val: "247",       note: "140 compounds hit 37 targets" },
  ];
  const tx0 = 0.75, ty0 = 3.4, tw = 3.95, th = 2.2, tg = 0.18;
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
      x: x + 0.3, y: ty0 + 0.25, w: tw - 0.6, h: 0.4,
      fontFace: F.body, fontSize: 11, bold: true, color: C.muted, charSpacing: 6, margin: 0,
    });
    s.addText(t.val, {
      x: x + 0.3, y: ty0 + 0.7, w: tw - 0.6, h: 1.0,
      fontFace: F.head, fontSize: 38, bold: true, color: i===2 ? C.amber : C.blue, margin: 0,
    });
    s.addText(t.note, {
      x: x + 0.3, y: ty0 + 1.7, w: tw - 0.6, h: 0.45,
      fontFace: F.body, fontSize: 11, color: C.muted, margin: 0, italic: true, valign: "top",
    });
  });

  s.addText("Of those hits, 148 (ad60) and 21 (ad90) reach a target whose disease annotations include T2D-related terms.", {
    x: 0.75, y: 5.85, w: 11.8, h: 0.6,
    fontFace: F.body, fontSize: 13, color: C.ink, italic: true, margin: 0, valign: "top",
  });
  addFooter(s, "8");
}

// ---------- Slide 9: Broad T2D landscape (fig1) ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "T2D-RELATED TARGETS — BROAD VIEW", "20 targets at ad60 mention T2D anywhere in their diseases");
  s.addImage({
    path: path.join(FIG, "fig1_t2d_targets.png"),
    x: 0.75, y: 1.6, w: 8.5, h: 5.4, sizing: { type: "contain", w: 8.5, h: 5.4 },
  });
  // Right-side commentary card
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.5, y: 1.6, w: 3.1, h: 5.4, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.5, y: 1.6, w: 3.1, h: 0.45, fill: { color: C.teal }, line: { color: C.teal },
  });
  s.addText("OBSERVATION", {
    x: 9.65, y: 1.66, w: 2.85, h: 0.32,
    fontFace: F.body, fontSize: 10, bold: true, color: C.white, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "P05177 (CYP1A2) ", options: { bold: true, breakLine: true } },
    { text: "dominates the broad list with 42 compounds — but it is a promiscuous metabolizer associated with many diseases beyond T2D.\n\n", options: {} },
    { text: "Eight targets survive the ad90 cutoff. Most novel-set compounds hit P05177; only a handful reach more T2D-specific targets.", options: {} },
  ], {
    x: 9.65, y: 2.15, w: 2.8, h: 4.7,
    fontFace: F.body, fontSize: 12, color: C.ink, lineSpacingMultiple: 1.3, margin: 0, valign: "top",
  });
  addFooter(s, "9");
}

// ---------- Slide 10: T2D-association strength (fig0) ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "T2D-ASSOCIATION STRENGTH", "A simple, reproducible per-target score");

  s.addText([
    { text: "Score = (count of T2D-related diseases on the target) + 2 if 'Diabetes Mellitus, Non-Insulin-Dependent' is named explicitly.", options: { color: C.ink, breakLine: true } },
    { text: "Filter applied throughout the rest of the deck: ", options: { color: C.ink } },
    { text: "strength ≥ 4", options: { bold: true, color: C.amber } },
    { text: " — keeps targets where the strict NIDDM term coexists with at least one comorbid condition (or four comorbid signals without the strict term).", options: {} },
  ], {
    x: 0.75, y: 1.55, w: 12, h: 1.0, fontFace: F.body, fontSize: 13, color: C.ink, margin: 0, lineSpacingMultiple: 1.3, valign: "top",
  });

  s.addImage({
    path: path.join(FIG, "fig0_t2d_target_strength.png"),
    x: 0.75, y: 2.65, w: 8.5, h: 4.4, sizing: { type: "contain", w: 8.5, h: 4.4 },
  });

  // Right summary
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 2.65, w: 3.05, h: 4.4, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 2.65, w: 3.05, h: 0.12, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("FILTER OUTCOME", {
    x: 9.7, y: 2.85, w: 2.85, h: 0.35,
    fontFace: F.body, fontSize: 10, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "6 / 20", options: { bold: true, color: C.white, fontSize: 36, fontFace: F.head, breakLine: true } },
    { text: "T2D-relevant targets retained\n", options: { color: "B0C4DE", fontSize: 12, breakLine: true } },
    { text: "\n26 ", options: { color: C.amber, bold: true, fontSize: 22, fontFace: F.head } },
    { text: "(compound, target) hits at ad60\n", options: { color: C.white, fontSize: 12, breakLine: true } },
    { text: "5 ", options: { color: C.amber, bold: true, fontSize: 22, fontFace: F.head } },
    { text: "(compound, target) hits at ad90", options: { color: C.white, fontSize: 12 } },
  ], { x: 9.7, y: 3.25, w: 2.85, h: 3.7, fontFace: F.body, valign: "top", margin: 0 });

  addFooter(s, "10");
}

// ---------- Slide 11: Strength-filtered targets (fig1s) ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "STRENGTH-FILTERED TARGETS", "After strength ≥ 4: six targets with strong T2D evidence");
  s.addImage({
    path: path.join(FIG, "fig1s_t2d_targets_strong.png"),
    x: 0.75, y: 1.6, w: 8.5, h: 4.0, sizing: { type: "contain", w: 8.5, h: 4.0 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 4.0, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 0.45, fill: { color: C.teal }, line: { color: C.teal },
  });
  s.addText("READING", {
    x: 9.7, y: 1.66, w: 2.85, h: 0.32,
    fontFace: F.body, fontSize: 10, bold: true, color: C.white, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "CYP19A1 / aromatase (P11511) ", options: { bold: true, breakLine: true } },
    { text: "leads with 12 compounds.\n", options: { breakLine: true } },
    { text: "ACE (P12821), MMP9 (P14780), VEGFA (P15692), PTP1B (P18031) and NF-κB p65 / RELA (Q04206) round out the six. PTP1B in particular is a flagship validated T2D drug target.", options: {} },
  ], {
    x: 9.7, y: 2.15, w: 2.85, h: 4.3,
    fontFace: F.body, fontSize: 12, color: C.ink, lineSpacingMultiple: 1.25, margin: 0, valign: "top",
  });

  // Bottom strip
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.75, y: 5.85, w: 11.85, h: 1.2, fill: { color: C.navy }, line: { color: C.navy },
  });
  s.addText([
    { text: "Key finding   ", options: { bold: true, color: C.amber, fontSize: 13, charSpacing: 5 } },
    { text: "None of the six strength-filtered T2D targets is hit by a compound from the 1,370 Pass-1 novel set.\n", options: { color: C.white, fontSize: 13 } },
    { text: "Implication   ", options: { bold: true, color: C.amber, fontSize: 13, charSpacing: 5 } },
    { text: "The strongest T2D leads from this pool come from already-known scaffolds; the novelty signal sits elsewhere.", options: { color: "CADCFC", fontSize: 13 } },
  ], { x: 0.95, y: 5.95, w: 11.5, h: 1.05, fontFace: F.body, valign: "top", margin: 0, lineSpacingMultiple: 1.3 });

  addFooter(s, "11");
}

// ---------- Slide 12: Probability vs AD ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "CONFIDENCE OF T2D PREDICTIONS", "ad90 promotions cluster in the high-confidence corner");
  s.addImage({
    path: path.join(FIG, "fig3_proba_vs_ad.png"),
    x: 0.75, y: 1.6, w: 8.5, h: 5.4, sizing: { type: "contain", w: 8.5, h: 5.4 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 5.4, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 0.45, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("READING", {
    x: 9.7, y: 1.66, w: 2.85, h: 0.32,
    fontFace: F.body, fontSize: 10, bold: true, color: C.white, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "Each dot is one (compound, T2D-related target) prediction at ad60. Red dots also pass ad90.\n\n", options: { breakLine: false } },
    { text: "Defensible hits", options: { bold: true, breakLine: true } },
    { text: "concentrate above probability 0.7 AND applicability-domain percentile 90 — exactly where you would want a clinical-grade prediction to live.", options: {} },
  ], {
    x: 9.7, y: 2.15, w: 2.85, h: 4.7,
    fontFace: F.body, fontSize: 12, color: C.ink, lineSpacingMultiple: 1.3, margin: 0, valign: "top",
  });
  addFooter(s, "12");
}

// ---------- Slide 13: Top compounds (filtered) ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "TOP T2D-ACTIVE COMPOUNDS", "Strength-filtered compound leaderboard");
  s.addImage({
    path: path.join(FIG, "fig4s_top_compounds_strong.png"),
    x: 0.75, y: 1.6, w: 8.5, h: 5.4, sizing: { type: "contain", w: 8.5, h: 5.4 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 5.4, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 0.45, fill: { color: C.teal }, line: { color: C.teal },
  });
  s.addText("LEADERBOARD", {
    x: 9.7, y: 1.66, w: 2.85, h: 0.32,
    fontFace: F.body, fontSize: 10, bold: true, color: C.white, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "Mol_12233 ", options: { bold: true, color: C.amber } },
    { text: "leads — hits 3 of the 6 strength-filtered T2D targets at ad60 with max p≈0.64.\n\n", options: {} },
    { text: "16 other compounds hit a single high-strength target each, with prediction probabilities ranging 0.54–0.91. None are in the 1,370 novel set.", options: {} },
  ], {
    x: 9.7, y: 2.15, w: 2.85, h: 4.7,
    fontFace: F.body, fontSize: 12, color: C.ink, lineSpacingMultiple: 1.3, margin: 0, valign: "top",
  });
  addFooter(s, "13");
}

// ---------- Slide 14: Compound × target heatmap (broad view) ----------
{
  const s = pres.addSlide();
  s.background = { color: C.cream };
  addContentTitle(s, "COMPOUND × TARGET DETAIL", "Probability heatmap, broad T2D-related set (ad60)");
  s.addImage({
    path: path.join(FIG, "fig2_compound_target_heatmap.png"),
    x: 1.2, y: 1.55, w: 8.0, h: 5.5, sizing: { type: "contain", w: 8.0, h: 5.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 5.5, fill: { color: C.white }, line: { color: C.divider, width: 0.5 },
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: 9.55, y: 1.6, w: 3.05, h: 0.45, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("READING", {
    x: 9.7, y: 1.66, w: 2.85, h: 0.32,
    fontFace: F.body, fontSize: 10, bold: true, color: C.white, charSpacing: 6, margin: 0,
  });
  s.addText([
    { text: "Cells coloured by max prediction probability (0.5 → 1.0).\n\n", options: {} },
    { text: "Compound rows ending '●'", options: { bold: true } },
    { text: " are members of the 1,370 Pass-1 novel set.\n\n", options: { breakLine: false } },
    { text: "Columns are the 12 most-hit T2D-related targets in the broad view.", options: {} },
  ], {
    x: 9.7, y: 2.15, w: 2.85, h: 4.8,
    fontFace: F.body, fontSize: 12, color: C.ink, lineSpacingMultiple: 1.3, margin: 0, valign: "top",
  });
  addFooter(s, "14");
}

// ---------- Slide 15: Summary & next ----------
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

  // Two columns
  const lx = 0.75, ly = 1.85, lw = 5.9;
  const rx = 7.0, ry = 1.85, rw = 5.6;

  s.addText("Headline numbers", {
    x: lx, y: ly, w: lw, h: 0.4, fontFace: F.body, fontSize: 13, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  const lines = [
    ["1,370",  "ANPDB compounds with no exact InChIKey match in COCONUT or ChEMBL"],
    ["1,014",  "remain truly novel under Tanimoto ≥ 0.85 fuzzy matching"],
    ["6",      "T2D-relevant targets pass the strength ≥ 4 filter"],
    ["26 / 5", "(compound, target) hits at ad60 / ad90 against those 6 targets"],
  ];
  let y = ly + 0.45;
  lines.forEach(([num, label]) => {
    s.addText(num, {
      x: lx, y, w: 2.0, h: 0.55, fontFace: F.head, fontSize: 26, bold: true, color: C.amber, margin: 0, align: "left",
    });
    s.addText(label, {
      x: lx + 2.05, y: y + 0.08, w: lw - 2.05, h: 0.55,
      fontFace: F.body, fontSize: 12, color: C.cream, margin: 0, valign: "top",
    });
    y += 0.65;
  });

  s.addText("Recommended next steps", {
    x: rx, y: ry, w: rw, h: 0.4, fontFace: F.body, fontSize: 13, bold: true, color: C.amber, charSpacing: 6, margin: 0,
  });
  const next = [
    "Run salt-strip + tautomer canonicalisation on ANPDB before keying — should cut the 1,370 number further and remove Mosher-ester artefacts.",
    "Cross-validate the 6 strength-filtered targets against an external resource (OpenTargets, DisGeNET) for genetic-evidence weighting.",
    "Re-run PIDGIN on the 1,014 truly-novel set explicitly — the current pre-filtered pool may have excluded some.",
    "For Mol_12233 specifically: literature search and / or experimental T2D assay against P11511, P14780, P12821.",
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
      fontFace: F.body, fontSize: 12, color: C.cream, margin: 0, valign: "top", lineSpacingMultiple: 1.25,
    });
    y2 += 1.0;
  });

  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: H - 0.35, w: W, h: 0.35, fill: { color: C.amber }, line: { color: C.amber },
  });
  s.addText("ANPDB novelty + T2D target prediction      ·      Pipeline summary      ·      15 / 15", {
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
