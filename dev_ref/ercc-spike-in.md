# ERCC spike-ins: how they work, what to expect, and where they bite

A developer reference for the `rnaseq_dashboard` project. Captures the domain
knowledge behind the **Spike-in (ERCC) titration QC** (P3d) and the **spike-in
filtering criteria** (P3e) so it doesn't get re-derived. Companion to
[normalization-scran-vs-deseq.md](normalization-scran-vs-deseq.md) and the
`rnaseq-bioc` skill.

## TL;DR

- ERCC spike-ins are **92 synthetic transcripts added at *known* concentrations**.
  Because the input is known, observed-vs-known is a **titration check**: a healthy
  library tracks the known concentration linearly (log-log **slope ≈ 1**, **high
  R²**) across a wide dynamic range.
- In this app spike-ins are **QC-only**. Size factors are estimated on endogenous
  features (`controlGenes = which(feature_class == "endogenous")`); spikes are
  **flagged, never dropped**, and **never** drive normalization (see "Why QC-only").
- **Bulk:** generally informative, *if* the spike amount is proportional to input
  RNA / cell number. **Single-cell:** unreliable in droplet platforms, and spike-in
  *normalization* is discouraged — prefer deconvolution (`scran`).

## What ERCC spike-ins are

The **External RNA Controls Consortium** (ERCC, hosted at NIST) defined a set of
**92 polyadenylated synthetic RNA transcripts**, ~250–2000 nt, with sequences
designed to be unlike any common reference genome so reads map unambiguously.
They are sold as two formulations:

- **Mix 1** and **Mix 2** contain the *same* 92 transcripts but at *different*
  concentrations, arranged into **4 subgroups** with fixed Mix1:Mix2 ratios
  (**4:1, 1:1, 0.67:1, 0.5:1**). Spiking Mix 1 into one condition and Mix 2 into
  another lets you also check **fold-change accuracy**, not just abundance.
- Within a mix, the nominal concentrations span roughly **2²⁰ (~10⁶, six orders of
  magnitude)** — from a few attomoles/µL down to tiny amounts — deliberately
  covering the range of real transcript abundances.
- Spikes are typically added at **~1–2% of total RNA mass** (protocol-dependent),
  *before* library prep, so they travel through reverse transcription, amplification
  and sequencing alongside the endogenous RNA.

The bundled reference (`inst/extdata/ercc_concentrations.csv`, provenance in
`ERCC_SOURCE.md`) holds the 92 ids with their Mix 1 / Mix 2 concentrations; the app
joins observed expression to these by id, or to a user-designated
`rowData$spike_concentration` column.

## The dose-response (what P3d/P3e compute)

Per sample, fit `lm(log10(observed) ~ log10(known concentration))` over the spike
features, where **observed** is a **linear, depth-normalized assay** (TPM/FPKM when
available, else CPM — *never* raw counts or logcounts). Zeros are dropped before the
log (no pseudo-count); the fit needs **≥ 3 detected points** or slope/R² are `NA`.
Read it as:

- **slope ≈ 1** — observed abundance scales proportionally with input. Slope well
  below 1 = compressed dynamic range (saturation / loss of low-abundance spikes);
  above 1 = over-dispersed.
- **high R² (≳ 0.9)** — tight linear relationship; a healthy titration.
- **lowest detected concentration (LOD)** — the smallest spiked concentration still
  observed; a **sensitivity** read. Its units depend on the source (attomoles/µL for
  the bundled Mix, arbitrary for a custom column), so it has **no portable fixed
  threshold** — judge it relative to the cohort.
- **% spike-in** and **# detected spikes** — content reads. Very **high** % spike-in
  → low endogenous complexity (degraded / low-input RNA, or over-spiking); very
  **low** % spike-in → under-spiking / pipetting loss, which also makes the fit
  unreliable. Hence P3e flags % spike-in with a **two-sided** fence.

## Target detection range to aim for

- A good bulk library detects spikes across **most of the ~6-order dynamic range**;
  the low-concentration spikes are the first to fall below detection, so a rising
  LOD across samples flags declining sensitivity.
- Linear region with **slope ≈ 1** and **R² ≳ 0.9**; the bottom of the curve flattens
  as counts hit shot-noise / zero — that floor is expected, not a defect.
- For **fold-change** checks (Mix1 vs Mix2 subgroups), observed log-ratios should
  track the nominal 4:1 / 1:1 / 0.67:1 / 0.5:1 ratios; compression toward 1:1 at low
  abundance is normal.

## Why QC-only here (not normalization)

Spike-in normalization assumes the spiked amount reflects a **constant per-cell /
per-sample reference**. That holds only when spikes are added proportionally to cell
number (or when a genuine **global** shift in total RNA content is the signal you
want to preserve — e.g. Lovén et al. 2012, c-Myc amplification). In the general case,
spike fraction varies with input mass, RNA quality and pipetting, so using spikes as
the normalization basis injects technical noise. This app therefore keeps size
factors on endogenous genes (DESeq2 median-of-ratios, `controlGenes`) and uses spikes
purely as an **independent QC witness**.

## Concerns when using ERCC spike-ins

### Bulk RNA-seq
- **Proportionality.** Spiking by RNA *mass* assumes comparable total RNA per sample;
  if biological total-RNA content differs, the spike fraction shifts for reasons
  unrelated to quality. Spiking by **cell number** is more defensible.
- **Composition bias.** ERCCs differ from endogenous mRNA in **GC content, length,
  and the absence of a 5′ cap / introns**, so their capture and amplification
  efficiency isn't identical — expect a systematic offset (intercept), which is why
  the **slope**, not the absolute level, is the robust quality metric.
- **Pipetting / dilution variability** dominates at the low-concentration end.
- **Small-n.** With 3–12 bulk samples, treat spike-derived flags as **advisory**, not
  auto-drop (consistent with the app's MAD-based sample flagging philosophy).

### Single-cell RNA-seq
- **Low, variable capture.** Per-cell input is tiny; ERCC capture efficiency is low
  and varies cell-to-cell, so spike counts are noisy.
- **Platform.** Spike-ins are workable in **plate-based** protocols (e.g. Smart-seq2)
  but **unreliable / rarely used in droplet** platforms (10x), where they add cost and
  behave inconsistently; many droplet workflows omit them entirely.
- **Normalization is discouraged.** Spike-in–based size factors for scRNA-seq were
  shown to be unreliable (Lun et al. 2017, *"Assessing the reliability of spike-in
  normalization..."*); the app's single-cell path (future P7) prefers **scran**
  pooling/deconvolution. Spikes remain useful for **technical-noise modelling**
  (Brennecke et al. 2013) and **cell-size / total-RNA** normalization *with caveats*.
- Ties into the bulk-first decision: when an `sce` is aggregated to **pseudobulk**,
  spike-in QC behaves like the bulk case again.

## Key references

- External RNA Controls Consortium / NIST — ERCC control design and concentrations
  (SRM 2374; Thermo ERCC Spike-In Mixes).
- Jiang L. et al. (2011) *Synthetic spike-in standards for RNA-seq experiments.*
  Genome Research.
- Lovén J. et al. (2012) *Revisiting global gene expression analysis.* Cell —
  spike-in normalization for global transcriptional shifts.
- Risso D. et al. (2014) *Normalization of RNA-seq data using factor analysis of
  control genes or samples (RUVg).* Nat. Biotechnol.
- Brennecke P. et al. (2013) *Accounting for technical noise in single-cell RNA-seq.*
  Nat. Methods.
- Lun A.T.L. et al. (2017) *Assessing the reliability of spike-in normalization for
  analyses of single-cell RNA sequencing data.* Genome Research.
