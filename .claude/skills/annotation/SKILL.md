---
name: annotation
description: Row-annotation (rowData) conventions for this RNA-seq dashboard — mapping feature ids to names/metadata from a Bioconductor OrgDb and from a user GTF. Covers id-type detection + Ensembl version-stripping, the OrgDb column→rowData mapping, GTF match-column resolution and attribute import (GTF authoritative over OrgDb), fill-where-matched-never-wipe overlay semantics, in_orgdb/in_gtf match flags + coverage banners, the feature-unit (meta$feature_type) that names the <unit>_name column, and the draft-composition wiring. Use before writing or reviewing any annotation/mapping logic on the Feature info page or in annotation_helpers.R / gtf_helpers.R.
metadata:
  type: project
  version: "1.0"
---

# Row annotation (OrgDb + GTF)

How the app turns opaque feature ids (`ENSMUSG…`, Entrez, symbols) into readable
`rowData` columns. The two sources are a Bioconductor **OrgDb** and a user-supplied
**GTF**. This skill is the *write/mapping* contract; `rnaseq-bioc` owns the adjacent
domain facts (`feature_class`, `feature_length` + the GTF union-exon length route, ERCC),
and `shiny-module` owns the draft/state mechanics. Pure logic lives in
`R/annotation_helpers.R` (OrgDb) and `R/gtf_helpers.R` (GTF); the Feature-info module
(`R/mod_feature.R`) wires them onto the draft.

## Feature ids and the feature unit

- **`detect_id_type(ids)`** majority-votes `"ensembl"` (`^ENS`), `"entrez"` (all digits),
  else `"symbol"`. The UI offers Auto / Ensembl / Entrez / Symbol.
- **Ensembl keys are version-stripped** (`sub("\\..*$","")` / `.strip_version()`) on *both*
  sides before matching — `ENSMUSG00000000001.5` matches `ENSMUSG00000000001`.
- **`meta$feature_type`** (gene/transcript/exon/feature; see `shiny-module` state model)
  decides the name column: symbols land in **`<feature_type>_name`** (e.g. `gene_name`).
  `lookup_feature()` and the adaptive labels read the same column.

## OrgDb mapping (`annotate_with_orgdb`)

Organism → package: `mouse`=`org.Mm.eg.db`, `human`=`org.Hs.eg.db` (Suggests; loaded via
`requireNamespace`). Maps the chosen `columns` (OrgDb keytypes) through
`AnnotationDbi::mapIds(multiVals = "first")`; columns absent from the OrgDb are skipped.
Each keytype has a **fixed** `rowData` target (`.orgdb_target()`):

| OrgDb column | rowData column        | default? |
|---|---|---|
| `SYMBOL`     | `<feature_type>_name` | yes |
| `GENENAME`   | `description`         | yes |
| `ENSEMBL`    | `ensembl_id`         | no |
| `ENTREZID`   | `entrez_id`          | no |
| `GENETYPE`   | `gene_biotype`       | no |

Chromosome is **not** taken from OrgDb (its `CHR` accessor is deprecated) — it comes from
the GTF (`seqnames`→`chromosome`). The default column set `c("SYMBOL","GENENAME")`
reproduces the original name + description behaviour.

## GTF attribute import (`annotate_with_gtf`)

Distinct from GTF length compute (that's `gtf_feature_lengths`, documented in `rnaseq-bioc`).
- **Match column** (`.resolve_match_col`): `"auto"` → Ensembl ids use `gene_id`, else
  `gene_name`, else `gene_id`; or the user picks an explicit column.
- **Attribute targets**: `gene_name`→`<feature_type>_name`, `seqnames`→`chromosome`, every
  other imported column keeps its own name.
- **GTF is authoritative over OrgDb** on matched rows (apply GTF *after* OrgDb).
- Build the one-row-per-feature table with **`gtf_attribute_table()`** (columnar, dedupes by
  the group key) — never `as.data.frame(gtf)` on the whole object (memory spike).

## Overlay semantics — fill, never wipe (`.fill`)

Both sources overlay onto existing `rowData`: write the mapped/matched value where the
source has a hit, **keep the prior value everywhere else**. Features absent from the source,
and extra source entries not in the `dds`, are simply ignored. So annotation only ever
*adds/refines*; it never blanks a column. `feature_class` and the raw `counts` assay are
never touched.

## Match flags & coverage

- Optional logical flag columns: **`in_orgdb`** / **`in_gtf`** (the `matched_col` argument),
  letting users filter unmatched rows. OrgDb's flag tracks the primary column (`SYMBOL` when
  requested, else the first requested column).
- **Coverage** for the UI banners: `orgdb_match_count()` / `gtf_match_count()` return
  `list(matched, total)` over *all* features (not just the preview). The shared
  `.coverage_banner()` colours it red (0) / amber (partial) / green (100%).
- **Previews** (`orgdb_annotation_preview`, the reader's `gtf_preview`) map only the first ~100
  ids — id (join key) first, then the columns that would be imported. They never commit.

## Draft-composition wiring (do not bypass)

Annotation edits the **editor draft**, not `state$working` (see `shiny-module` → composable
sub-modules). `mod_feature` calls `editor$set(annotate_*(editor$draft(), …))`; the user
commits with one `state_mutate` on **Save**. Writing straight to `state$working` would
re-seed the draft and silently drop unsaved edits — a real bug we hit and fixed. An
overwrite-confirm modal warns before replacing populated target columns (with a
"don't warn again this session" opt-out).

## Review checklist

- [ ] Ensembl keys version-stripped on both sides before matching.
- [ ] OrgDb columns land in their fixed `rowData` targets; chromosome only from GTF.
- [ ] GTF applied after OrgDb (authoritative); attribute import vs length compute kept separate.
- [ ] Overlay fills where matched and never wipes existing values; `counts`/`feature_class` untouched.
- [ ] Goes through the editor draft + Save (one `state_mutate`), not `state$working` directly.
- [ ] No whole-GTF `as.data.frame()`; uses `gtf_attribute_table()` / first-rows previews.
