# ERCC spike-in reference concentrations

`ercc_concentrations.csv` holds the nominal input concentrations of the 92 ERCC
spike-in control transcripts (External RNA Controls Consortium / Thermo Fisher
ERCC RNA Spike-In Mixes, cat. 4456740 / 4456739), used for the spike-in
dose-response QC (`R/ercc_helpers.R`).

## Columns
- `ercc_id`   — ERCC transcript id (e.g. `ERCC-00002`), matches `^ERCC-` rownames.
- `subgroup`  — ExFold subgroup `A`/`B`/`C`/`D` (the four expected Mix1:Mix2 ratios).
- `conc_mix1` — nominal concentration in **Mix 1** (attomoles/µL).
- `conc_mix2` — nominal concentration in **Mix 2** (attomoles/µL).

## Provenance
Derived verbatim from the official Thermo Fisher **"ERCC Controls Analysis"**
table (`cms_095046.txt`), columns *ERCC ID*, *subgroup*, *concentration in Mix 1
(attomoles/ul)*, *concentration in Mix 2 (attomoles/ul)*. The `Re-sort ID`,
`expected fold-change ratio`, and `log2(Mix 1/Mix 2)` columns are omitted (the
last two are derivable from the two concentrations). Values are unmodified.

Source URL: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt
