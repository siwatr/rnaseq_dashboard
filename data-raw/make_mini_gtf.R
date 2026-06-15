# Generate inst/extdata/demo_annotation.gtf: a tiny GTF fixture used by the GTF
# annotation tests and as a try-it demo. Three Ensembl-style genes with gene /
# transcript / exon rows. GeneA's exons overlap and are split across the locus,
# so union-exon length (reduce -> sum width) differs from the gene span:
#
#   union-exon length (type = "exon"):  GeneA 202, GeneB 201, GeneC 501
#   gene span        (type = "gene"):   GeneA 351, GeneB 201, GeneC 501
#
# (GTF is 1-based inclusive, so width = end - start + 1.) Hand-checkable values
# are asserted in tests/testthat/test-gtf_helpers.R. Pure text writer so the
# fixture can be regenerated without rtracklayer.

attrs <- function(...) paste0(paste0(sprintf('%s "%s";', c(...)[c(TRUE, FALSE)],
                                             c(...)[c(FALSE, TRUE)]), collapse = " "))

row <- function(seqname, type, start, end, strand, ...) {
  paste(seqname, "mini", type, start, end, ".", strand, ".", attrs(...), sep = "\t")
}

lines <- c(
  # GeneA (chr1, +): exons 100-200, 150-250 (overlap -> 100-250), 400-450
  row("chr1", "gene",       100, 450, "+", "gene_id", "ENSG00000000001", "gene_name", "GeneA", "gene_biotype", "protein_coding"),
  row("chr1", "transcript", 100, 450, "+", "gene_id", "ENSG00000000001", "transcript_id", "ENST00000000001", "gene_name", "GeneA", "gene_biotype", "protein_coding"),
  row("chr1", "exon",       100, 200, "+", "gene_id", "ENSG00000000001", "transcript_id", "ENST00000000001", "gene_name", "GeneA", "gene_biotype", "protein_coding"),
  row("chr1", "exon",       150, 250, "+", "gene_id", "ENSG00000000001", "transcript_id", "ENST00000000001", "gene_name", "GeneA", "gene_biotype", "protein_coding"),
  row("chr1", "exon",       400, 450, "+", "gene_id", "ENSG00000000001", "transcript_id", "ENST00000000001", "gene_name", "GeneA", "gene_biotype", "protein_coding"),
  # GeneB (chr2, -): exons 1000-1100, 1100-1200 (touch -> 1000-1200)
  row("chr2", "gene",       1000, 1200, "-", "gene_id", "ENSG00000000002", "gene_name", "GeneB", "gene_biotype", "protein_coding"),
  row("chr2", "transcript", 1000, 1200, "-", "gene_id", "ENSG00000000002", "transcript_id", "ENST00000000002", "gene_name", "GeneB", "gene_biotype", "protein_coding"),
  row("chr2", "exon",       1000, 1100, "-", "gene_id", "ENSG00000000002", "transcript_id", "ENST00000000002", "gene_name", "GeneB", "gene_biotype", "protein_coding"),
  row("chr2", "exon",       1100, 1200, "-", "gene_id", "ENSG00000000002", "transcript_id", "ENST00000000002", "gene_name", "GeneB", "gene_biotype", "protein_coding"),
  # GeneC (chr3, +): single exon
  row("chr3", "gene",       5000, 5500, "+", "gene_id", "ENSG00000000003", "gene_name", "GeneC", "gene_biotype", "lincRNA"),
  row("chr3", "transcript", 5000, 5500, "+", "gene_id", "ENSG00000000003", "transcript_id", "ENST00000000003", "gene_name", "GeneC", "gene_biotype", "lincRNA"),
  row("chr3", "exon",       5000, 5500, "+", "gene_id", "ENSG00000000003", "transcript_id", "ENST00000000003", "gene_name", "GeneC", "gene_biotype", "lincRNA")
)

out <- file.path("inst", "extdata", "demo_annotation.gtf")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, out)
message("Wrote ", out)
