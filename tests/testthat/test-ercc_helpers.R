test_that("ercc_concentrations() returns the 92-transcript reference", {
  ref <- ercc_concentrations()
  expect_setequal(colnames(ref), c("ercc_id", "subgroup", "conc_mix1", "conc_mix2"))
  expect_equal(nrow(ref), 92L)
  expect_true(all(ref$conc_mix1 > 0) && all(ref$conc_mix2 > 0))
  expect_setequal(unique(ref$subgroup), c("A", "B", "C", "D"))
  expect_true(all(grepl("^ERCC-", ref$ercc_id)))
})

test_that("resolve_spike_concentration reads a column and each bundled mix", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 4, seed = 1)
  spike <- rownames(dds)[ddsdashboard:::.detect_spike_features(dds)]

  # From the rowData column the mock populates.
  v_col <- resolve_spike_concentration(dds, "column")
  expect_setequal(names(v_col), spike)
  expect_true(all(is.finite(v_col)))

  # From the bundled reference: mock ids ERCC-0000N exist in the 92-id table.
  v1 <- resolve_spike_concentration(dds, "mix1")
  v2 <- resolve_spike_concentration(dds, "mix2")
  ref <- ercc_concentrations()
  matched <- spike %in% ref$ercc_id
  expect_equal(unname(v1[matched]), ref$conc_mix1[match(spike[matched], ref$ercc_id)])
  expect_true(all(is.na(v1[!matched])))            # unmatched ids -> NA
  expect_false(isTRUE(all.equal(unname(v1), unname(v2))))  # mixes differ
})

test_that("spike_dose_response recovers slope ~ 1 / high R^2 on a clean titration", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 12, seed = 2)
  spike <- rownames(dds)[ddsdashboard:::.detect_spike_features(dds)]
  conc <- resolve_spike_concentration(dds, "column")
  # Make CPM-scale counts exactly proportional to concentration (clean dose-response).
  cnt <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  for (s in seq_len(ncol(cnt))) cnt[spike, s] <- round(conc[spike] * 5)
  SummarizedExperiment::assay(dds, "counts") <- cnt

  dr <- spike_dose_response(dds, assay = "CPM", source = "column")
  expect_setequal(colnames(dr$per_sample),
                  c("sample", "pct_spike", "n_spike_detected", "n_points",
                    "slope", "r_squared", "lod"))
  expect_equal(nrow(dr$per_sample), ncol(dds))
  ok <- !is.na(dr$per_sample$slope)
  expect_true(all(ok))
  expect_true(all(abs(dr$per_sample$slope[ok] - 1) < 0.1))      # ~ proportional -> slope ~ 1
  expect_true(all(dr$per_sample$r_squared[ok] > 0.99))
  expect_equal(unname(dr$per_sample$lod), rep(min(conc[spike]), ncol(dds)))
})

test_that("spike_dose_response drops zeros and is NA-safe below 3 points", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 5, seed = 3)
  spike <- rownames(dds)[ddsdashboard:::.detect_spike_features(dds)]
  cnt <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  cnt[spike, ] <- 0L; cnt[spike[1:2], ] <- 100L          # only 2 detected -> < 3 points
  SummarizedExperiment::assay(dds, "counts") <- cnt
  dr <- spike_dose_response(dds, source = "column")
  expect_true(all(is.na(dr$per_sample$slope)))           # < 3 usable points
  expect_true(all(dr$per_sample$n_spike_detected == 2L))
})

test_that("spike_dose_response is graceful with no spike-in features", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 30, n_per_group = 2, n_spike = 0, seed = 4)
  dr <- spike_dose_response(dds, source = "column")
  expect_equal(nrow(dr$long), 0L)
  expect_equal(nrow(dr$per_sample), ncol(dds))
  expect_true(all(dr$per_sample$n_spike_detected == 0L))
})

test_that("set_spike_concentration writes spike_concentration; missing checker flags gaps", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 4, seed = 5)
  spike <- rownames(dds)[ddsdashboard:::.detect_spike_features(dds)]
  dds <- add_meta_column(dds, "rowData", "my_conc", "numeric", NA)
  rd <- SummarizedExperiment::rowData(dds)
  rd$my_conc[match(spike, rownames(dds))] <- c(10, 20, 0, NA)   # one zero, one NA
  SummarizedExperiment::rowData(dds) <- rd
  out <- set_spike_concentration(dds, "my_conc")
  sc <- SummarizedExperiment::rowData(out)$spike_concentration
  expect_equal(sc[match(spike[1], rownames(out))], 10)
  expect_true(all(is.na(sc[!ddsdashboard:::.detect_spike_features(out)])))  # non-spike -> NA
  miss <- spike_features_missing_conc(out)
  expect_setequal(miss, spike[3:4])                                         # the 0 and NA
})
