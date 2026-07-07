library(testthat)

.defaultQuantifySpecPath <- function() {
  system.file("extdata", "quantify_default_spec.tsv", package = "Rdiscvr")
}

.makeSyntheticCellTable <- function() {
  tibble::tibble(
    sourceOutputFileId = rep(111L, 30),
    cDNA_ID = rep(c(1001L, 1002L), each = 15),
    CellType = c(
      rep(c("CD8_T", "NK", "B"), length.out = 15),
      rep(c("CD14_Mono", "CD8_T", "B"), length.out = 15)
    ),
    ScoreCol = c(
      rep(c(0.1, 0.5, 0.9), length.out = 15),
      rep(c(0.2, 0.6, 1.0), length.out = 15)
    ),
    ImmuneClass = c(
      rep(c("T_NK", "T_NK", "Myeloid"), length.out = 15),
      rep(c("Myeloid", "T_NK", "Myeloid"), length.out = 15)
    ),
    TNKClass = c(
      rep(c("CD8_T", "CD4_T", "CD8_T"), length.out = 15),
      rep(c("CD8_T", "CD8_T", "CD4_T"), length.out = 15)
    )
  )
}

.makeCompactRiraCellTable <- function() {
  tibble::tibble(
    sourceOutputFileId = rep(111L, 80),
    cDNA_ID = c(
      rep(1001L, 20), rep(1001L, 8), rep(1001L, 10),
      rep(1002L, 15), rep(1002L, 12), rep(1002L, 10), rep(1002L, 5)
    ),
    `RIRA_Immune_v2.cellclass` = c(
      rep("T_NK", 20), rep("Myeloid", 8), rep("Bcell", 10),
      rep("T_NK", 15), rep("Myeloid", 12), rep("Bcell", 10), rep("Epithelial", 5)
    ),
    `RIRA_TNK_v2.cellclass` = c(
      rep(c("CD8+ T Cells", "NK Cells"), c(8, 6)),
      rep(c("CD4+ T Cells", "Other"), c(4, 2)),
      rep("Unassigned", 8),
      rep("Unassigned", 10),
      rep(c("CD8+ T Cells", "NK Cells"), c(10, 5)),
      rep("Unassigned", 12),
      rep("Unassigned", 10),
      rep("Unassigned", 5)
    ),
    `RIRA_Myeloid_v3.cellclass` = c(
      rep("Unassigned", 20),
      rep("CD14+ Monocytes", 8),
      rep("Unassigned", 10),
      rep("Unassigned", 15),
      rep("CD16+ Monocytes", 12),
      rep("Unassigned", 10),
      rep("Unassigned", 5)
    ),
    `RIRA_Myeloid_v3.coarseclass` = c(
      rep("Unassigned", 20),
      rep("Monocytes", 8),
      rep("Unassigned", 10),
      rep("Unassigned", 15),
      rep("Monocytes", 12),
      rep("Unassigned", 10),
      rep("Unassigned", 5)
    ),
    Cytotoxicity_UCell = c(
      rep(c(0.3, 0.5, 0.7, 0.9), length.out = 20),
      rep(0.25, 8), rep(0.5, 10),
      rep(0.5, 15), rep(0.4, 12), rep(0.5, 10), rep(0.5, 5)
    ),
    Interferon_Response_UCell = rep(0.5, 80),
    TNK_Type = c(
      rep("Alpha/Beta", 14), rep("NK (CD3-/TCR-)", 6),
      rep("Unassigned", 8), rep("Unassigned", 10),
      rep("Alpha/Beta", 10), rep("NK (CD3-/TCR-)", 5),
      rep("Unassigned", 12), rep("Unassigned", 10), rep("Unassigned", 5)
    ),
    TRD = rep(NA_character_, 80),
    HasCDR3Data = c(
      rep(TRUE, 14), rep(FALSE, 6), rep(TRUE, 18),
      rep(TRUE, 10), rep(FALSE, 5), rep(TRUE, 27)
    ),
    HasCD3 = c(
      rep(TRUE, 14), rep(FALSE, 6), rep(TRUE, 18),
      rep(TRUE, 10), rep(FALSE, 5), rep(TRUE, 27)
    ),
    Tcell_EffectorDifferentiation = rep(4, 80)
  )
}

.mockDownload <- function(synthetic) {
  function(output_file_id) {
    lane <- synthetic[synthetic$sourceOutputFileId == output_file_id, , drop = FALSE]
    lane$sourceOutputFileId <- output_file_id
    as.data.frame(lane)
  }
}

.runQuantifyWithMock <- function(
  synthetic,
  spec,
  failure_log_file = tempfile(fileext = ".txt"),
  classify_tnk = TRUE
) {
  with_mocked_bindings(
    `.DownloadCellMetadataPerLane` = .mockDownload(synthetic),
    .package = "Rdiscvr",
    {
      Quantify10xData(
        outputFileIds = unique(synthetic$sourceOutputFileId),
        spec = spec,
        failureLogFile = failure_log_file,
        classifyTNK = classify_tnk
      )
    }
  )
}

test_that("Quantify10xData quantifies sum, score, and scope rules", {
  synthetic <- .makeSyntheticCellTable()

  spec <- tibble::tibble(
    GroupingVariable = rep("cDNA_ID", 7),
    TargetField = c(
      "CD8_count", "B_count", "CD8_score",
      "TNK_CD8_count", "TNK_CD8_score",
      "missing_count", "missing_score"
    ),
    SourceField = c(
      "CellType", "CellType", "CellType",
      "TNKClass", "TNKClass",
      "CellType", "CellType"
    ),
    MatchValue = c(
      "CD8_T", "B", "CD8_T",
      "CD8_T", "CD8_T",
      "NotPresent", "NotPresent"
    ),
    QuantificationType = c("sum", "sum", "score", "sum", "score", "sum", "score"),
    QuantificationSourceField = c("", "", "ScoreCol", "", "ScoreCol", "", "ScoreCol"),
    QuantificationScoreType = c("", "", "quantiles", "", "quantiles", "", "quantiles"),
    ScopeField = c("", "", "", "ImmuneClass", "ImmuneClass", "", ""),
    ScopeMatchValue = c("", "", "", "T_NK", "T_NK", "", "")
  )

  result <- .runQuantifyWithMock(synthetic, spec)

  # check shape
  expect_true(all(c("countsWide", "failures", "laneSummary") %in% names(result)))
  expect_equal(nrow(result$countsWide), 2)
  expect_equal(result$laneSummary$nCellsIngested, 30L)
  expect_equal(nrow(result$failures), 0)

  # sum per cDNA_ID
  expect_equal(result$countsWide$CD8_count[result$countsWide$cDNA_ID == 1001], 5L)
  expect_equal(result$countsWide$B_count[result$countsWide$cDNA_ID == 1002], 5L)

  # score quantiles
  cd8_scores <- synthetic$ScoreCol[synthetic$CellType == "CD8_T" & synthetic$cDNA_ID == 1001]
  expect_true(all(c("CD8_score__p05", "CD8_score__median", "CD8_score__p95") %in% names(result$countsWide)))
  expect_equal(
    result$countsWide$CD8_score__median[result$countsWide$cDNA_ID == 1001],
    as.numeric(stats::quantile(cd8_scores, probs = 0.5, na.rm = TRUE, names = FALSE))
  )

  # ScopeField restricts matching to parent class
  expect_equal(result$countsWide$TNK_CD8_count[result$countsWide$cDNA_ID == 1001], 5L)
  expect_equal(
    result$countsWide$TNK_CD8_score__median[result$countsWide$cDNA_ID == 1001],
    as.numeric(stats::quantile(rep(0.1, 5), probs = 0.5, na.rm = TRUE, names = FALSE))
  )

  # zero-match sum is 0 and zero-match quantiles are NA
  expect_equal(result$countsWide$missing_count, rep(0L, 2))
  expect_true(all(is.na(result$countsWide$missing_score__median)))
})

test_that("Default RIRA spec integrates with synthetic metadata", {
  spec_path <- .defaultQuantifySpecPath()
  expect_true(nzchar(spec_path) && file.exists(spec_path))

  synthetic <- .makeCompactRiraCellTable()
  result <- .runQuantifyWithMock(synthetic, spec_path)

  expect_equal(nrow(result$failures), 0)
  expect_equal(nrow(result$countsWide), 2)

  row_1001 <- result$countsWide[result$countsWide$cDNA_ID == 1001, ]
  row_1002 <- result$countsWide[result$countsWide$cDNA_ID == 1002, ]

  expect_equal(row_1001$Immune__T_NK, 20L)
  expect_equal(row_1001$TNK__CD8plus_T_Cells, 8L)
  expect_equal(row_1001$Myeloid__coarse__Monocytes, 8L)
  expect_false(is.na(row_1001$Cytotoxic__TNK__CD8plus_T_Cells__median))

  expect_equal(row_1002$Immune__T_NK, 15L)
  expect_equal(row_1002$TNK__CD8plus_T_Cells, 10L)
  expect_equal(row_1002$Myeloid__coarse__Monocytes, 12L)
})

test_that("TCR TNK override and failure tolerance", {
  nk_synthetic <- tibble::tibble(
    sourceOutputFileId = rep(111L, 5),
    cDNA_ID = rep(1001L, 5),
    `RIRA_Immune_v2.cellclass` = rep("T_NK", 5),
    `RIRA_TNK_v2.cellclass` = rep("CD8+ T Cells", 5),
    `RIRA_Myeloid_v3.cellclass` = rep("Unassigned", 5),
    `RIRA_Myeloid_v3.coarseclass` = rep("Unassigned", 5),
    Cytotoxicity_UCell = rep(0.5, 5),
    Interferon_Response_UCell = rep(0.5, 5),
    TNK_Type = c(
      "NK (CD3-/TCR-)", "NK (CD3-/TCR-)", "Alpha/Beta",
      "NK (CD3-/TCR-)", "Alpha/Beta"
    ),
    TRD = rep(NA_character_, 5),
    HasCDR3Data = c(FALSE, FALSE, TRUE, FALSE, TRUE),
    HasCD3 = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  )

  nk_spec <- tibble::tibble(
    GroupingVariable = "cDNA_ID",
    TargetField = "TNK__NK_Cells",
    SourceField = "RIRA_TNK_v2.cellclass",
    MatchValue = "NK Cells",
    QuantificationType = "sum",
    QuantificationSourceField = "",
    QuantificationScoreType = "",
    ScopeField = "RIRA_Immune_v2.cellclass",
    ScopeMatchValue = "T_NK"
  )

  # NK cells: TCR evidence overrides RIRA CD8 label when classifyTNK=TRUE
  nk_result <- .runQuantifyWithMock(nk_synthetic, nk_spec, classify_tnk = TRUE)
  expect_equal(nk_result$countsWide$TNK__NK_Cells, 3L)

  gd_synthetic <- tibble::tibble(
    sourceOutputFileId = rep(111L, 4),
    cDNA_ID = rep(1001L, 4),
    `RIRA_Immune_v2.cellclass` = rep("T_NK", 4),
    `RIRA_TNK_v2.cellclass` = c(
      "CD8+ T Cells", "CD4+ T Cells", "Unassigned", "Unassigned"
    ),
    `RIRA_Myeloid_v3.cellclass` = rep("Unassigned", 4),
    `RIRA_Myeloid_v3.coarseclass` = rep("Unassigned", 4),
    Cytotoxicity_UCell = rep(0.5, 4),
    Interferon_Response_UCell = rep(0.5, 4),
    TNK_Type = c("Gamma/Delta", "Gamma Chain-Only", "Gamma/Delta", "Alpha/Beta"),
    TRD = c("CALWELQYAGNVLTF", NA_character_, "CALWEVSNADKLTF", NA_character_),
    HasCDR3Data = rep(TRUE, 4),
    HasCD3 = rep(TRUE, 4)
  )

  gd_spec <- tibble::tibble(
    GroupingVariable = "cDNA_ID",
    TargetField = "TNK__Gamma_Delta_Cells",
    SourceField = "RIRA_TNK_v2.cellclass",
    MatchValue = "Gamma/Delta Cells",
    QuantificationType = "sum",
    QuantificationSourceField = "",
    QuantificationScoreType = "",
    ScopeField = "RIRA_Immune_v2.cellclass",
    ScopeMatchValue = "T_NK"
  )

  # Gamma/Delta: counts Gamma/Delta, not Gamma Chain-Only
  gd_result <- .runQuantifyWithMock(gd_synthetic, gd_spec, classify_tnk = TRUE)
  expect_equal(gd_result$countsWide$TNK__Gamma_Delta_Cells, 2L)

  # classifyTNK=FALSE: RIRA labels preserved (NK count stays 0)
  rira_only <- .runQuantifyWithMock(nk_synthetic, nk_spec, classify_tnk = FALSE)
  expect_equal(rira_only$countsWide$TNK__NK_Cells, 0L)

  # Bad spec row: missing metadata column logs warning, records failure, continues
  bad_spec <- tibble::tibble(
    GroupingVariable = "cDNA_ID",
    TargetField = "bad_count",
    SourceField = "MissingColumn",
    MatchValue = "X",
    QuantificationType = "sum",
    QuantificationSourceField = "",
    QuantificationScoreType = ""
  )
  failure_log <- tempfile(fileext = ".txt")
  expect_warning(
    bad_result <- .runQuantifyWithMock(nk_synthetic, bad_spec, failure_log_file = failure_log),
    "quantification failure"
  )
  expect_false("bad_count" %in% names(bad_result$countsWide))
  expect_true(any(bad_result$failures$field == "SourceField"))
  expect_gt(length(readLines(failure_log)), 0)
})

test_that("EDS subset columns split T cells by score cutpoints", {
  eds_synthetic <- tibble::tibble(
    sourceOutputFileId = rep(111L, 3),
    cDNA_ID = rep(1001L, 3),
    `RIRA_Immune_v2.cellclass` = rep("T_NK", 3),
    `RIRA_TNK_v2.cellclass` = rep("CD4+ T Cells", 3),
    `RIRA_Myeloid_v3.cellclass` = rep("Unassigned", 3),
    `RIRA_Myeloid_v3.coarseclass` = rep("Unassigned", 3),
    Tcell_EffectorDifferentiation = c(1, 4, 8),
    HasCDR3Data = rep(TRUE, 3),
    HasCD3 = rep(TRUE, 3)
  )

  eds_spec <- tibble::tibble(
    GroupingVariable = rep("cDNA_ID", 3),
    TargetField = c(
      "TNK__CD4plus_T_Cells__Naive",
      "TNK__CD4plus_T_Cells__MemoryLike",
      "TNK__CD4plus_T_Cells__Effector"
    ),
    SourceField = rep("RIRA_TNK_v2.cellclass", 3),
    MatchValue = rep("CD4+ T Cells", 3),
    QuantificationType = rep("sum", 3),
    QuantificationSourceField = rep("", 3),
    QuantificationScoreType = rep("", 3),
    ScopeField = rep("RIRA_Immune_v2.cellclass", 3),
    ScopeMatchValue = rep("T_NK", 3),
    EffectorDifferentiationScoreField = rep("Tcell_EffectorDifferentiation", 3),
    SubsetCutpointLow = rep("2", 3),
    SubsetCutpointHigh = rep("6", 3),
    SubsetPhenotype = c("Naive", "MemoryLike", "Effector")
  )

  eds_result <- .runQuantifyWithMock(eds_synthetic, eds_spec)
  expect_equal(nrow(eds_result$failures), 0)
  expect_equal(eds_result$countsWide$TNK__CD4plus_T_Cells__Naive, 1L)
  expect_equal(eds_result$countsWide$TNK__CD4plus_T_Cells__MemoryLike, 1L)
  expect_equal(eds_result$countsWide$TNK__CD4plus_T_Cells__Effector, 1L)

  missing_eds_spec <- tibble::tibble(
    GroupingVariable = "cDNA_ID",
    TargetField = "TNK__CD4plus_T_Cells__Naive",
    SourceField = "RIRA_TNK_v2.cellclass",
    MatchValue = "CD4+ T Cells",
    QuantificationType = "sum",
    QuantificationSourceField = "",
    QuantificationScoreType = "",
    ScopeField = "RIRA_Immune_v2.cellclass",
    ScopeMatchValue = "T_NK",
    EffectorDifferentiationScoreField = "Tcell_EffectorDifferentiation",
    SubsetCutpointLow = "2",
    SubsetCutpointHigh = "6",
    SubsetPhenotype = "Naive"
  )
  eds_no_col <- eds_synthetic[, setdiff(names(eds_synthetic), "Tcell_EffectorDifferentiation")]
  failure_log <- tempfile(fileext = ".txt")
  expect_warning(
    missing_result <- .runQuantifyWithMock(
      eds_no_col,
      missing_eds_spec,
      failure_log_file = failure_log
    ),
    "quantification failure"
  )
  expect_true(any(missing_result$failures$field == "EffectorDifferentiationScoreField"))
  expect_false("TNK__CD4plus_T_Cells__Naive" %in% names(missing_result$countsWide))
})
