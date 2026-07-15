library(testthat)

.rhesusQuantifySpecPath <- function() {
  system.file("extdata", "quantify_rhesus_spec.tsv", package = "Rdiscvr")
}

.humanImmuneQuantifySpecPath <- function() {
  system.file("extdata", "quantify_human_immune_spec.tsv", package = "Rdiscvr")
}

# Metadata shaped like a RIRA-classified rhesus lane for quantify_rhesus_spec.tsv
.makeRhesusRiraCellTable <- function() {
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
    Tcell_EffectorDifferentiation = rep(4, 80),
    TandNK_Activation2_UCell = rep(0.2, 80),
    MHCII_UCell = rep(0.1, 80),
    Proliferation_UCell = rep(0.05, 80),
    Perforin_UCell = rep(0.15, 80),
    Is_TCR_Stimulated = rep(c(TRUE, FALSE), length.out = 80),
    TRA = ifelse(seq_len(80) %% 3 == 0, NA_character_, paste0("TRA", seq_len(80) %% 10)),
    TRB = ifelse(seq_len(80) %% 3 == 0, NA_character_, paste0("TRB", seq_len(80) %% 8)),
    TRA_V = ifelse(seq_len(80) %% 3 == 0, NA_character_, "TRAV1"),
    TRB_V = ifelse(seq_len(80) %% 3 == 0, NA_character_, "TRBV1"),
    PDCD1 = as.integer(seq_len(80) %% 4 == 0),
    KLRK1 = as.integer(seq_len(80) %% 5 == 0),
    KLRB1 = as.integer(seq_len(80) %% 6 == 0),
    HAVCR2 = as.integer(seq_len(80) %% 7 == 0),
    TIGIT = as.integer(seq_len(80) %% 3 == 0),
    PRF1 = as.integer(seq_len(80) %% 4 == 1),
    GZMB = as.integer(seq_len(80) %% 5 == 1),
    KLRC1 = as.integer(seq_len(80) %% 6 == 1),
    FCGR3 = as.integer(seq_len(80) %% 7 == 1),
    FCGR3A = as.integer(seq_len(80) %% 8 == 1),
    FOXP3 = as.integer(seq_len(80) %% 9 == 1),
    IL2RA = as.integer(seq_len(80) %% 10 == 1)
  )
}

# minimal metadata for QuantificationType (sum/score/pct_positive) and ScopeField
.makeTypedQuantifyCellTable <- function() {
  tibble::tibble(
    sourceOutputFileId = rep(111L, 9),
    cDNA_ID = c(rep(1001L, 6), rep(1002L, 3)),
    CellType = c(rep("CD8_T", 3), rep("B", 3), rep("B", 3)),
    ScoreCol = c(0.1, 0.5, 0.9, 0.2, 0.4, 0.6, 0.3, 0.5, 0.7),
    ImmuneClass = c(rep("T_NK", 3), rep("Myeloid", 3), rep("Myeloid", 3)),
    TNKClass = c("CD8_T", "CD8_T", "CD4_T", "CD8_T", "CD8_T", "CD4_T", "CD8_T", "CD8_T", "CD4_T"),
    PDCD1 = c(1, 0, 1, 0, 0, 0, 0, 0, 0)
  )
}

.mockPrepareLaneCellTable <- function(cell_table) {
  function(
    output_file_id,
    spec_table,
    classify_tnk,
    apply_role_map = FALSE,
    role_map = NULL,
    rhesus_spec_table = NULL
  ) {
    lane <- cell_table[cell_table$sourceOutputFileId == output_file_id, , drop = FALSE]
    lane$sourceOutputFileId <- output_file_id
    lane <- tibble::as_tibble(lane)

    detected_species <- if ("celltypist.Immune_All_High.cellclass" %in% names(lane)) {
      "human"
    } else if ("RIRA_Immune_v2.cellclass" %in% names(lane)) {
      "rhesus"
    } else {
      NULL
    }

    lane_apply_map <- isTRUE(apply_role_map) && identical(detected_species, "human")
    if (lane_apply_map && !is.null(role_map) && nrow(role_map)) {
      lane <- Rdiscvr:::.ApplyHumanToRhesusRoleMap(lane, role_map)
    }

    lane_classify_tnk <- isTRUE(classify_tnk) && identical(detected_species, "rhesus")
    if (lane_classify_tnk) {
      override_result <- Rdiscvr:::.ApplyTcrTNKOverridesToCellTable(lane)
      lane <- override_result$cell_table
      failures <- override_result$failures
    } else {
      failures <- Rdiscvr:::.EmptyFailuresTibble()
    }

    list(
      cell_table = lane,
      failures = failures,
      detected_species = detected_species
    )
  }
}

.quantifyWithMockedPrepare <- function(
  cell_table,
  spec,
  failure_log_file = tempfile(fileext = ".txt"),
  classify_tnk = TRUE,
  species = NULL,
  coerceToRIRA = FALSE,
  outputPrefix = NULL,
  apply_role_map = coerceToRIRA,
  role_map = if (apply_role_map) LoadQuantifyRoleMap() else NULL
) {
  with_mocked_bindings(
    `.PrepareLaneCellTable` = .mockPrepareLaneCellTable(cell_table),
    .package = "Rdiscvr",
    {
      Quantify10xData(
        outputFileIds = unique(cell_table$sourceOutputFileId),
        spec = spec,
        failureLogFile = failure_log_file,
        classifyTNK = classify_tnk,
        species = species,
        coerceToRIRA = coerceToRIRA,
        outputPrefix = outputPrefix
      )
    }
  )
}

test_that("Quantify10xData applies QuantificationType and ScopeField rules", {
  cell_table <- .makeTypedQuantifyCellTable()
  spec <- tibble::tibble(
    GroupingVariable = rep("cDNA_ID", 8),
    TargetField = c(
      "CD8_count", "B_count", "CD8_score",
      "TNK_CD8_count", "TNK_CD8_score",
      "missing_count", "missing_score",
      "CD8_PDCD1_pct"
    ),
    SourceField = c(
      "CellType", "CellType", "CellType",
      "TNKClass", "TNKClass",
      "CellType", "CellType",
      "CellType"
    ),
    MatchValue = c(
      "CD8_T", "B", "CD8_T",
      "CD8_T", "CD8_T",
      "NotPresent", "NotPresent",
      "CD8_T"
    ),
    QuantificationType = c(
      "sum", "sum", "score", "sum", "score", "sum", "score", "pct_positive"
    ),
    QuantificationSourceField = c(
      "", "", "ScoreCol", "", "ScoreCol", "", "ScoreCol", "PDCD1"
    ),
    QuantificationScoreType = c(
      "", "", "quantiles", "", "quantiles", "", "quantiles", ""
    ),
    ScopeField = c("", "", "", "ImmuneClass", "ImmuneClass", "", "", ""),
    ScopeMatchValue = c("", "", "", "T_NK", "T_NK", "", "", "")
  )

  result <- .quantifyWithMockedPrepare(cell_table, spec)
  expect_equal(nrow(result$failures), 0)
  expect_equal(result$countsWide$CD8_count[result$countsWide$cDNA_ID == 1001], 3L)
  expect_equal(result$countsWide$B_count[result$countsWide$cDNA_ID == 1002], 3L)
  expect_equal(result$countsWide$CD8_score__median[result$countsWide$cDNA_ID == 1001], 0.5)
  expect_equal(result$countsWide$TNK_CD8_count[result$countsWide$cDNA_ID == 1001], 2L)
  expect_equal(result$countsWide$missing_count, rep(0L, 2))
  expect_true(all(is.na(result$countsWide$missing_score__median)))
  expect_equal(
    result$countsWide$CD8_PDCD1_pct[result$countsWide$cDNA_ID == 1001],
    100 * (2 / 3)
  )
  # groups with no matched cells get 0 (same as sum), not NA from left_join
  expect_equal(
    result$countsWide$CD8_PDCD1_pct[result$countsWide$cDNA_ID == 1002],
    0
  )
})

test_that("Quantify10xData records SourceField failures and continues", {
  cell_table <- tibble::tibble(
    sourceOutputFileId = 111L,
    cDNA_ID = 1001L,
    CellType = "T"
  )
  spec <- tibble::tibble(
    GroupingVariable = rep("cDNA_ID", 2),
    TargetField = c("T_count", "bad_count"),
    SourceField = c("CellType", "MissingColumn"),
    MatchValue = c("T", "X"),
    QuantificationType = rep("sum", 2),
    QuantificationSourceField = rep("", 2),
    QuantificationScoreType = rep("", 2)
  )
  expect_warning(
    result <- .quantifyWithMockedPrepare(cell_table, spec),
    "quantification failure"
  )
  expect_equal(result$countsWide$T_count, 1L)
  expect_false("bad_count" %in% names(result$countsWide))
  expect_true(any(result$failures$field == "SourceField"))
})

test_that("classifyTNK overrides RIRA labels for NK and gamma-delta cells", {
  nk_cells <- tibble::tibble(
    sourceOutputFileId = rep(111L, 5),
    cDNA_ID = rep(1001L, 5),
    `RIRA_Immune_v2.cellclass` = rep("T_NK", 5),
    `RIRA_TNK_v2.cellclass` = rep("CD8+ T Cells", 5),
    TNK_Type = c(
      "NK (CD3-/TCR-)", "NK (CD3-/TCR-)", "Alpha/Beta",
      "NK (CD3-/TCR-)", "Alpha/Beta"
    ),
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
  expect_equal(
    .quantifyWithMockedPrepare(nk_cells, nk_spec, classify_tnk = TRUE)$countsWide$TNK__NK_Cells,
    3L
  )
  expect_equal(
    .quantifyWithMockedPrepare(nk_cells, nk_spec, classify_tnk = FALSE)$countsWide$TNK__NK_Cells,
    0L
  )

  gamma_delta_cells <- tibble::tibble(
    sourceOutputFileId = rep(111L, 4),
    cDNA_ID = rep(1001L, 4),
    `RIRA_Immune_v2.cellclass` = rep("T_NK", 4),
    `RIRA_TNK_v2.cellclass` = c(
      "CD8+ T Cells", "CD4+ T Cells", "Unassigned", "Unassigned"
    ),
    TNK_Type = c("Gamma/Delta", "Gamma Chain-Only", "Gamma/Delta", "Alpha/Beta"),
    HasCDR3Data = rep(TRUE, 4),
    HasCD3 = rep(TRUE, 4)
  )
  gamma_delta_spec <- tibble::tibble(
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
  expect_equal(
    .quantifyWithMockedPrepare(
      gamma_delta_cells,
      gamma_delta_spec,
      classify_tnk = TRUE
    )$countsWide$TNK__Gamma_Delta_Cells,
    2L
  )
})

test_that("EffectorDifferentiationScore cutpoints split Naive/MemoryLike/Effector", {
  cell_table <- tibble::tibble(
    sourceOutputFileId = rep(111L, 3),
    cDNA_ID = rep(1001L, 3),
    `RIRA_Immune_v2.cellclass` = rep("T_NK", 3),
    `RIRA_TNK_v2.cellclass` = rep("CD4+ T Cells", 3),
    Tcell_EffectorDifferentiation = c(1, 4, 8),
    HasCDR3Data = rep(TRUE, 3),
    HasCD3 = rep(TRUE, 3)
  )
  spec <- tibble::tibble(
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
    EffectorDifferentiationCutpointLow = rep("2", 3),
    EffectorDifferentiationCutpointHigh = rep("6", 3),
    SubsetPhenotypeOutputFieldName = c("Naive", "MemoryLike", "Effector")
  )
  result <- .quantifyWithMockedPrepare(cell_table, spec)
  expect_equal(result$countsWide$TNK__CD4plus_T_Cells__Naive, 1L)
  expect_equal(result$countsWide$TNK__CD4plus_T_Cells__MemoryLike, 1L)
  expect_equal(result$countsWide$TNK__CD4plus_T_Cells__Effector, 1L)
})

.makeHumanCellTypistCellTable <- function() {
  tibble::tibble(
    sourceOutputFileId = rep(222L, 12),
    cDNA_ID = rep(2001L, 12),
    `celltypist.Immune_All_High.cellclass` = c(
      rep("T cells", 6), rep("B cells", 3), rep("Monocytes", 3)
    ),
    `celltypist.Immune_All_Low.cellclass` = c(
      rep("Tcm/Naive helper T cells", 3),
      rep("Tem/Temra cytotoxic T cells", 3),
      rep("Naive B cells", 3),
      rep("Classical monocytes", 3)
    ),
    Cytotoxicity_UCell = seq(0.1, 0.9, length.out = 12),
    Interferon_Response_UCell = rep(0.4, 12),
    TandNK_Activation2_UCell = rep(0.2, 12),
    MHCII_UCell = rep(0.15, 12),
    Proliferation_UCell = rep(0.05, 12),
    Perforin_UCell = rep(0.12, 12),
    Tcell_EffectorDifferentiation = c(1, 4, 8, 2, 5, 7, rep(4, 6)),
    Is_TCR_Stimulated = rep(c(TRUE, FALSE), length.out = 12),
    TRA = ifelse(seq_len(12) %% 3 == 0, NA_character_, paste0("TRA", seq_len(12))),
    TRB = ifelse(seq_len(12) %% 3 == 0, NA_character_, paste0("TRB", seq_len(12))),
    PDCD1 = as.integer(seq_len(12) %% 2 == 0),
    KLRK1 = as.integer(seq_len(12) %% 3 == 0),
    KLRB1 = 0L,
    HAVCR2 = 0L,
    TIGIT = 0L,
    PRF1 = 0L,
    GZMB = as.integer(seq_len(12) %% 4 == 0),
    KLRC1 = 0L,
    FCGR3 = 0L,
    FCGR3A = 0L,
    FOXP3 = 0L,
    IL2RA = 0L
  )
}

test_that("ApplyHumanToRhesusRoleMap prefers Low over High labels", {
  cell_table <- tibble::tibble(
    `celltypist.Immune_All_High.cellclass` = c("T cells", "B cells"),
    `celltypist.Immune_All_Low.cellclass` = c(
      "Tcm/Naive helper T cells",
      "Naive B cells"
    )
  )
  mapped <- Rdiscvr:::.ApplyHumanToRhesusRoleMap(cell_table, LoadQuantifyRoleMap())
  expect_equal(mapped$`RIRA_Immune_v2.cellclass`, c("T_NK", "Bcell"))
  expect_equal(mapped$`RIRA_TNK_v2.cellclass`[1], "CD4+ T Cells")
  expect_equal(mapped$`RIRA_TNK_v2.cellclass`[2], "Unknown")
})

test_that("role map puts helper ILCs on TNK Other and NKT on CD8", {
  role_map <- LoadQuantifyRoleMap()
  cell_table <- tibble::tibble(
    `celltypist.Immune_All_High.cellclass` = c(
      "ILC", "ILC", "T cells", "ILC", "Cycling cells", "Cycling cells"
    ),
    `celltypist.Immune_All_Low.cellclass` = c(
      "ILC1", "ILC2", "NKT cells", "NK cells",
      "Cycling NK cells", "Cycling B cells"
    )
  )
  mapped <- Rdiscvr:::.ApplyHumanToRhesusRoleMap(cell_table, role_map)

  expect_equal(
    mapped$`RIRA_Immune_v2.cellclass`,
    c("T_NK", "T_NK", "T_NK", "T_NK", "T_NK", "Bcell")
  )
  expect_equal(
    mapped$`RIRA_TNK_v2.cellclass`,
    c("Other", "Other", "CD8+ T Cells", "NK Cells", "NK Cells", "Unknown")
  )

  ilc_roles <- role_map[
    role_map$humanSourceField == "Immune_All_Low" &
      role_map$humanLabel %in% c("ILC1", "ILC2", "ILC3", "ILC", "ILC precursor") &
      role_map$rhesusTargetField == "RIRA_TNK_v2.cellclass",
    ,
    drop = FALSE
  ]
  expect_true(nrow(ilc_roles) >= 5)
  expect_true(all(ilc_roles$rhesusLabel == "Other"))
  expect_true(all(ilc_roles$lineageRole == "NK_ILC"))
})

test_that("human coerceToRIRA returns rhesus-shaped countsWideRIRA columns", {
  human_path <- .humanImmuneQuantifySpecPath()
  rhesus_path <- .rhesusQuantifySpecPath()
  cell_table <- .makeHumanCellTypistCellTable()

  result <- .quantifyWithMockedPrepare(
    cell_table,
    human_path,
    species = "human",
    coerceToRIRA = TRUE
  )

  expect_true("countsWideRIRA" %in% names(result))
  expect_true(any(grepl("^ImmuneHigh__", names(result$countsWide))))
  expect_true("Immune__T_NK" %in% names(result$countsWideRIRA))
  expect_true("TNK__CD4plus_T_Cells" %in% names(result$countsWideRIRA))

  expect_equal(nrow(result$failures), 0)

  rhesus_only <- .quantifyWithMockedPrepare(
    Rdiscvr:::.ApplyHumanToRhesusRoleMap(cell_table, LoadQuantifyRoleMap()),
    rhesus_path,
    species = "rhesus",
    coerceToRIRA = FALSE
  )
  expect_equal(
    sort(names(result$countsWideRIRA)),
    sort(names(rhesus_only$countsWide))
  )
})

test_that("outputPrefix writes native and coerced counts TSV files", {
  cell_table <- .makeHumanCellTypistCellTable()
  prefix <- tempfile()
  on.exit(unlink(c(
    paste0(prefix, "_counts_wide.tsv"),
    paste0(prefix, "_counts_wide_coercedToRIRA.tsv")
  )))
  .quantifyWithMockedPrepare(
    cell_table,
    .humanImmuneQuantifySpecPath(),
    species = "human",
    coerceToRIRA = TRUE,
    outputPrefix = prefix
  )
  expect_true(file.exists(paste0(prefix, "_counts_wide.tsv")))
  expect_true(file.exists(paste0(prefix, "_counts_wide_coercedToRIRA.tsv")))
})

test_that("bundled quantify specs parse properly; rhesus RIRA panel quantifies", {
  rhesus_path <- .rhesusQuantifySpecPath()
  human_path <- .humanImmuneQuantifySpecPath()
  expect_gt(nrow(Rdiscvr:::.ParseQuantifySpec(rhesus_path)$spec), 150)
  expect_gt(nrow(Rdiscvr:::.ParseQuantifySpec(human_path)$spec), 700)

  result <- .quantifyWithMockedPrepare(.makeRhesusRiraCellTable(), rhesus_path)
  expect_equal(nrow(result$failures), 0)
  row_1001 <- result$countsWide[result$countsWide$cDNA_ID == 1001, ]
  expect_equal(row_1001$Immune__T_NK, 20L)
  expect_equal(row_1001$TNK__CD8plus_T_Cells, 8L)
  expect_equal(row_1001$Myeloid__coarse__Monocytes, 8L)
  expect_false(is.na(row_1001$TNK__CD8plus_T_Cells__CytotoxicityScore__median))
  expect_true("TNK__CD8plus_T_Cells__Diversity" %in% names(result$countsWide))
})
