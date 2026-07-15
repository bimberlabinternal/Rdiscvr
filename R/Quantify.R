#' @include Utils.R
#' @include SingleCell.R
#' @include Studies.R

#' @title Quantify10xData
#' @description Download full Seurat objects for LabKey output file IDs and quantify
#' grouped metrics from a TSV or data.frame spec.
#'
#' @param outputFileIds Integer vector of LabKey sequenceanalysis output file row IDs.
#' @param spec Path to a tab-separated spec file, or a data.frame with the required
#'   columns (see Details). A rhesus RIRA spec is bundled at
#'   \code{system.file("extdata", "quantify_rhesus_spec.tsv", package = "Rdiscvr")}.
#'   A human Immune_All_High/Low CellTypist spec is at
#'   \code{system.file("extdata", "quantify_human_immune_spec.tsv", package = "Rdiscvr")}.
#' @param failureLogFile Path for a text log of validation and processing failures.
#'   Defaults to \code{quantification_failures_<timestamp>.txt}. Also returned in
#'   the \code{failures} tibble.
#' @param classifyTNK If \code{TRUE} (default), per lane run TCR-aware TNK
#'   classification when RIRA immune and TNK metadata columns are present, then
#'   apply gamma-delta and NK overrides on \code{RIRA_TNK_v2.cellclass} within the
#'   \code{T_NK} compartment. Soft-skips when RIRA fields are absent. Set
#'   \code{FALSE} to quantify from downloaded RIRA labels only.
#' @param species \code{"rhesus"} or \code{"human"}. When \code{NULL} (default),
#'   inferred from lane metadata (\code{RIRA_Immune_v2.cellclass} vs
#'   \code{celltypist.Immune_All_High.cellclass}).
#' @param coerceToRIRA For human datasets, when \code{TRUE} (default), map CellTypist
#'   High/Low labels to RIRA fields via
#'   \code{quantify_human_to_rhesus_role_map.tsv} and also compute
#'   \code{countsWideRIRA} using the bundled rhesus spec. Ignored for rhesus.
#' @param outputPrefix Optional path prefix for writing wide count tables. When
#'   set, writes \code{<prefix>_counts_wide.tsv} and, for human with
#'   \code{coerceToRIRA = TRUE}, \code{<prefix>_counts_wide_coercedToRIRA.tsv}.
#' @return A named list with:
#' \describe{
#'   \item{\code{countsWide}}{Native species spec output (RIRA for rhesus;
#'     CellTypist High/Low for human).}
#'   \item{\code{countsWideRIRA}}{Human-only RIRA-shaped table when
#'     \code{coerceToRIRA = TRUE}; otherwise omitted.}
#'   \item{\code{failures}}{Tibble of deferred validation and processing failures
#'     (\code{specRow}, \code{outputFileId}, \code{field}, \code{reason}).}
#'   \item{\code{laneSummary}}{Per-\code{sourceOutputFileId} ingest counts
#'     (\code{nCellsIngested}).}
#' }
#' @details
#' Each spec row defines a grouping, cell filter, and a quantification type:
#' \code{sum} (cell count), \code{score} (quantiles of a numeric column),
#' \code{pct_positive} (percent of matched cells with a gene detected in counts),
#' or \code{diversity} (TCR clone Hill number q=2 via rarefaction with
#' \pkg{iNEXT}).
#'
#' Required columns: \code{GroupingVariable} (pipe-separated metadata column
#' names), \code{TargetField}, \code{SourceField}, \code{MatchValue},
#' \code{QuantificationType}, \code{QuantificationSourceField} (required for
#' \code{score} and \code{pct_positive}), and \code{QuantificationScoreType}
#' (\code{quantiles} when type is \code{score}).
#'
#' Optional columns \code{ScopeField} and \code{ScopeMatchValue} restrict
#' matching before the \code{SourceField} filter. Optional effector-
#' differentiation subset columns further restrict matched cells by a numeric
#' score bin.
#'
#' Per lane, missing standard UCell score columns (names ending in \code{_UCell})
#' are computed with \code{RIRA::CalculateUCellScores}; \code{Proliferation_UCell}
#' is added separately when needed. \code{Is_TCR_Stimulated} is predicted with
#' \code{RIRA::PredictTcellActivation} when absent.
#'
#' All spec rows must share the same \code{GroupingVariable}. Processing
#' continues after failures; see \code{failures} and \code{failureLogFile}.
#'
#' Bundled-spec \code{TargetField} names use \code{__} between prefix segments.
#' Regenerate the rhesus RIRA spec with \code{inst/scripts/generate_quantify_rhesus_spec.py},
#' the human immune spec with \code{inst/scripts/generate_quantify_human_spec.py},
#' and the human→RIRA role map with
#' \code{inst/scripts/generate_quantify_human_to_rhesus_role_map.py}.
#' @export
#' @importFrom dplyr %>% group_by summarize mutate filter distinct left_join bind_rows count pull across
#' @importFrom tibble as_tibble tibble
#' @importFrom stats quantile
#' @importFrom Seurat DefaultAssay GetAssayData
Quantify10xData <- function(
  outputFileIds,
  spec,
  failureLogFile = paste0(
    "quantification_failures_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".txt"
  ),
  classifyTNK = TRUE,
  species = NULL,
  coerceToRIRA = TRUE,
  outputPrefix = NULL
) {
  if (missing(outputFileIds) || !length(outputFileIds)) {
    stop("outputFileIds must be a non-empty integer vector")
  }

  if (missing(spec)) {
    stop("spec must be a file path or data.frame")
  }

  output_file_ids <- as.integer(unname(outputFileIds))
  if (any(is.na(output_file_ids))) {
    stop("outputFileIds must be coercible to integer")
  }

  failures <- .EmptyFailuresTibble()
  parsed_spec <- .ParseQuantifySpec(spec)
  failures <- dplyr::bind_rows(failures, parsed_spec$failures)
  spec_table <- parsed_spec$spec

  resolved_species <- if (!is.null(species) && length(species) == 1 && nzchar(species)) {
    match.arg(tolower(species), c("rhesus", "human"))
  } else {
    NULL
  }

  apply_role_map <- isTRUE(coerceToRIRA) &&
    (identical(resolved_species, "human") || is.null(resolved_species))
  role_map <- if (apply_role_map) {
    LoadQuantifyRoleMap()
  } else {
    NULL
  }
  rhesus_spec_table <- if (apply_role_map) {
    .ParseQuantifySpec(.BundledRhesusQuantifySpecPath())$spec
  } else {
    NULL
  }

  lane_tables <- list()
  lane_summary_pieces <- list()

  for (output_file_id in output_file_ids) {
    lane_result <- tryCatch(
      {
        prepare_result <- .PrepareLaneCellTable(
          output_file_id = output_file_id,
          spec_table = spec_table,
          classify_tnk = classifyTNK,
          apply_role_map = apply_role_map,
          role_map = role_map,
          rhesus_spec_table = rhesus_spec_table
        )
        list(
          lane_table = prepare_result$cell_table,
          lane_summary = tibble::tibble(
            sourceOutputFileId = output_file_id,
            nCellsIngested = nrow(prepare_result$cell_table)
          ),
          failures = prepare_result$failures,
          detected_species = prepare_result$detected_species
        )
      },
      error = function(error_condition) {
        list(
          lane_table = NULL,
          lane_summary = tibble::tibble(
            sourceOutputFileId = output_file_id,
            nCellsIngested = 0L
          ),
          failures = tibble::tibble(
            specRow = NA_integer_,
            outputFileId = output_file_id,
            field = NA_character_,
            reason = paste0("download failed: ", conditionMessage(error_condition))
          ),
          detected_species = NULL
        )
      }
    )

    failures <- dplyr::bind_rows(failures, lane_result$failures)
    lane_summary_pieces[[length(lane_summary_pieces) + 1]] <- lane_result$lane_summary
    if (!is.null(lane_result$detected_species) && is.null(resolved_species)) {
      resolved_species <- lane_result$detected_species
    }
    if (!is.null(lane_result$lane_table)) {
      lane_tables[[length(lane_tables) + 1]] <- lane_result$lane_table
    }
  }

  lane_summary <- dplyr::bind_rows(lane_summary_pieces)

  if (length(lane_tables)) {
    cell_table <- tibble::as_tibble(do.call(rbind, lane_tables))
    rownames(cell_table) <- NULL
  } else {
    cell_table <- tibble::tibble()
  }

  if (is.null(resolved_species) && nrow(cell_table) && isTRUE(coerceToRIRA)) {
    resolved_species <- tryCatch(
      .DetectQuantifySpeciesFromCellTable(cell_table),
      error = function(error_condition) NULL
    )
  }

  native_loop <- .RunQuantifyMetricLoop(
    spec_table = spec_table,
    cell_table = cell_table,
    failures = failures
  )
  counts_wide <- native_loop$countsWide
  failures <- native_loop$failures

  counts_wide_rira <- NULL
  if (identical(resolved_species, "human") && isTRUE(coerceToRIRA)) {
    coerced_loop <- .RunQuantifyMetricLoop(
      spec_table = rhesus_spec_table,
      cell_table = cell_table,
      failures = failures
    )
    counts_wide_rira <- coerced_loop$countsWide
    failures <- coerced_loop$failures
  }

  .WriteFailureLog(failures, failure_log_file = failureLogFile)
  if (nrow(failures) > 0) {
    warning(
      nrow(failures),
      " quantification failure(s) recorded; see failures and ",
      failureLogFile
    )
  }

  if (!is.null(outputPrefix) && nzchar(outputPrefix)) {
    .WriteQuantifyCountsTsv(
      counts_wide,
      paste0(outputPrefix, "_counts_wide.tsv")
    )
    if (!is.null(counts_wide_rira)) {
      .WriteQuantifyCountsTsv(
        counts_wide_rira,
        paste0(outputPrefix, "_counts_wide_coercedToRIRA.tsv")
      )
    }
  }

  result <- list(
    countsWide = counts_wide,
    failures = failures,
    laneSummary = lane_summary
  )
  if (!is.null(counts_wide_rira)) {
    result$countsWideRIRA <- counts_wide_rira
  }
  result
}

## Helper Functions ##

#' @title LoadQuantifyRoleMap
#' @description Load the bundled human CellTypist → RIRA role map used by
#'   \code{Quantify10xData(coerceToRIRA = TRUE)}.
#' @param mapFile Optional path to a role-map TSV. Defaults to the bundled
#'   \code{quantify_human_to_rhesus_role_map.tsv}.
#' @return A data.frame with columns \code{humanSourceField}, \code{humanLabel},
#'   \code{rhesusTargetField}, \code{rhesusLabel}, and \code{lineageRole}.
#' @export
LoadQuantifyRoleMap <- function(mapFile = NULL) {
  if (is.null(mapFile)) {
    mapFile <- system.file(
      "extdata",
      "quantify_human_to_rhesus_role_map.tsv",
      package = "Rdiscvr"
    )
  }
  if (!file.exists(mapFile)) {
    stop("role map file does not exist: ", mapFile)
  }
  utils::read.delim(
    mapFile,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.BundledRhesusQuantifySpecPath <- function() {
  system.file("extdata", "quantify_rhesus_spec.tsv", package = "Rdiscvr")
}

.DetectQuantifySpeciesFromCellTable <- function(cell_table) {
  if ("RIRA_Immune_v2.cellclass" %in% names(cell_table)) {
    return("rhesus")
  }
  if ("celltypist.Immune_All_High.cellclass" %in% names(cell_table)) {
    return("human")
  }
  stop("Cannot auto-detect quantify species from cell table metadata")
}

.DetectQuantifySpeciesFromSeurat <- function(seurat_obj) {
  meta_cols <- names(seurat_obj@meta.data)
  if ("RIRA_Immune_v2.cellclass" %in% meta_cols) {
    return("rhesus")
  }
  if ("celltypist.Immune_All_High.cellclass" %in% meta_cols) {
    return("human")
  }
  NULL
}

.ApplyHumanToRhesusRoleMap <- function(cell_table, role_map) {
  high_field <- "celltypist.Immune_All_High.cellclass"
  low_field <- "celltypist.Immune_All_Low.cellclass"

  if (!high_field %in% names(cell_table) && !low_field %in% names(cell_table)) {
    return(cell_table)
  }

  low_values <- if (low_field %in% names(cell_table)) {
    as.character(cell_table[[low_field]])
  } else {
    rep(NA_character_, nrow(cell_table))
  }
  high_values <- if (high_field %in% names(cell_table)) {
    as.character(cell_table[[high_field]])
  } else {
    rep(NA_character_, nrow(cell_table))
  }

  for (target_field in unique(role_map$rhesusTargetField)) {
    low_rows <- role_map[
      role_map$humanSourceField == "Immune_All_Low" &
        role_map$rhesusTargetField == target_field,
      ,
      drop = FALSE
    ]
    high_rows <- role_map[
      role_map$humanSourceField == "Immune_All_High" &
        role_map$rhesusTargetField == target_field,
      ,
      drop = FALSE
    ]
    low_lookup <- stats::setNames(low_rows$rhesusLabel, low_rows$humanLabel)
    high_lookup <- stats::setNames(high_rows$rhesusLabel, high_rows$humanLabel)

    mapped <- unname(low_lookup[low_values])
    use_high <- is.na(mapped) | !nzchar(mapped)
    if (any(use_high)) {
      mapped[use_high] <- unname(high_lookup[high_values[use_high]])
    }
    mapped[is.na(mapped) | !nzchar(mapped)] <- "Unknown"
    cell_table[[target_field]] <- mapped
  }

  cell_table
}

.ApplyHumanToRhesusRoleMapToSeurat <- function(seurat_obj, role_map) {
  cell_table <- tibble::as_tibble(seurat_obj@meta.data)
  mapped <- .ApplyHumanToRhesusRoleMap(cell_table, role_map)
  for (target_field in unique(role_map$rhesusTargetField)) {
    seurat_obj[[target_field]] <- mapped[[target_field]]
  }
  seurat_obj
}

.RunQuantifyMetricLoop <- function(spec_table, cell_table, failures) {
  counts_wide <- tibble::tibble()
  target_fields_seen <- character(0)

  if (nrow(spec_table) == 0) {
    return(list(countsWide = counts_wide, failures = failures))
  }

  grouping_specs <- unique(spec_table$GroupingVariable)
  mixed_grouping <- length(grouping_specs) > 1
  if (mixed_grouping) {
    failures <- dplyr::bind_rows(
      failures,
      tibble::tibble(
        specRow = NA_integer_,
        outputFileId = NA_integer_,
        field = "GroupingVariable",
        reason = paste0(
          "mixed grouping across spec rows: ",
          paste(grouping_specs, collapse = "; ")
        )
      )
    )
    return(list(countsWide = counts_wide, failures = failures))
  }

  grouping_columns <- trimws(strsplit(as.character(grouping_specs[[1]]), "|", fixed = TRUE)[[1]])
  grouping_columns <- grouping_columns[nzchar(grouping_columns)]

  if (nrow(cell_table) > 0 && length(grouping_columns)) {
    counts_wide <- cell_table %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(grouping_columns)))
  } else if (length(grouping_columns)) {
    empty_group_table <- as.data.frame(
      stats::setNames(
        rep(list(integer(0)), length(grouping_columns)),
        grouping_columns
      ),
      stringsAsFactors = FALSE
    )
    counts_wide <- tibble::as_tibble(empty_group_table)
  }

  rarefaction_level <- .ComputeSharedRarefactionLevel(
    spec_table = spec_table,
    cell_table = cell_table,
    grouping_columns = grouping_columns
  )

  for (spec_row_idx in seq_len(nrow(spec_table))) {
    spec_row <- spec_table[spec_row_idx, , drop = FALSE]
    original_spec_row <- as.integer(spec_row$specRow[[1]])
    row_failures <- .ValidateSpecRow(spec_row, cell_table)
    failures <- dplyr::bind_rows(failures, row_failures)
    if (nrow(row_failures) > 0) {
      next
    }

    target_field <- as.character(spec_row$TargetField[[1]])
    quantification_type <- tolower(trimws(as.character(spec_row$QuantificationType[[1]])))
    output_columns <- if (identical(quantification_type, "score")) {
      .BuildScoreQuantileColumnNames(target_field)
    } else {
      target_field
    }
    duplicate_targets <- intersect(output_columns, target_fields_seen)
    if (length(duplicate_targets)) {
      failures <- dplyr::bind_rows(
        failures,
        tibble::tibble(
          specRow = original_spec_row,
          outputFileId = NA_integer_,
          field = "TargetField",
          reason = paste0(
            "duplicate TargetField output column(s): ",
            paste(duplicate_targets, collapse = ", ")
          )
        )
      )
      next
    }

    metrics_table <- .ComputeMetricsForSpecRow(
      spec_row = spec_row,
      cell_table = cell_table,
      grouping_columns = grouping_columns,
      rarefaction_level = rarefaction_level
    )

    if (nrow(counts_wide) == 0 && nrow(metrics_table) > 0) {
      counts_wide <- metrics_table
    } else if (nrow(metrics_table) > 0) {
      counts_wide <- dplyr::left_join(
        counts_wide,
        metrics_table,
        by = grouping_columns
      )
    }

    target_fields_seen <- c(target_fields_seen, output_columns)
  }

  list(countsWide = counts_wide, failures = failures)
}

.WriteQuantifyCountsTsv <- function(counts_wide, output_path) {
  utils::write.table(
    counts_wide,
    file = output_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  invisible(output_path)
}

# register Quantify pct-positive gene panels into RIRA's gene-set registry
.RegisterQuantifyGeneSets <- function() {
  if (!requireNamespace("RIRA", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  register_if_absent <- function(name, genes) {
    tryCatch(
      RIRA:::.RegisterGeneSet(name, genes),
      error = function(error_condition) invisible(NULL)
    )
  }

  register_if_absent(
    "Quantify.PctPositive",
    c(
      "PDCD1",
      "KLRK1",
      "KLRB1",
      "HAVCR2",
      "TIGIT",
      "PRF1",
      "GZMB",
      "KLRC1",
      "FCGR3",
      "FCGR3A",
      "FOXP3",
      "IL2RA"
    )
  )
  register_if_absent(
    "Quantify.PctPositive.Inhibitory",
    c("PDCD1", "HAVCR2", "TIGIT")
  )
  register_if_absent(
    "Quantify.PctPositive.NK",
    c("KLRK1", "KLRC1", "FCGR3A", "FCGR3", "KLRB1")
  )
  register_if_absent(
    "Quantify.PctPositive.Treg",
    c("FOXP3", "IL2RA")
  )
  register_if_absent(
    "Quantify.PctPositive.Effector",
    c("GZMB", "KLRK1")
  )

  invisible(TRUE)
}

.RegisterQuantifyGeneSets()

.QuantifyNameDelim <- "__"

.BuildScoreQuantileColumnNames <- function(target_field) {
  paste0(target_field, .QuantifyNameDelim, c("p05", "median", "p95"))
}

#returns an empty tibble with columns for spec row, output file, field name, and reason
.EmptyFailuresTibble <- function() {
  tibble::tibble(
    specRow = integer(0),
    outputFileId = integer(0),
    field = character(0),
    reason = character(0)
  )
}

#assumes output_file_id refers to a valid LabKey output file row
.DownloadSeuratPerLane <- function(output_file_id) {
  out_file <- tempfile(fileext = ".rds")
  on.exit(unlink(out_file), add = TRUE)
  DownloadOutputFile(
    outputFileId = output_file_id,
    outFile = out_file,
    overwrite = TRUE
  )
  readRDS(out_file)
}

# download seurat, ensure scores/activation, optional human→RIRA map, classify TNK, extract cell table
.PrepareLaneCellTable <- function(
  output_file_id,
  spec_table,
  classify_tnk,
  apply_role_map = FALSE,
  role_map = NULL,
  rhesus_spec_table = NULL
) {
  seurat_obj <- .DownloadSeuratPerLane(output_file_id)
  detected_species <- .DetectQuantifySpeciesFromSeurat(seurat_obj)
  ensure_result <- .EnsureUCellAndActivationMetadata(seurat_obj, spec_table)
  seurat_obj <- ensure_result$seurat_obj
  failures <- ensure_result$failures

  lane_apply_map <- isTRUE(apply_role_map) && identical(detected_species, "human")
  if (lane_apply_map && !is.null(role_map) && nrow(role_map)) {
    seurat_obj <- .ApplyHumanToRhesusRoleMapToSeurat(seurat_obj, role_map)
  }

  lane_classify_tnk <- isTRUE(classify_tnk) && identical(detected_species, "rhesus")
  if (lane_classify_tnk) {
    classify_result <- .ClassifyTNKWithTcrOverrides(seurat_obj)
    seurat_obj <- classify_result$seurat_obj
    failures <- dplyr::bind_rows(failures, classify_result$failures)
  }

  combined_spec <- spec_table
  if (lane_apply_map && !is.null(rhesus_spec_table)) {
    combined_spec <- dplyr::bind_rows(combined_spec, rhesus_spec_table)
  }

  cell_table <- .BuildCellTableFromSeurat(seurat_obj, output_file_id, combined_spec)
  list(
    cell_table = cell_table,
    failures = failures,
    detected_species = detected_species
  )
}

.UCellColumnNeedsComputation <- function(seurat_obj, column_name) {
  if (!column_name %in% names(seurat_obj@meta.data)) {
    return(TRUE)
  }
  any(is.na(seurat_obj@meta.data[[column_name]]))
}

# ensure UCell scores and TCR activation metadata exist on the Seurat object
.EnsureUCellAndActivationMetadata <- function(seurat_obj, spec_table) {
  failures <- .EmptyFailuresTibble()
  append_failure <- function(field_name, reason_text) {
    failures <<- dplyr::bind_rows(
      failures,
      tibble::tibble(
        specRow = NA_integer_,
        outputFileId = NA_integer_,
        field = field_name,
        reason = reason_text
      )
    )
  }

  score_rows <- spec_table[
    tolower(trimws(spec_table$QuantificationType)) == "score",
    ,
    drop = FALSE
  ]
  ucell_fields <- unique(trimws(as.character(score_rows$QuantificationSourceField)))
  ucell_fields <- ucell_fields[grepl("_UCell$", ucell_fields, perl = TRUE)]
  standard_ucells <- setdiff(ucell_fields, "Proliferation_UCell")

  #RIRA::CalculateUCellScores covers standard *_UCell columns; Proliferation_UCell is handled below
  needs_standard_ucells <- standard_ucells[
    vapply(standard_ucells, .UCellColumnNeedsComputation, logical(1), seurat_obj = seurat_obj)
  ]

  if (length(needs_standard_ucells)) {
    if (!requireNamespace("RIRA", quietly = TRUE)) {
      append_failure(
        "RIRA",
        "RIRA package unavailable; skipping standard UCell score computation"
      )
    } else {
      seurat_obj <- tryCatch(
        RIRA::CalculateUCellScores(seurat_obj, plotCor = FALSE),
        error = function(error_condition) {
          append_failure(
            "RIRA::CalculateUCellScores",
            conditionMessage(error_condition)
          )
          seurat_obj
        }
      )
    }
  }

  #Proliferation_UCell is not in CalculateUCellScores; compute separately when needed
  if ("Proliferation_UCell" %in% ucell_fields &&
      .UCellColumnNeedsComputation(seurat_obj, "Proliferation_UCell")) {
    if (!requireNamespace("UCell", quietly = TRUE)) {
      append_failure(
        "UCell",
        "UCell package unavailable; skipping Proliferation_UCell computation"
      )
    } else {
      proliferation_genes <- tryCatch(
        RIRA::GetGeneSet("Proliferation"),
        error = function(error_condition) {
          append_failure(
            "RIRA::GetGeneSet",
            paste0(
              "GetGeneSet(Proliferation) failed: ",
              conditionMessage(error_condition),
              "; using MKI67/TOP2A fallback"
            )
          )
          c("MKI67", "TOP2A")
        }
      )
      seurat_obj <- tryCatch(
        UCell::AddModuleScore_UCell(
          seurat_obj,
          features = list(Proliferation = proliferation_genes),
          missing_genes = "skip"
        ),
        error = function(error_condition) {
          append_failure(
            "UCell::AddModuleScore_UCell",
            conditionMessage(error_condition)
          )
          seurat_obj
        }
      )
    }
  }

  #predict TCR stimulation status when the metadata column is absent
  if (!"Is_TCR_Stimulated" %in% names(seurat_obj@meta.data)) {
    if (!requireNamespace("RIRA", quietly = TRUE)) {
      append_failure(
        "RIRA",
        "RIRA package unavailable; skipping Is_TCR_Stimulated prediction"
      )
    } else {
      seurat_obj <- tryCatch(
        RIRA::PredictTcellActivation(seurat_obj),
        error = function(error_condition) {
          append_failure(
            "RIRA::PredictTcellActivation",
            conditionMessage(error_condition)
          )
          seurat_obj
        }
      )
    }
  }

  list(seurat_obj = seurat_obj, failures = failures)
}

# run ClassifyTNKByExpression when possible, then apply TCR TNK overrides
.ClassifyTNKWithTcrOverrides <- function(seurat_obj) {
  failures <- .EmptyFailuresTibble()
  immune_field <- "RIRA_Immune_v2.cellclass"
  tnk_field <- "RIRA_TNK_v2.cellclass"

  if (!immune_field %in% names(seurat_obj@meta.data) ||
      !tnk_field %in% names(seurat_obj@meta.data)) {
    return(list(seurat_obj = seurat_obj, failures = failures))
  }

  meta_cols <- names(seurat_obj@meta.data)
  has_cdr3 <- "HasCDR3Data" %in% meta_cols &&
    any(as.logical(seurat_obj$HasCDR3Data), na.rm = TRUE)
  has_tra <- "TRA" %in% meta_cols && any(!is.na(seurat_obj$TRA))
  has_trb <- "TRB" %in% meta_cols && any(!is.na(seurat_obj$TRB))

  if (has_cdr3 || has_tra || has_trb) {
    seurat_obj <- tryCatch(
      ClassifyTNKByExpression(seurat_obj),
      error = function(error_condition) seurat_obj
    )
  }

  cell_table <- tibble::as_tibble(seurat_obj@meta.data)
  override_result <- .ApplyTcrTNKOverridesToCellTable(cell_table)
  seurat_obj[[tnk_field]] <- override_result$cell_table[[tnk_field]]

  list(
    seurat_obj = seurat_obj,
    failures = dplyr::bind_rows(failures, override_result$failures)
  )
}

# extract metadata and materialize pct_positive gene columns from counts
.BuildCellTableFromSeurat <- function(seurat_obj, output_file_id, spec_table) {
  cell_table <- tibble::as_tibble(seurat_obj@meta.data)
  cell_table$sourceOutputFileId <- output_file_id

  pct_genes <- spec_table %>%
    dplyr::filter(tolower(trimws(.data$QuantificationType)) == "pct_positive") %>%
    dplyr::pull(.data$QuantificationSourceField) %>%
    unique()
  pct_genes <- trimws(as.character(pct_genes))
  pct_genes <- pct_genes[nzchar(pct_genes)]

  if (length(pct_genes)) {
    assay_name <- Seurat::DefaultAssay(seurat_obj)
    counts <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
    gene_rownames <- rownames(counts)
    for (gene in pct_genes) {
      if (gene %in% gene_rownames) {
        cell_table[[gene]] <- as.numeric(counts[gene, ] > 0)
      }
    }
  }

  cell_table
}

#assumes spec is a file path or data.frame containing the required columns
#returns the cleaned spec table and any failures from blank or malformed rows
.ParseQuantifySpec <- function(spec) {
  required_columns <- c(
    "GroupingVariable",
    "TargetField",
    "SourceField",
    "MatchValue",
    "QuantificationType",
    "QuantificationSourceField",
    "QuantificationScoreType"
  )
  optional_columns <- c(
    "ScopeField",
    "ScopeMatchValue",
    "EffectorDifferentiationScoreField",
    "EffectorDifferentiationCutpointLow",
    "EffectorDifferentiationCutpointHigh",
    "SubsetPhenotypeOutputFieldName"
  )

  failures <- .EmptyFailuresTibble()

  if (is.character(spec) && length(spec) == 1) {
    if (!file.exists(spec)) {
      stop("spec file does not exist: ", spec)
    }
    spec_table <- utils::read.delim(
      spec,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (is.data.frame(spec)) {
    spec_table <- as.data.frame(spec, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("spec must be a file path or data.frame")
  }

  spec_table <- tibble::as_tibble(spec_table)
  missing_columns <- setdiff(required_columns, names(spec_table))
  if (length(missing_columns)) {
    stop(
      "spec is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  for (optional_column in optional_columns) {
    if (!optional_column %in% names(spec_table)) {
      spec_table[[optional_column]] <- ""
    } else {
      spec_table[[optional_column]] <- ifelse(
        is.na(spec_table[[optional_column]]),
        "",
        as.character(spec_table[[optional_column]])
      )
    }
  }

  spec_table <- spec_table[, c(required_columns, optional_columns), drop = FALSE]
  spec_table$specRow <- seq_len(nrow(spec_table))

  blank_row_idx <- which(
    !nzchar(trimws(spec_table$GroupingVariable)) &
      !nzchar(trimws(spec_table$TargetField)) &
      !nzchar(trimws(spec_table$SourceField)) &
      !nzchar(trimws(spec_table$MatchValue)) &
      !nzchar(trimws(spec_table$QuantificationType))
  )
  if (length(blank_row_idx)) {
    failures <- dplyr::bind_rows(
      failures,
      tibble::tibble(
        specRow = blank_row_idx,
        outputFileId = NA_integer_,
        field = NA_character_,
        reason = "malformed spec row: all required fields are blank"
      )
    )
    spec_table <- spec_table[-blank_row_idx, , drop = FALSE]
  }

  list(
    spec = spec_table,
    failures = failures
  )
}

#assumes cell_table is seurat object metadata with RIRA immune and T/NK class columns
#overrides only gamma-delta and NK in RIRA_TNK_v2.cellclass within T_NK scope
.ApplyTcrTNKOverridesToCellTable <- function(cell_table) {
  failures <- .EmptyFailuresTibble()
  immune_field <- "RIRA_Immune_v2.cellclass"
  tnk_field <- "RIRA_TNK_v2.cellclass"

  if (!immune_field %in% names(cell_table) || !tnk_field %in% names(cell_table)) {
    return(list(cell_table = cell_table, failures = failures))
  }

  has_cdr3_data <- if ("HasCDR3Data" %in% names(cell_table)) {
    as.logical(cell_table$HasCDR3Data)
  } else if ("CDR3s" %in% names(cell_table)) {
    !is.na(cell_table$CDR3s)
  } else {
    rep(FALSE, nrow(cell_table))
  }

  if (!any(has_cdr3_data, na.rm = TRUE)) {
    failures <- dplyr::bind_rows(
      failures,
      tibble::tibble(
        specRow = NA_integer_,
        outputFileId = NA_integer_,
        field = "HasCDR3Data",
        reason = "no TCR data available for TNK classification; GD/NK override skipped"
      )
    )
    return(list(cell_table = cell_table, failures = failures))
  }

  #column fallbacks mirror ClassifyTNKByExpression metadata fields
  is_gamma_delta <- if ("TNK_Type" %in% names(cell_table)) {
    as.character(cell_table$TNK_Type) == "Gamma/Delta"
  } else if ("IsGammaDelta" %in% names(cell_table)) {
    as.logical(cell_table$IsGammaDelta)
  } else if ("TRD" %in% names(cell_table)) {
    !is.na(cell_table$TRD)
  } else {
    rep(FALSE, nrow(cell_table))
  }

  is_nk_cell <- if ("TNK_Type" %in% names(cell_table)) {
    as.character(cell_table$TNK_Type) == "NK (CD3-/TCR-)"
  } else if ("IsNKCell" %in% names(cell_table)) {
    as.logical(cell_table$IsNKCell)
  } else {
    has_cd3 <- if ("HasCD3" %in% names(cell_table)) {
      as.logical(cell_table$HasCD3)
    } else {
      cd3_genes <- intersect(c("CD3D", "CD3E", "CD3G"), names(cell_table))
      if (!length(cd3_genes)) {
        NULL
      } else {
        has_cd3_values <- rep(FALSE, nrow(cell_table))
        for (gene in cd3_genes) {
          counts <- suppressWarnings(as.numeric(cell_table[[gene]]))
          has_cd3_values <- has_cd3_values | (!is.na(counts) & counts > 0)
        }
        has_cd3_values
      }
    }

    if (is.null(has_cd3)) {
      NULL
    } else {
      !has_cdr3_data & !has_cd3
    }
  }

  if (is.null(is_nk_cell)) {
    failures <- dplyr::bind_rows(
      failures,
      tibble::tibble(
        specRow = NA_integer_,
        outputFileId = NA_integer_,
        field = "HasCD3",
        reason = "HasCD3 unavailable for TNK classification; NK override skipped"
      )
    )
    is_nk_cell <- rep(FALSE, nrow(cell_table))
  }

  tnk_scope <- as.character(cell_table[[immune_field]]) == "T_NK"
  updated_tnk <- as.character(cell_table[[tnk_field]])

  updated_tnk[tnk_scope & is_gamma_delta] <- "Gamma/Delta Cells"
  updated_tnk[tnk_scope & is_nk_cell] <- "NK Cells"
  cell_table[[tnk_field]] <- updated_tnk

  list(cell_table = cell_table, failures = failures)
}

# scope, source/match, and EDS filters shared by quantify row handlers
.FilterCellsBySpecRow <- function(cell_table, spec_row) {
  if (nrow(cell_table) == 0) {
    return(cell_table)
  }

  scoped_cells <- cell_table
  scope_field <- trimws(as.character(spec_row$ScopeField[[1]]))
  scope_match_value <- as.character(spec_row$ScopeMatchValue[[1]])
  if (nzchar(scope_field)) {
    scope_values <- as.character(scoped_cells[[scope_field]])
    scoped_cells <- scoped_cells[scope_values == scope_match_value, , drop = FALSE]
  }

  source_field <- as.character(spec_row$SourceField[[1]])
  match_value <- as.character(spec_row$MatchValue[[1]])
  source_values <- as.character(scoped_cells[[source_field]])
  filtered_cells <- scoped_cells[source_values == match_value, , drop = FALSE]
  .FilterCellsByEffectorDifferentiationBin(filtered_cells, spec_row)
}

#assumes cell_table column names match ingested metadata when cells were materialized
#returns validation failures for this rule without stopping the rest of the run
.ValidateSpecRow <- function(spec_row, cell_table) {
  failures <- .EmptyFailuresTibble()
  append_failure <- function(field_name, reason_text) {
    failures <<- dplyr::bind_rows(
      failures,
      tibble::tibble(
        specRow = as.integer(spec_row$specRow[[1]]),
        outputFileId = NA_integer_,
        field = field_name,
        reason = reason_text
      )
    )
  }

  target_field <- trimws(as.character(spec_row$TargetField[[1]]))
  source_field <- trimws(as.character(spec_row$SourceField[[1]]))
  match_value <- as.character(spec_row$MatchValue[[1]])
  quantification_type <- tolower(trimws(as.character(spec_row$QuantificationType[[1]])))
  score_field <- trimws(as.character(spec_row$QuantificationSourceField[[1]]))
  score_type <- tolower(trimws(as.character(spec_row$QuantificationScoreType[[1]])))
  grouping_variable <- trimws(as.character(spec_row$GroupingVariable[[1]]))
  scope_field <- trimws(as.character(spec_row$ScopeField[[1]]))
  scope_match_value <- as.character(spec_row$ScopeMatchValue[[1]])
  eds_score_field <- trimws(as.character(spec_row$EffectorDifferentiationScoreField[[1]]))
  eds_cutpoint_low <- trimws(as.character(spec_row$EffectorDifferentiationCutpointLow[[1]]))
  eds_cutpoint_high <- trimws(as.character(spec_row$EffectorDifferentiationCutpointHigh[[1]]))
  phenotype_output_field <- trimws(as.character(spec_row$SubsetPhenotypeOutputFieldName[[1]]))

  if (!nzchar(target_field)) {
    append_failure("TargetField", "malformed spec row: TargetField is blank")
  }
  if (!nzchar(source_field)) {
    append_failure("SourceField", "malformed spec row: SourceField is blank")
  }
  if (!nzchar(match_value)) {
    append_failure("MatchValue", "malformed spec row: MatchValue is blank")
  }
  if (!nzchar(quantification_type)) {
    append_failure("QuantificationType", "malformed spec row: QuantificationType is blank")
  } else if (!quantification_type %in% c("sum", "score", "pct_positive", "diversity")) {
    append_failure(
      "QuantificationType",
      paste0("invalid QuantificationType: ", quantification_type)
    )
  }

  grouping_columns <- trimws(strsplit(as.character(grouping_variable), "|", fixed = TRUE)[[1]])
  grouping_columns <- grouping_columns[nzchar(grouping_columns)]
  if (!length(grouping_columns)) {
    append_failure(
      "GroupingVariable",
      "malformed spec row: GroupingVariable is blank"
    )
  }

  metadata_columns <- if (nrow(cell_table)) names(cell_table) else character(0)
  missing_grouping <- setdiff(grouping_columns, metadata_columns)
  if (length(missing_grouping)) {
    append_failure(
      "GroupingVariable",
      paste0(
        "grouping column(s) not found in metadata: ",
        paste(missing_grouping, collapse = ", ")
      )
    )
  }

  if (nzchar(source_field) && !source_field %in% metadata_columns) {
    append_failure(
      "SourceField",
      paste0("SourceField not found in metadata: ", source_field)
    )
  }

  if (xor(nzchar(scope_field), nzchar(scope_match_value))) {
    append_failure(
      "ScopeField",
      "malformed spec row: ScopeField and ScopeMatchValue must both be set or both be blank"
    )
  } else if (nzchar(scope_field) && !scope_field %in% metadata_columns) {
    append_failure(
      "ScopeField",
      paste0("ScopeField not found in metadata: ", scope_field)
    )
  }

  if (identical(quantification_type, "score")) {
    if (!nzchar(score_field)) {
      append_failure(
        "QuantificationSourceField",
        "malformed spec row: QuantificationSourceField is required when QuantificationType is score"
      )
    } else if (!score_field %in% metadata_columns) {
      append_failure(
        "QuantificationSourceField",
        paste0(
          "QuantificationSourceField not found in metadata: ",
          score_field
        )
      )
    }

    if (!nzchar(score_type)) {
      append_failure(
        "QuantificationScoreType",
        "malformed spec row: QuantificationScoreType is required when QuantificationType is score"
      )
    } else if (!identical(score_type, "quantiles")) {
      append_failure(
        "QuantificationScoreType",
        paste0("invalid QuantificationScoreType: ", score_type)
      )
    }
  }

  if (identical(quantification_type, "pct_positive")) {
    if (!nzchar(score_field)) {
      append_failure(
        "QuantificationSourceField",
        "malformed spec row: QuantificationSourceField is required when QuantificationType is pct_positive"
      )
    } else if (!score_field %in% metadata_columns) {
      append_failure(
        "QuantificationSourceField",
        paste0(
          "QuantificationSourceField not found in cell table: ",
          score_field
        )
      )
    }
  }

  if (identical(quantification_type, "diversity")) {
    if (!"TRA" %in% metadata_columns) {
      append_failure("TRA", "TRA not found in metadata")
    }
    if (!"TRB" %in% metadata_columns) {
      append_failure("TRB", "TRB not found in metadata")
    }
  }

  eds_subset_columns_set <- c(
    nzchar(eds_score_field),
    nzchar(eds_cutpoint_low),
    nzchar(eds_cutpoint_high),
    nzchar(phenotype_output_field)
  )
  if (any(eds_subset_columns_set) && !all(eds_subset_columns_set)) {
    if (!nzchar(phenotype_output_field)) {
      append_failure(
        "SubsetPhenotypeOutputFieldName",
        "malformed spec row: SubsetPhenotypeOutputFieldName is required when other effector-differentiation subset columns are set"
      )
    }
    if (!nzchar(eds_score_field)) {
      append_failure(
        "EffectorDifferentiationScoreField",
        "malformed spec row: EffectorDifferentiationScoreField is required when SubsetPhenotypeOutputFieldName is set"
      )
    }
    if (!nzchar(eds_cutpoint_low)) {
      append_failure(
        "EffectorDifferentiationCutpointLow",
        "malformed spec row: EffectorDifferentiationCutpointLow is required when SubsetPhenotypeOutputFieldName is set"
      )
    }
    if (!nzchar(eds_cutpoint_high)) {
      append_failure(
        "EffectorDifferentiationCutpointHigh",
        "malformed spec row: EffectorDifferentiationCutpointHigh is required when SubsetPhenotypeOutputFieldName is set"
      )
    }
  }

  if (nzchar(phenotype_output_field)) {
    valid_phenotypes <- c("naive", "memorylike", "effector")
    if (!tolower(phenotype_output_field) %in% valid_phenotypes) {
      append_failure(
        "SubsetPhenotypeOutputFieldName",
        paste0("invalid SubsetPhenotypeOutputFieldName: ", phenotype_output_field)
      )
    }

    low_value <- suppressWarnings(as.numeric(eds_cutpoint_low))
    high_value <- suppressWarnings(as.numeric(eds_cutpoint_high))
    if (is.na(low_value)) {
      append_failure(
        "EffectorDifferentiationCutpointLow",
        paste0("invalid EffectorDifferentiationCutpointLow: ", eds_cutpoint_low)
      )
    }
    if (is.na(high_value)) {
      append_failure(
        "EffectorDifferentiationCutpointHigh",
        paste0("invalid EffectorDifferentiationCutpointHigh: ", eds_cutpoint_high)
      )
    }
    if (!is.na(low_value) && !is.na(high_value) && low_value > high_value) {
      append_failure(
        "EffectorDifferentiationCutpointLow",
        "EffectorDifferentiationCutpointLow must be less than or equal to EffectorDifferentiationCutpointHigh"
      )
    }

    if (nzchar(eds_score_field) && !eds_score_field %in% metadata_columns) {
      append_failure(
        "EffectorDifferentiationScoreField",
        paste0(
          "EffectorDifferentiationScoreField not found in metadata: ",
          eds_score_field
        )
      )
    }
  }

  failures
}

#restrict matched cells to Naive/MemoryLike/Effector bin of EffectorDifferentiationScoreField
.FilterCellsByEffectorDifferentiationBin <- function(cells, spec_row) {
  phenotype_output_field <- trimws(as.character(spec_row$SubsetPhenotypeOutputFieldName[[1]]))
  if (!nzchar(phenotype_output_field) || !nrow(cells)) {
    return(cells)
  }

  score_field <- trimws(as.character(spec_row$EffectorDifferentiationScoreField[[1]]))
  low <- as.numeric(spec_row$EffectorDifferentiationCutpointLow[[1]])
  high <- as.numeric(spec_row$EffectorDifferentiationCutpointHigh[[1]])
  scores <- suppressWarnings(as.numeric(cells[[score_field]]))
  keep <- switch(
    tolower(phenotype_output_field),
    naive = !is.na(scores) & scores < low,
    memorylike = !is.na(scores) & scores >= low & scores <= high,
    effector = !is.na(scores) & scores > high,
    rep(FALSE, nrow(cells))
  )
  cells[keep, , drop = FALSE]
}

# minimum TCR+ cell count per group across diversity spec rows for shared rarefaction
.ComputeSharedRarefactionLevel <- function(spec_table, cell_table, grouping_columns) {
  diversity_rows <- spec_table[
    tolower(trimws(spec_table$QuantificationType)) == "diversity",
    ,
    drop = FALSE
  ]
  if (!nrow(diversity_rows) || !nrow(cell_table) || !length(grouping_columns)) {
    return(NA_integer_)
  }
  if (!all(c("TRA", "TRB") %in% names(cell_table))) {
    return(NA_integer_)
  }

  min_counts <- integer(0)
  for (diversity_row_idx in seq_len(nrow(diversity_rows))) {
    spec_row <- diversity_rows[diversity_row_idx, , drop = FALSE]
    filtered_cells <- .FilterCellsBySpecRow(cell_table, spec_row)
    tra <- as.character(filtered_cells$TRA)
    trb <- as.character(filtered_cells$TRB)
    tcr_plus <- !is.na(tra) & nzchar(tra) & !is.na(trb) & nzchar(trb)
    filtered_tcr <- filtered_cells[tcr_plus, , drop = FALSE]

    if (nrow(filtered_tcr)) {
      group_counts <- filtered_tcr %>%
        dplyr::count(dplyr::across(dplyr::all_of(grouping_columns)), name = "n_tcr")
      min_counts <- c(min_counts, min(group_counts$n_tcr))
    }
  }

  if (!length(min_counts)) {
    NA_integer_
  } else {
    as.integer(min(min_counts))
  }
}

# computes Hill number q=2 at shared rarefaction level per grouping key
.ComputeRarefiedHillDiversity <- function(
  filtered_cells,
  all_groups,
  grouping_columns,
  target_field,
  rarefaction_level
) {
  all_groups[[target_field]] <- NA_real_

  if (is.na(rarefaction_level) || rarefaction_level < 1L || nrow(filtered_cells) == 0) {
    return(all_groups)
  }

  if (!requireNamespace("iNEXT", quietly = TRUE)) {
    return(all_groups)
  }

  tra <- as.character(filtered_cells$TRA)
  trb <- as.character(filtered_cells$TRB)
  tcr_plus <- !is.na(tra) & nzchar(tra) & !is.na(trb) & nzchar(trb)
  tcr_cells <- filtered_cells[tcr_plus, , drop = FALSE]

  if (!nrow(tcr_cells)) {
    return(all_groups)
  }

  #clone identity is TRA/TRB, plus TRA_V/TRB_V when present
  clone_cols <- c("TRA", "TRB")
  if ("TRA_V" %in% names(tcr_cells)) {
    clone_cols <- c(clone_cols, "TRA_V")
  }
  if ("TRB_V" %in% names(tcr_cells)) {
    clone_cols <- c(clone_cols, "TRB_V")
  }

  clone_sizes <- tcr_cells %>%
    dplyr::count(dplyr::across(dplyr::all_of(c(grouping_columns, clone_cols))), name = "clone_size")

  group_key <- interaction(clone_sizes[, grouping_columns, drop = FALSE], drop = TRUE)
  split_clones <- split(clone_sizes$clone_size, group_key)

  diversity_by_key <- vapply(
    names(split_clones),
    function(key_name) {
      abund <- as.integer(split_clones[[key_name]])
      if (sum(abund) < rarefaction_level) {
        return(NA_real_)
      }
      est <- tryCatch(
        iNEXT::estimateD(
          list(Group = abund),
          q = 2,
          datatype = "abundance",
          base = "size",
          level = as.integer(rarefaction_level)
        ),
        error = function(error_condition) NULL
      )
      if (is.null(est) || !is.data.frame(est) || !nrow(est)) {
        return(NA_real_)
      }
      qd_values <- est$qD[est$Order.q == 2]
      if (!length(qd_values)) {
        NA_real_
      } else {
        as.numeric(qd_values[[1]])
      }
    },
    FUN.VALUE = numeric(1)
  )

  all_group_keys <- interaction(all_groups[, grouping_columns, drop = FALSE], drop = TRUE)
  matched_idx <- match(as.character(all_group_keys), names(diversity_by_key))
  all_groups[[target_field]] <- as.numeric(diversity_by_key[matched_idx])
  all_groups
}

#assumes the rule passed validation and grouping columns are present in cell_table
#returns one row per group with the target metric column(s) for this rule
.ComputeMetricsForSpecRow <- function(
  spec_row,
  cell_table,
  grouping_columns,
  rarefaction_level = NA_integer_
) {
  target_field <- as.character(spec_row$TargetField[[1]])
  quantification_type <- tolower(trimws(as.character(spec_row$QuantificationType[[1]])))

  if (nrow(cell_table) == 0 || !length(grouping_columns)) {
    return(tibble::tibble())
  }

  all_groups <- cell_table %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(grouping_columns)))

  filtered_cells <- .FilterCellsBySpecRow(cell_table, spec_row)

  #sum: cell counts per group
  if (identical(quantification_type, "sum")) {
    if (nrow(filtered_cells) == 0) {
      count_values <- rep(0L, nrow(all_groups))
    } else {
      count_values <- filtered_cells %>%
        dplyr::count(dplyr::across(dplyr::all_of(grouping_columns)), name = target_field) %>%
        dplyr::right_join(all_groups, by = grouping_columns) %>%
        dplyr::mutate(!!target_field := ifelse(is.na(.data[[target_field]]), 0L, .data[[target_field]])) %>%
        dplyr::arrange(dplyr::across(dplyr::all_of(grouping_columns)))
      return(count_values)
    }

    all_groups[[target_field]] <- as.integer(count_values)
    return(all_groups)
  }

  #pct_positive: percent of matched cells with the gene detected
  if (identical(quantification_type, "pct_positive")) {
    gene_field <- as.character(spec_row$QuantificationSourceField[[1]])

    if (nrow(filtered_cells) == 0) {
      all_groups[[target_field]] <- 0
      return(all_groups)
    }

    gene_values <- filtered_cells[[gene_field]]
    is_positive <- if (is.logical(gene_values)) {
      gene_values
    } else {
      numeric_values <- suppressWarnings(as.numeric(gene_values))
      !is.na(numeric_values) & numeric_values > 0
    }
    filtered_cells$pct_positive_flag <- is_positive

    pct_metrics <- filtered_cells %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) %>%
      dplyr::summarize(
        !!target_field := 100 * mean(.data[["pct_positive_flag"]], na.rm = TRUE),
        .groups = "drop"
      )

    dplyr::left_join(all_groups, pct_metrics, by = grouping_columns) %>%
      dplyr::mutate(!!target_field := ifelse(is.na(.data[[target_field]]), 0, .data[[target_field]]))
  } else if (identical(quantification_type, "diversity")) {
    #diversity: rarefied Hill number q=2 at the shared TCR+ depth
    .ComputeRarefiedHillDiversity(
      filtered_cells = filtered_cells,
      all_groups = all_groups,
      grouping_columns = grouping_columns,
      target_field = target_field,
      rarefaction_level = rarefaction_level
    )
  } else {
    #score: p05/median/p95 of QuantificationSourceField among matched cells
    score_field <- as.character(spec_row$QuantificationSourceField[[1]])
    quantile_columns <- .BuildScoreQuantileColumnNames(target_field)
    p05_col <- quantile_columns[[1]]
    median_col <- quantile_columns[[2]]
    p95_col <- quantile_columns[[3]]

    if (nrow(filtered_cells) == 0) {
      all_groups[[p05_col]] <- NA_real_
      all_groups[[median_col]] <- NA_real_
      all_groups[[p95_col]] <- NA_real_
      return(all_groups)
    }

    score_values <- suppressWarnings(as.numeric(filtered_cells[[score_field]]))
    filtered_cells[[score_field]] <- score_values

    quantile_metrics <- filtered_cells %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) %>%
      dplyr::summarize(
        !!p05_col := {
          score_vector <- .data[[score_field]]
          score_vector <- score_vector[!is.na(score_vector)]
          if (!length(score_vector)) {
            NA_real_
          } else {
            as.numeric(stats::quantile(score_vector, probs = 0.05, na.rm = TRUE, names = FALSE))
          }
        },
        !!median_col := {
          score_vector <- .data[[score_field]]
          score_vector <- score_vector[!is.na(score_vector)]
          if (!length(score_vector)) {
            NA_real_
          } else {
            as.numeric(stats::quantile(score_vector, probs = 0.5, na.rm = TRUE, names = FALSE))
          }
        },
        !!p95_col := {
          score_vector <- .data[[score_field]]
          score_vector <- score_vector[!is.na(score_vector)]
          if (!length(score_vector)) {
            NA_real_
          } else {
            as.numeric(stats::quantile(score_vector, probs = 0.95, na.rm = TRUE, names = FALSE))
          }
        },
        .groups = "drop"
      )

    dplyr::left_join(all_groups, quantile_metrics, by = grouping_columns)
  }
}

#assumes failure_log_file is a path where we can write a text log
.WriteFailureLog <- function(failures, failure_log_file) {
  if (!nrow(failures)) {
    writeLines(character(0), con = failure_log_file)
    return(invisible(failure_log_file))
  }

  log_lines <- apply(failures, 1, function(failure_row) {
    paste0(
      "row=", failure_row[["specRow"]],
      " outputFileId=", failure_row[["outputFileId"]],
      " field=", failure_row[["field"]],
      " reason=", failure_row[["reason"]]
    )
  })

  writeLines(log_lines, con = failure_log_file)
  invisible(failure_log_file)
}
