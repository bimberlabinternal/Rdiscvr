#' @include Utils.R
#' @include SingleCell.R
#' @include Studies.R

#' @title Quantify10xData
#' @description Download Seurat Object metadata for output file IDs and quantify
#' metrics from a TSV spec. Each spec row defines a grouping, filter, and either a
#' raw count (\code{sum}) or quantiles of a score column (\code{score}).
#'
#' @param outputFileIds Integer vector of output file row IDs.
#' @param spec Path to a tab-separated spec file or a data.frame with the required
#'   columns (see Details).
#' @param failureLogFile Path for a text log of validation and processing failures.
#'   Defaults to \code{quantification_failures_<timestamp>.txt}.
#' @param classifyTNK If \code{TRUE} (default), gamma-delta and NK labels in
#'   \code{RIRA_TNK_v2.cellclass} are overridden within the \code{T_NK} immune
#'   compartment using TCR-based rules matching \code{ClassifyTNKByExpression}.
#'   CD4/CD8 and other T/NK subtypes remain RIRA CellTypist labels. Set
#'   \code{FALSE} to quantify from RIRA labels only.
#' @return A list with \code{countsWide} (grouping columns plus one column per
#'   \code{TargetField}, or \code{__p05}/\code{__median}/\code{__p95} for quantile
#'   rules), \code{failures} (tibble of deferred failures), and \code{laneSummary}
#'   (per-\code{sourceOutputFileId} ingest summary).
#' @details
#' Required spec columns: \code{GroupingVariable} (pipe-separated metadata
#' column names, e.g. \code{cDNA_ID} or \code{SubjectId|Tissue|TimepointLabel}),
#' \code{TargetField}, \code{SourceField}, \code{MatchValue}, \code{QuantificationType}
#' (\code{sum} or \code{score}), \code{QuantificationSourceField} (required when
#' type is \code{score}), and \code{QuantificationScoreType} (\code{quantiles}
#' only, when type is \code{score}). Optional columns \code{ScopeField} and
#' \code{ScopeMatchValue} restrict matching to cells where
#' \code{ScopeField == ScopeMatchValue} before applying the \code{SourceField}
#' filter. Optional columns \code{EffectorDifferentiationScoreField},
#' \code{EffectorDifferentiationCutpointLow},
#' \code{EffectorDifferentiationCutpointHigh}, and
#' \code{SubsetPhenotypeOutputFieldName} further restrict matched
#' cells by a numeric score bin (\code{Naive}: score below low cutpoint;
#' \code{MemoryLike}: score between cutpoints inclusive; \code{Effector}: score
#' above high cutpoint). The score column must already exist in metadata (e.g.
#' \code{Tcell_EffectorDifferentiation} from RIRA \code{ScoreUsingSavedComponent}).
#'
#' \code{TargetField} names in the default spec use \code{__} as a structural
#' delimiter between prefix segments and sanitized class labels; quantile rules
#' append \code{__p05}, \code{__median}, and \code{__p95} to \code{TargetField}.
#'
#' All spec rows must share the same \code{GroupingVariable}. Sum rules count
#' cells where \code{SourceField == MatchValue} per group (0 when no matches).
#' Score rules compute p05, median, and p95 of \code{QuantificationSourceField}
#' among matched cells (NA when no matches). Processing continues after failures;
#' see \code{failures} and \code{failureLogFile}.
#'
#' A default spec is bundled at
#' \code{system.file("extdata", "quantify_default_spec.tsv", package = "Rdiscvr")}.
#' Regenerate it from RIRA CellTypist pkl models with
#' \code{inst/scripts/generate_quantify_default_spec.py} (see script header).
#'
#' When \code{classifyTNK = TRUE}, gamma-delta and NK quantification rows use
#' TCR evidence (\code{TNK_Type}, \code{TRD}, \code{HasCDR3Data}, \code{HasCD3})
#' from the downloaded metadata rather than RIRA CellTypist for those two
#' subtypes. Original RIRA labels are retained in
#' \code{RIRA_TNK_v2.cellclass_rira}.
#' @export
#' @importFrom dplyr %>% group_by summarize mutate filter distinct left_join bind_rows count
#' @importFrom tibble as_tibble tibble
#' @importFrom stats quantile
Quantify10xData <- function(
  outputFileIds,
  spec,
  failureLogFile = paste0(
    "quantification_failures_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".txt"
  ),
  classifyTNK = TRUE
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
  #santize spec
  parsed_spec <- .ParseQuantifySpec(spec)
  failures <- dplyr::bind_rows(failures, parsed_spec$failures)
  #sanitized spec table
  spec_table <- parsed_spec$spec

  lane_tables <- list()
  lane_summary_pieces <- list()

  #download cell metadata for each output file, and record failures when a download fails
  for (output_file_id in output_file_ids) {
    lane_result <- tryCatch(
      {
        lane_table <- .DownloadCellMetadataPerLane(output_file_id)
        list(
          lane_table = lane_table,
          lane_summary = tibble::tibble(
            sourceOutputFileId = output_file_id,
            nCellsIngested = nrow(lane_table)
          ),
          failures = .EmptyFailuresTibble()
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
          )
        )
      }
    )

    failures <- dplyr::bind_rows(failures, lane_result$failures)
    lane_summary_pieces[[length(lane_summary_pieces) + 1]] <- lane_result$lane_summary
    if (!is.null(lane_result$lane_table)) {
      lane_tables[[length(lane_tables) + 1]] <- lane_result$lane_table
    }
  }

  #summarize how many cells were ingested from each output file
  lane_summary <- dplyr::bind_rows(lane_summary_pieces)

  #stack all downloaded lane tables into one table with one row per cell
  if (length(lane_tables)) {
    cell_table <- tibble::as_tibble(do.call(rbind, lane_tables))
    rownames(cell_table) <- NULL
  } else {
    cell_table <- tibble::tibble()
  }

  #override gamma-delta and NK labels in RIRA_TNK using TCR-based classification
  if (classifyTNK && nrow(cell_table) > 0) {
    override_result <- .ApplyTcrTNKOverridesToCellTable(cell_table)
    cell_table <- override_result$cell_table
    failures <- dplyr::bind_rows(failures, override_result$failures)
  }

  #initialize the wide results table; metric columns are added one spec rule at a time
  counts_wide <- tibble::tibble()
  target_fields_seen <- character(0)

  if (nrow(spec_table) > 0) {
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
    }

    #metadata columns that define each group, shared across all spec rules
    grouping_columns <- trimws(strsplit(as.character(grouping_specs[[1]]), "|", fixed = TRUE)[[1]])
    grouping_columns <- grouping_columns[nzchar(grouping_columns)]
    if (!mixed_grouping) {
      if (nrow(cell_table) > 0 && length(grouping_columns)) {
        #start the wide result with one row per distinct grouping key before joining metrics
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

      #iterate through each spec rule, validate columns, compute metrics, and join onto the wide table
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
          .QuantifyQuantileColumns(target_field)
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

        metrics_table <- .QuantifySpecRow(
          spec_row = spec_row,
          cell_table = cell_table,
          grouping_columns = grouping_columns
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
    }
  }

  #write any validation and processing failures to the failure log file
  .WriteFailureLog(failures, failure_log_file = failureLogFile)
  if (nrow(failures) > 0) {
    warning(
      nrow(failures),
      " quantification failure(s) recorded; see failures and ",
      failureLogFile
    )
  }

  list(
    countsWide = counts_wide,
    failures = failures,
    laneSummary = lane_summary
  )
}

## Helper Functions ##
.QuantifyNameDelim <- "__"

.QuantifyQuantileColumns <- function(target_field) {
  paste0(target_field, .QuantifyNameDelim, c("p05", "median", "p95"))
}

.EmptyFailuresTibble <- function() {
  #returns an empty tibble with columns for spec row, output file, field name, and reason
  tibble::tibble(
    specRow = integer(0),
    outputFileId = integer(0),
    field = character(0),
    reason = character(0)
  )
}

#assumes output_file_id refers to a valid LabKey output file row
#assumes DownloadMetadataForSeuratObject returns one row per cell for that file
.DownloadCellMetadataPerLane <- function(output_file_id) {
  df <- DownloadMetadataForSeuratObject(
    outputFileId = output_file_id,
    outFile = tempfile(fileext = ".txt.gz"),
    returnDataFrame = TRUE
  )
  df$sourceOutputFileId <- output_file_id
  tibble::as_tibble(df)
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

  #column fallbacks mirror ClassifyTNKByExpression metadata fields;
  # upgrade path is shared helper if a third caller appears.
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

  cell_table[[paste0(tnk_field, "_rira")]] <- as.character(cell_table[[tnk_field]])
  tnk_scope <- as.character(cell_table[[immune_field]]) == "T_NK"
  updated_tnk <- as.character(cell_table[[tnk_field]])

  updated_tnk[tnk_scope & is_gamma_delta] <- "Gamma/Delta Cells"
  updated_tnk[tnk_scope & is_nk_cell] <- "NK Cells"
  cell_table[[tnk_field]] <- updated_tnk

  list(cell_table = cell_table, failures = failures)
}
#assumes cell_table column names match the downloaded metadata when cells were ingested
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
  } else if (!quantification_type %in% c("sum", "score")) {
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

.ApplySubsetScoreFilter <- function(cells, spec_row) {
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

#this function does the bulk of the work for populating the output table.
#assumes the rule passed validation and grouping columns are present in cell_table
#returns one row per group with the target metric column(s) for this rule
.QuantifySpecRow <- function(spec_row, cell_table, grouping_columns) {
  target_field <- as.character(spec_row$TargetField[[1]])
  source_field <- as.character(spec_row$SourceField[[1]])
  match_value <- as.character(spec_row$MatchValue[[1]])
  quantification_type <- tolower(trimws(as.character(spec_row$QuantificationType[[1]])))

  if (nrow(cell_table) == 0 || !length(grouping_columns)) {
    return(tibble::tibble())
  }

  all_groups <- cell_table %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(grouping_columns)))

  scoped_cells <- cell_table
  scope_field <- trimws(as.character(spec_row$ScopeField[[1]]))
  scope_match_value <- as.character(spec_row$ScopeMatchValue[[1]])
  if (nzchar(scope_field)) {
    scope_values <- as.character(scoped_cells[[scope_field]])
    scoped_cells <- scoped_cells[scope_values == scope_match_value, , drop = FALSE]
  }

  source_values <- as.character(scoped_cells[[source_field]])
  matched_cells <- scoped_cells[source_values == match_value, , drop = FALSE]
  matched_cells <- .ApplySubsetScoreFilter(matched_cells, spec_row)

  if (identical(quantification_type, "sum")) {
    if (nrow(matched_cells) == 0) {
      count_values <- rep(0L, nrow(all_groups))
    } else {
      count_values <- matched_cells %>%
        dplyr::count(dplyr::across(dplyr::all_of(grouping_columns)), name = target_field) %>%
        dplyr::right_join(all_groups, by = grouping_columns) %>%
        dplyr::mutate(!!target_field := ifelse(is.na(.data[[target_field]]), 0L, .data[[target_field]])) %>%
        dplyr::arrange(dplyr::across(dplyr::all_of(grouping_columns)))
      return(count_values)
    }

    all_groups[[target_field]] <- as.integer(count_values)
    return(all_groups)
  }

  score_field <- as.character(spec_row$QuantificationSourceField[[1]])
  quantile_columns <- .QuantifyQuantileColumns(target_field)
  p05_col <- quantile_columns[[1]]
  median_col <- quantile_columns[[2]]
  p95_col <- quantile_columns[[3]]

  if (nrow(matched_cells) == 0) {
    all_groups[[p05_col]] <- NA_real_
    all_groups[[median_col]] <- NA_real_
    all_groups[[p95_col]] <- NA_real_
    return(all_groups)
  }

  score_values <- suppressWarnings(as.numeric(matched_cells[[score_field]]))
  matched_cells[[score_field]] <- score_values

  quantile_metrics <- matched_cells %>%
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
