#' @include LabKeySettings.R
#' @include Utils.R
#' @import utils

utils::globalVariables(
  names = c('sortOrder', 'SampleName', 'SubjectId', 'c_gene', 'cdna', 'count', 'd_gene', 'j_gene', 'population', 'raw_clonotype_id', 'raw_consensus_id', 'v_gene', 'CellBarcode', 'IsActive', 'IsActiveLabel', 'GroupField', 'Fraction', 'Label', 'IsShared', 'Tissue', 'metricName',
            'TRA_Segments', 'TRB_Segments', 'TRD_Segments', 'TRG_Segments',
  'TRA_WithProductive', 'TRB_WithProductive', 'TRD_WithProductive', 'TRG_WithProductive'),
  package = 'Rdiscvr',
  add = TRUE
)

#' @title DownloadAndAppendTcrClonotypes
#' @description Download And Append TCR Clonotypes data from Prime-Seq
#' @param seuratObj A Seurat object
#' @param outPath The output filepath
#' @param dropExisting If true, any existing clonotype data will be replaced
#' @param overwriteTcrTable If true, any existing table(s) of TCR clones will be overwritten and re-downloaded
#' @param allowMissing If true, samples missing data will be skipped. Otherwise, the function will fail.
#' @param dropConflictingVJSegments If true, any TRB rows with a TRA/D V/J segments will as dropped, as are TRA rows with TRB/G segments
#' @importFrom magrittr %>%
#' @return A modified Seurat object.
#' @export
DownloadAndAppendTcrClonotypes <- function(seuratObj, outPath = tempdir(), dropExisting = T, overwriteTcrTable = F, allowMissing = FALSE, dropConflictingVJSegments = TRUE){
  if (!'BarcodePrefix' %in% names(seuratObj@meta.data) || any(is.na(seuratObj$BarcodePrefix))){
    seuratObj <- .AttemptToRestoreBarcodePrefix(seuratObj)
  }

  if (!'BarcodePrefix' %in% names(seuratObj@meta.data) || any(is.na(seuratObj$BarcodePrefix))){
    stop('Seurat object lacks BarcodePrefix column or has missing values')
  }

  i <- 0
  allPrefixes <- unique(unlist(seuratObj[['BarcodePrefix']]))
  hasTcrs <- FALSE
  for (barcodePrefix in allPrefixes) {
    i <- i + 1
    print(paste0('Adding TCR clonotypes for prefix: ', barcodePrefix, '. ', i, ' of ', length(allPrefixes)))

    vloupeId <- .FindMatchedVloupe(barcodePrefix)
    if (is.na(vloupeId)){
      hasTcr <- .HasTcrLibrary(barcodePrefix)
      if (is.null(hasTcr)) {
        stop(paste0('Unable to determine if GEX library has a TCR library: ', barcodePrefix))
      } else if (!hasTcr) {
        message(paste0('GEX library lacks a TCR library, skipping: ', barcodePrefix))
        next
      }

      if (allowMissing) {
        warning(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
        next
      } else {
        stop(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
      }
    }

    clonotypeFile <- file.path(outPath, paste0(barcodePrefix, '.', vloupeId, '.clonotypes.csv'))
    .DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = clonotypeFile, overwrite = overwriteTcrTable)
    if (!file.exists(clonotypeFile)){
      stop(paste0('Unable to download clonotype file for prefix: ', barcodePrefix))
    }

    doDropExisting <- i == 1 && dropExisting
    hasTcrs <- TRUE
    seuratObj <- .AppendTcrClonotypes(seuratObj, clonotypeFile, barcodePrefix = barcodePrefix, dropExisting = doDropExisting, dropConflictingVJSegments)
  }

  if (hasTcrs == FALSE) {
    seuratObj$HasCDR3Data <- FALSE
  }

  return(seuratObj)
}

.HasTcrLibrary <- function(loupeDataId, allowNonBlankStatus = FALSE) {
  loupeRows <- suppressWarnings(labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    colSelect="readset",
    colFilter=makeFilter(c("rowid", "EQUAL", loupeDataId)),
    containerFilter=NULL,
    colNameOpt="rname"
  ))

  if (nrow(loupeRows) == 0) {
    return(NULL)
  }

  rs <- loupeRows$readset[1]
  tcrRows <- suppressWarnings(labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="singlecell",
    queryName="cdna_libraries",
    colSelect="tcrreadsetid,tcrreadsetid/status",
    colFilter=makeFilter(c("readsetId", "EQUAL", rs)),
    containerFilter=NULL,
    colNameOpt="rname"
  ))

  foundLibrary <- nrow(tcrRows) > 0
  tcrRows <- tcrRows[!is.na(tcrRows$tcrreadsetid),]

  if (foundLibrary && nrow(tcrRows) == 0) {
    return(FALSE)
  }

  tcrRowsPassing <- tcrRows[is.na(tcrRows$tcrreadsetid_status),]
  if (allowNonBlankStatus || nrow(tcrRowsPassing) > 0) {
    return(TRUE)
  }

  return(FALSE)
}

#' @title CreateMergedTcrClonotypeFile
#' @description Download TCR Clonotypes for all datasets in a seuratObj, update their cellbarcodes with barcodePrefix, and create a merged table
#' @param seuratObj A Seurat object
#' @param outputFile The path where the merged CSV will be written
#' @param overwriteTcrTable If true, any existing table(s) of TCR clones will be overwritten and re-downloaded
#' @param downloadPath The output filepath for per-dataset files
#' @param allowMissing If true, samples missing data will be skipped. Otherwise, the function will fail.
#' @param cellRangerType The type of cellranger data to download. all_contig_annotations_combined is a special case for data processed by prime-seq, which merged the a/b and g/d calls. filtered_contig_annotations_combined.csv is a second option.
#' @param dropConflictingVJSegments If true, any TRB rows with a TRA/D V/J segments will as dropped, as are TRA rows with TRB/G segments
#' @export
CreateMergedTcrClonotypeFile <- function(seuratObj, outputFile, overwriteTcrTable = F, downloadPath = tempdir(), allowMissing = FALSE, cellRangerType = 'all_contig_annotations_combined.csv', dropConflictingVJSegments = TRUE){
  if (all(is.null(seuratObj[['BarcodePrefix']]))){
    stop('Seurat object lacks BarcodePrefix column')
  }

  i <- 0
  allPrefixes <- unique(unlist(seuratObj[['BarcodePrefix']]))
  for (barcodePrefix in allPrefixes) {
    i <- i + 1
    print(paste0('Adding TCR clonotypes for prefix: ', barcodePrefix, '. ', i, ' of ', length(allPrefixes)))

    vloupeId <- .FindMatchedVloupe(barcodePrefix)
    if (is.na(vloupeId)){
      if (allowMissing) {
        warning(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
        next
      } else {
        stop(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
      }
    }

    clonotypeFile <- file.path(downloadPath, paste0(barcodePrefix, '.', vloupeId, '.clonotypes.csv'))
    .DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = clonotypeFile, overwrite = overwriteTcrTable, fileName = cellRangerType)
    if (!file.exists(clonotypeFile)){
      stop(paste0('Unable to download clonotype file for prefix: ', barcodePrefix))
    }

    # Read input, update barcodes:
    dat <- read.table(clonotypeFile, header = T, sep = ',')
    dat$barcode <- gsub("-1", "", dat$barcode)
    dat$barcode <- paste0(barcodePrefix, '_', dat$barcode)
    dat$raw_clonotype_id <- ifelse(dat$is_cell == "true" & dat$productive == "true", yes = paste0(barcodePrefix,'_', dat$raw_clonotype_id), no = "")

    # Check for TRA/B segments that dont match the chain:
    if (dropConflictingVJSegments) {
      dat <- .DropConflictingVJSegments(dat)
    }

    write.table(dat,
                file = outputFile,
                append = i != 1,
                sep = ",",
                row.names = F,
                col.names = i == 1
    )
  }
}

.DropConflictingVJSegments <- function(dat) {
  sel <- dat$chain == 'TRA' & grepl(dat$v_gene, pattern = 'BV|GV')
  if (sum(sel) > 0) {
    print(paste0('Dropping TRA rows with a BV/GV genes, total: ', sum(sel)))
    dat <- dat[!sel,]
  }

  sel <- dat$chain == 'TRA' & grepl(dat$j_gene, pattern = 'BJ|GJ')
  if (sum(sel) > 0) {
    print(paste0('Dropping TRA rows with a BJ/GJ genes, total: ', sum(sel)))
    dat <- dat[!sel,]
  }

  sel <- dat$chain == 'TRB' & grepl(dat$v_gene, pattern = 'AV|DV')
  if (sum(sel) > 0) {
    print(paste0('Dropping TRB rows with a AV/DV genes, total: ', sum(sel)))
    dat <- dat[!sel,]
  }

  sel <- dat$chain == 'TRB' & grepl(dat$j_gene, pattern = 'AJ|DJ')
  if (sum(sel) > 0) {
    print(paste0('Dropping TRB rows with a AJ/DJ genes, total: ', sum(sel)))
    dat <- dat[!sel,]
  }

  return(dat)
}
.AppendTcrClonotypes <- function(seuratObj = NA, clonotypeFile = NA, barcodePrefix = NULL, dropExisting = F, dropConflictingVJSegments = FALSE){
  tcr <- .ProcessAndAggregateTcrClonotypes(clonotypeFile, dropConflictingVJSegments)
  if (!is.null(barcodePrefix)){
    tcr$barcode <- as.character(tcr$barcode)
    tcr$barcode <- paste0(barcodePrefix, '_', tcr$barcode)
    tcr$barcode <- as.factor(tcr$barcode)
    tcr$raw_clonotype_id <- as.character(tcr$raw_clonotype_id)
    tcr$raw_clonotype_id <- ifelse(!is.na(tcr$raw_clonotype_id), paste0(barcodePrefix, '_', tcr$raw_clonotype_id), NA)
    tcr$raw_clonotype_id <- as.factor(tcr$raw_clonotype_id)
  }

  origRows <- nrow(tcr)

  datasetSelect <- seuratObj$BarcodePrefix == barcodePrefix
  gexBarcodes <- colnames(seuratObj)[datasetSelect]

  tcrIntersect <- tcr[tcr$barcode %in% gexBarcodes,]
  pct <- round(nrow(tcrIntersect) / origRows * 100, 2)
  pct2 <- round(nrow(tcrIntersect) / length(gexBarcodes) * 100, 2)

  print(paste0('Barcodes with clonotypes: ', origRows, ', intersecting with GEX data (total ', length(gexBarcodes),'): ', nrow(tcrIntersect), " (", pct, "% of TCR / ", pct2, "% of GEX)"))
  print(paste0('Fraction of TCR records with a TRA: ', round(sum(!is.na(tcr$TRA)) / nrow(tcr), 2)))
  print(paste0('Fraction of TCR records with a TRB: ', round(sum(!is.na(tcr$TRB)) / nrow(tcr), 2)))
  print(paste0('Fraction of TCR records with a TRG: ', round(sum(!is.na(tcr$TRG)) / nrow(tcr), 2)))
  print(paste0('Fraction of TCR records with a TRD: ', round(sum(!is.na(tcr$TRD)) / nrow(tcr), 2)))

  if (nrow(tcrIntersect) == 0) {
    print('no barcodes shared')
    print(paste0('first GEX barcodes:'))
    print(head(gexBarcodes))
    print(paste0('first TCR barcodes:'))
    print(head(tcr$barcode))
  }

  tcr <- tcrIntersect

  merged <- merge(data.frame(barcode = gexBarcodes, sortOrder = seq_along(gexBarcodes)), tcr, by = 'barcode', all.x = T)
  rownames(merged) <- merged$barcode
  merged <- dplyr::arrange(merged, sortOrder)
  merged <- merged[colnames(merged) != 'sortOrder']

  # Check barcodes match before merge
  if (sum(merged$barcode != gexBarcodes) > 0) {
    stop(paste0('Seurat and TCR barcodes do not match after merge, total different: ', sum(merged$barcode != gexBarcodes )))
  }

  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    toAdd <- as.character(merged[[colName]])
    names(toAdd) <- merged[['barcode']]

    if ((colName %in% names(seuratObj@meta.data)) && dropExisting) {
      seuratObj@meta.data[colName] <- NULL
    }

    # Handle legacy columns:
    if (grepl(pattern = '_', x = colName)) {
      legacyName <- gsub(colName, pattern = '_', replacement = '')
      if (legacyName %in% names(seuratObj@meta.data)) {
        print(paste0('Dropping legacy column: ', legacyName))
        seuratObj@meta.data[[legacyName]] <- NULL
      }
    }

    if (!(colName %in% names(seuratObj@meta.data))) {
      toUpdate <- rep(NA, ncol(seuratObj))
    } else {
      toUpdate <- unlist(seuratObj[[colName]])
    }

    # Convert to string in case levels do not match:
    if (is.factor(toUpdate)) {
      toUpdate <- as.character(toUpdate)
    }

    names(toUpdate) <- colnames(seuratObj)
    toUpdate[datasetSelect] <- toAdd
    seuratObj[[colName]] <- as.factor(toUpdate)
  }

  seuratObj$HasCDR3Data <- !is.na(seuratObj$CDR3s)
  seuratObj$HasTRAorB <- !is.na(seuratObj$TRA) | !is.na(seuratObj$TRB)
  seuratObj$HasTRGorD <- !is.na(seuratObj$TRG) | !is.na(seuratObj$TRD)

  return(seuratObj)
}

.FindMatchedVloupe <- function(loupeDataId) {
  rows <- suppressWarnings(labkey.selectRows(
		baseUrl=.getBaseUrl(),
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		viewName="",
		colSort="-rowid",
		colSelect="readset/cdna/tcrReadsetId",
		colFilter=makeFilter(c("rowid", "EQUAL", loupeDataId)),
		containerFilter=NULL,
		colNameOpt="rname"
  ))

  if (nrow(rows) == 0) {
    translatedRows <- .ResolveLoupeIdFromDeleted(loupeDataId, throwOnError = FALSE)
    if (all(is.null(translatedRows))) {
      print(paste0("Loupe File ID: ", loupeDataId, " not found"))
      return(NA)
    }

    if (nrow(translatedRows) > 0) {
      loupeDataId <- unique(translatedRows$loupefileid)
      if (length(loupeDataId) > 1) {
        loupeDataId <- max(loupeDataId)
        print(paste0('more than one matching loupeDataId found, using most recent: ', loupeDataId))
      }

      gexRows <- suppressWarnings(labkey.selectRows(
        baseUrl=.getBaseUrl(),
        folderPath=.getLabKeyDefaultFolder(),
        schemaName="singlecell",
        queryName="singlecellDatasets",
        colSelect="readsetId",
        colFilter=makeFilter(
          c("loupeFileId", "EQUAL", loupeDataId)
        ),
        containerFilter=NULL,
        colNameOpt="rname"
      )) |> unique()

      if (nrow(gexRows) > 0) {
        gexReadset <- sort(unique(gexRows$readsetid), decreasing = TRUE)[1]

        cdnaRows <- suppressWarnings(labkey.selectRows(
          baseUrl=.getBaseUrl(),
          folderPath=.getLabKeyDefaultFolder(),
          schemaName="singlecell",
          queryName="cdna_libraries",
          colSelect="tcrReadsetId",
          colFilter=makeFilter(
            c("readsetid", "EQUAL", gexReadset)
          ),
          containerFilter=NULL,
          colNameOpt="rname"
        )) |> unique()

        if (nrow(gexRows) > 0) {
          tcrReadsets <- sort(unique(cdnaRows$tcrreadsetid), decreasing = TRUE)
          if (length(tcrReadsets) > 1) {
            stop(paste0('More than one TCR readset found associated with GEX readset: ', gexReadset, ' and barcode: ', loupeDataId))
          }

          print(paste0('Found TCR readset in deleted loupe files table: ', tcrReadsets[1]))
          rows <- cdnaRows
          names(rows) <- 'readset_cdna_tcrreadsetid'
        }
      }
    }
  }

  if (nrow(rows) != 1) {
    return(NA)
  }

  tcrReadsetId <- rows[['readset_cdna_tcrreadsetid']]
  if (is.na(tcrReadsetId) || is.null(tcrReadsetId)) {
    return(NA)
  }

  rows <- labkey.selectRows(
		baseUrl=.getBaseUrl(),
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		viewName="",
		colSort="-rowid",
		colSelect="rowid,",
		colFilter=makeFilter(c("readset", "EQUAL", tcrReadsetId), c("category", "EQUAL", "10x VLoupe")),
		containerFilter=NULL,
		colNameOpt="rname"
  ) |> unique()

  if (nrow(rows) > 1){
    print(paste0('more than one matching VLoupe file found, using most recent for TCR readset: ', tcrReadsetId))
    rows <- rows[1,,drop = F]
  } else if (nrow(rows) == 0) {
    print(paste0('Vloupe file not found for TCR readset: ', tcrReadsetId))
    return(NA)
  }

  return(rows$rowid[1])
}

.DownloadCellRangerClonotypes <- function(vLoupeId, outFile, overwrite = T, fileName = 'all_contig_annotations_combined.csv') {
  #The file will be in the same directory as the VLoupe file
  rows <- labkey.selectRows(
		baseUrl=.getBaseUrl(),
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		viewName="",
		colSort="-rowid",
		colSelect="rowid,workbook/workbookid,dataid/webdavurlrelative",
		colFilter=makeFilter(c("rowid", "EQUAL", vLoupeId)),
		containerFilter=NULL,
		colNameOpt="rname"
  )

  if (nrow(rows) != 1) {
    return(NA)
  }

  wb <- rows[['workbook_workbookid']]
  if (is.na(wb) || is.null(wb)){
    wb <- ''
  }

  remotePath <- rows[['dataid_webdavurlrelative']]
  remotePath <- paste0(dirname(remotePath), '/', fileName)

  fe <- labkey.webdav.pathExists(
    baseUrl=.getBaseUrl(),
    folderPath=paste0(.getLabKeyDefaultFolder(),wb),
    remoteFilePath = remotePath
  )

  if (!fe) {
    stop(paste0('Remote file does not exist: ', remotePath))
  }

  success <- labkey.webdav.get(
		baseUrl=.getBaseUrl(),
		folderPath=paste0(.getLabKeyDefaultFolder(),wb),
		remoteFilePath = remotePath,
		overwrite = overwrite,
		localFilePath = outFile
  )

  if (!success | !file.exists(outFile)) {
    return(NA)
  }

  return(outFile)
}

utils::globalVariables(
	names = c('chain', 'cdr3', 'LabelCol', 'barcode', 'ChainCDR3s', 'TRA', 'TRB', 'TRD', 'TRG', 'TRAV', 'TRBV', 'TRDV', 'TRGV', 'CloneName'),
	package = 'Rdiscvr',
	add = TRUE
)

.ProcessTcrClonotypes <- function(clonotypeFile, dropConflictingVJSegments = FALSE){
  tcr <- utils::read.table(clonotypeFile, header=T, sep = ',', fill = TRUE)
  tcr <- tcr[tcr$cdr3 != 'None' & tcr$cdr3 != '',]

  # drop cellranger '-1' suffix
  tcr$barcode <- gsub("-1", "", tcr$barcode)

  # Many TRDV genes can be used as either alpha or delta TCRs.  10x classifies and TRDV/TRAJ/TRAC clones as 'Multi'.  Re-classify these:
  tcr$chain[tcr$chain == 'Multi' & grepl(pattern = 'TRD', x = tcr$v_gene) & grepl(pattern = 'TRAJ', x = tcr$j_gene) & grepl(pattern = 'TRAC', x = tcr$c_gene)] <- 'TRA'

  if (dropConflictingVJSegments) {
    tcr <- .DropConflictingVJSegments(tcr)
  }

  tcr$cdr3WithSuffix <- paste0(tcr$cdr3, ifelse(tcr$productive == 'True', yes = '', no = ' (NP)'))

  #Download named clonotypes and merge:
  # Add clone names:
  labelDf <- suppressWarnings(labkey.selectRows(
		baseUrl=.getBaseUrl(),
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="tcrdb",
		queryName="clones",
		showHidden=TRUE,
		colSelect=c('clonename','chain','cdr3','displayname'),
		containerFilter=NULL,
		colNameOpt='rname'
  ))

  labelDf$LabelCol <- coalesce(as.character(labelDf$displayname), as.character(labelDf$clonename))

  labelDf <- labelDf %>%
    group_by(chain, cdr3) %>%
    summarise(CloneName = paste0(sort(unique(LabelCol)), collapse = ","))

  tcr <- merge(tcr, labelDf, by.x = c('chain', 'cdr3'), by.y = c('chain', 'cdr3'), all.x = TRUE, all.y = FALSE)

  # TODO: smarter join accounting for A/B and also commas

  # Add chain-specific columns:
  tcr$ChainCDR3s <- paste0(tcr$chain, ':', tcr$cdr3)
  for (l in c('TRA', 'TRB', 'TRD', 'TRG')){
    tcr[[l]] <- NA
    tcr[[l]][tcr$chain == l] <- as.character(tcr$cdr3[tcr$chain == l])

    target <- paste0(l, 'V')
    tcr[[target]] <- NA
    tcr[[target]][tcr$chain == l] <- as.character(tcr$v_gene[tcr$chain == l])

    target <- paste0(l, 'J')
    tcr[[target]] <- NA
    tcr[[target]][tcr$chain == l] <- as.character(tcr$j_gene[tcr$chain == l])

    target <- paste0(l, '_WithProductive')
    tcr[[target]] <- NA
    tcr[[target]][tcr$chain == l] <- as.character(tcr$cdr3WithSuffix[tcr$chain == l])

    if (l %in% c('TRB', 'TRD')) {
      target <- paste0(l, 'D')
      tcr[[target]] <- NA
      tcr[[target]][tcr$chain == l] <- as.character(tcr$d_gene[tcr$chain == l])
    }

    target <- paste0(l, 'C')
    tcr[[target]] <- NA
    tcr[[target]][tcr$chain == l] <- as.character(tcr$c_gene[tcr$chain == l])
  }

  return(tcr)
}

#' @import Rlabkey
#' @importFrom dplyr %>% coalesce group_by summarise
#' @importFrom naturalsort naturalsort
.ProcessAndAggregateTcrClonotypes <- function(clonotypeFile, dropConflictingVJSegments = FALSE){
  tcr <- .ProcessTcrClonotypes(clonotypeFile, dropConflictingVJSegments = dropConflictingVJSegments)

  # Summarise, grouping by barcode
  tcr <- tcr %>%
    mutate(
      TRA_Segments = ifelse(TRA == '', yes = NA, no = paste0(dplyr::coalesce(TRA, '-'), '|', dplyr::coalesce(TRAV, '-'), '|', dplyr::coalesce(TRAJ, '-'))),
      TRB_Segments = ifelse(TRB == '', yes = NA, no = paste0(dplyr::coalesce(TRB, '-'), '|', dplyr::coalesce(TRBV, '-'), '|', dplyr::coalesce(TRBJ, '-'))),
      TRD_Segments = ifelse(TRD == '', yes = NA, no = paste0(dplyr::coalesce(TRD, '-'), '|', dplyr::coalesce(TRDV, '-'), '|', dplyr::coalesce(TRDJ, '-'))),
      TRG_Segments = ifelse(TRG == '', yes = NA, no = paste0(dplyr::coalesce(TRG, '-'), '|', dplyr::coalesce(TRGV, '-'), '|', dplyr::coalesce(TRGJ, '-')))
    ) %>%
    group_by(barcode) %>%
    summarise(
      ChainCDR3s = paste0(sort(unique(ChainCDR3s[ChainCDR3s != ''])), collapse = ","),
      CDR3s = paste0(sort(unique(cdr3[cdr3 != ''])), collapse = ","),
      TRA = paste0(sort(unique(as.character(TRA[TRA != '']))), collapse = ","),
      TRB = paste0(sort(unique(as.character(TRB[TRB != '']))), collapse = ","),
      TRD = paste0(sort(unique(as.character(TRD[TRD != '']))), collapse = ","),
      TRG = paste0(sort(unique(as.character(TRG[TRG != '']))), collapse = ","),
      TRA_V = paste0(sort(unique(as.character(TRAV[TRAV != '']))), collapse = ","),
      TRB_V = paste0(sort(unique(as.character(TRBV[TRBV != '']))), collapse = ","),
      TRD_V = paste0(sort(unique(as.character(TRDV[TRDV != '']))), collapse = ","),
      TRG_V = paste0(sort(unique(as.character(TRGV[TRGV != '']))), collapse = ","),
      TRA_J = paste0(sort(unique(as.character(TRAJ[TRAJ != '']))), collapse = ","),
      TRB_J = paste0(sort(unique(as.character(TRBJ[TRBJ != '']))), collapse = ","),
      TRD_J = paste0(sort(unique(as.character(TRDJ[TRDJ != '']))), collapse = ","),
      TRG_J = paste0(sort(unique(as.character(TRGJ[TRGJ != '']))), collapse = ","),
      TRA_C = paste0(sort(unique(as.character(TRAC[TRAC != '']))), collapse = ","),
      TRB_C = paste0(sort(unique(as.character(TRBC[TRBC != '']))), collapse = ","),
      TRD_C = paste0(sort(unique(as.character(TRDC[TRDC != '']))), collapse = ","),
      TRG_C = paste0(sort(unique(as.character(TRGC[TRGC != '']))), collapse = ","),
      TRB_D = paste0(sort(unique(as.character(TRBD[TRBD != '']))), collapse = ","),
      TRD_D = paste0(sort(unique(as.character(TRDD[TRDD != '']))), collapse = ","),
      TRA_WithProductive = paste0(sort(unique(as.character(TRA_WithProductive[TRA_WithProductive != '']))), collapse = ","),
      TRB_WithProductive = paste0(sort(unique(as.character(TRB_WithProductive[TRB_WithProductive != '']))), collapse = ","),
      TRD_WithProductive = paste0(sort(unique(as.character(TRD_WithProductive[TRD_WithProductive != '']))), collapse = ","),
      TRG_WithProductive = paste0(sort(unique(as.character(TRG_WithProductive[TRG_WithProductive != '']))), collapse = ","),
      TRA_Segments = paste0(sort(unique(as.character(TRA_Segments))), collapse = ","),
      TRB_Segments = paste0(sort(unique(as.character(TRB_Segments))), collapse = ","),
      TRD_Segments = paste0(sort(unique(as.character(TRD_Segments))), collapse = ","),
      TRG_Segments = paste0(sort(unique(as.character(TRG_Segments))), collapse = ","),
      raw_clonotype_id = paste0(sort(unique(as.character(raw_clonotype_id[raw_clonotype_id != '']))), collapse = ","),
      CloneNames = paste0(sort(unique(CloneName[CloneName != ''])), collapse = ",")  #this is imprecise b/c we count a hit if we match any chain, but this is probably what we often want
    )

  # Note: we should attempt to produce a more specfic call, assuming we have data from multiple chains
  # The intent of this was to allow a A- or B-only hit to produce a call, but if we have both A/B, take their intersect.

  tcr$CloneNames <- sapply(strsplit(as.character(tcr$CloneNames), ",", fixed = TRUE), function(x) paste0(naturalsort(unique(x)), collapse = ","))

  tcr$barcode <- as.factor(tcr$barcode)
  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    v <- tcr[[colName]]
    v <- as.character(v)
    v[v == ''] <- NA

    tcr[[colName]] <- as.factor(v)
  }

  return(tcr)
}

#' @title Download10xRawDataForLoupeFile
#' @description Downloads the raw_feature_bc_matrix folder associated with the provided Loupe file
#' @param outputFileId The outputfile Id of the loupe file
#' @param outFile The local path to write this file. The data will be written to a subfolder named raw_feature_bc_matrix
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param countType Either raw_feature_bc_matrix or filtered_feature_bc_matrix
#' @export
#'
#' @import Rlabkey
Download10xRawDataForLoupeFile <- function(outputFileId, outFile, overwrite = T, countType = 'raw_feature_bc_matrix') {
	if (!dir.exists(outFile)) {
    print('Creating output folder')
		dir.create(outFile)
	}

  if (substr(outFile, nchar(outFile), nchar(outFile)) != '/'){
    outFile <- paste0(outFile, '/')
  }

	expectedDir <- paste0(outFile, '/', countType)
	if (!overwrite && dir.exists(outFile)) {
		print('File exists, will not overwrite')
    return(outFile)
	} else if (overwrite && dir.exists(outFile)) {
		print('File exists, deleting')
		unlink(expectedDir, recursive = T)
	}

	return(DownloadOutputDirectoryFromOutputFile(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = function(x){
		return(paste0(dirname(x), '/', countType))
	}))
}

#' @title Run CoNGA
#'
#' @description Runs CoNGA on a seurat object
#' @param seuratObj The Seurat object containing the data to be run using CoNGA.
#' @param organism The organism for the TCR references, either 'human' or 'rhesus'
#' @param seuratToCongaDir The directory to store the results of SeuratToCoNGA() (the input files for the python call to run_CoNGA()).
#' @param assayName Pass-through variable for accessing the assay name within the Seurat object for SeuratToCoNGA().
#' @param runCongaOutputFilePrefix prefix for the output files from the python call to run_CoNGA().
#' @param gexDatatype This should be "10x_h5" since we're using DropletUtils to write the counts, although "h5ad" is supported.
#' @param runCongaOutfilePrefixForQcPlots Prefix for the qc output files to be generated by run_CoNGA().
#' @param runCongaOutputDirectory The directory that will store the many files created during run_CoNGA().
#' @param congaMetadataPrefix A prefix to be added to the columns of the metadata that will be added from CoNGA within the returned Seurat object.
#' @param pngConversionTool specify the conversion tool used by CoNGA for svg to png conversion. Can be any of: convert, inkscape, rsvg
#' @return A Seurat object with the conga metadata relevant to TCR + gene expression clustering appended.
#' @export
RunCoNGA <- function(seuratObj,
                     organism = NULL,
                     seuratToCongaDir = "./seuratToConga",
                     assayName = "RNA",
                     runCongaOutputFilePrefix = "conga_output",
                     gexDatatype = "10x_h5",
                     runCongaOutfilePrefixForQcPlots = "qc_plots",
                     runCongaOutputDirectory = "./conga_output",
                     congaMetadataPrefix = "conga_",
                     pngConversionTool = NULL) {

  clonotypeFile <- tempfile(fileext = '.csv')
  CreateMergedTcrClonotypeFile(seuratObj, outputFile = clonotypeFile, overwriteTcrTable = T)

  return(CellMembrane::RunCoNGA(
    seuratObj = seuratObj,
    tcrClonesFile = clonotypeFile,
    organism = organism,
    seuratToCongaDir = seuratToCongaDir,
    assayName = assayName,
    runCongaOutputFilePrefix = runCongaOutputFilePrefix,
    gexDatatype = gexDatatype,
    runCongaOutfilePrefixForQcPlots = runCongaOutfilePrefixForQcPlots,
    runCongaOutputDirectory = runCongaOutputDirectory,
    congaMetadataPrefix = congaMetadataPrefix,
    pngConversionTool = pngConversionTool
  ))
}

#' @title Classify T and NK By Expression
#' @description Classify T and NK By Expression and TCR clonotype, using best available evidence of ground-truth
#' @param seuratObj The seurat object
#' @param assayName The name of the RNA assay
#' @param constantRegionCountThreshold Any cell with a TCR constant region raw count above this threshold is considered positive for that gene
#' @param includeConstantRegionExpression If includeConstantRegionExpression is true, RNA expression of TRA/B/G constant regions will be considered in the classification. This is disabled by default because we find expression of these genes may be leaky and not necessarily linked to a functional TCR
#' @param includeDeltaConstantRegionExpression Similar to includeConstantRegionExpression, but applies to the delta constant region alone
#' @param collapseGOnlyToGD By default, cells with gamma-chain alone are reported as a separate category (since A/B T cells can encode a gamma-chain). If TRUE, these will be collapsed into the Gamma/Delta category.
#' @import Seurat ggplot2
#' @export
ClassifyTNKByExpression <- function(seuratObj, assayName = 'RNA', constantRegionCountThreshold = 1.5, includeConstantRegionExpression = FALSE, includeDeltaConstantRegionExpression = FALSE, collapseGOnlyToGD = FALSE) {
  if (!'HasCDR3Data' %in% names(seuratObj@meta.data)) {
    # Note: this will error if data are missing, but not if TCR is not used
    seuratObj <- Rdiscvr::DownloadAndAppendTcrClonotypes(seuratObj, allowMissing = FALSE)
  }

  if (any(is.na(seuratObj$HasCDR3Data))) {
    print('Fixing HasCDR3Data field')
    seuratObj$HasCDR3Data <- !is.na(seuratObj$CDR3s)
  }

  if (all(seuratObj$HasCDR3Data == FALSE)) {
    print('None of the cells had CDR3 data, skipping')
    return(seuratObj)
  }

  ad <- Seurat::GetAssayData(seuratObj, layer = 'counts', assay = assayName)

  testGeneGt0 <- function(ad, featureName, defaultValue = FALSE, threshold = 0){
    if (!featureName %in% rownames(ad)){
      return(defaultValue)
    }

    return(as.numeric(ad[featureName,]) > threshold)
  }

  # LOC711031 = TRDC
  seuratObj$IsGammaDelta <- !is.na(seuratObj$TRD) | (includeDeltaConstantRegionExpression & testGeneGt0(ad, 'LOC711031'))
  print(DimPlot(seuratObj, group.by = 'IsGammaDelta'))

  seuratObj$HasCD3 <- testGeneGt0(ad, 'CD3D') | testGeneGt0(ad, 'CD3E') | testGeneGt0(ad, 'CD3G')
  print(DimPlot(seuratObj, group.by = 'HasCD3'))

  seuratObj$HasTCRConstant <- testGeneGt0(ad, 'LOC711031', threshold = constantRegionCountThreshold) |
    testGeneGt0(ad, 'LOC720538', threshold = constantRegionCountThreshold) |
    testGeneGt0(ad, 'LOC705095', threshold = constantRegionCountThreshold) |
    testGeneGt0(ad, 'LOC710951', threshold = constantRegionCountThreshold) |
    testGeneGt0(ad, 'LOC114677140', threshold = constantRegionCountThreshold)

  print(DimPlot(seuratObj, group.by = 'HasTCRConstant'))

  seuratObj$IsNKCell <- !seuratObj$HasCDR3Data & !seuratObj$HasCD3 & (!includeConstantRegionExpression | !seuratObj$HasTCRConstant)
  print(DimPlot(seuratObj, group.by = 'IsNKCell'))

  seuratObj$HasGammaChain <- !is.na(seuratObj$TRG) | (includeConstantRegionExpression & (testGeneGt0(ad, 'LOC720538') | testGeneGt0(ad, 'LOC705095')))
  print(DimPlot(seuratObj, group.by = 'HasGammaChain'))

  # As above, allow either TCR data or constant chaine expression:
  seuratObj$IsAlphaBeta <- FALSE
  seuratObj$IsAlphaBeta[!is.na(seuratObj$TRA) | !is.na(seuratObj$TRB)] <- TRUE
  if (includeConstantRegionExpression) {
    seuratObj$IsAlphaBeta[testGeneGt0(ad, 'LOC710951') | testGeneGt0(ad, 'LOC114677140')] <- TRUE
  }
  print(DimPlot(seuratObj, group.by = 'IsAlphaBeta'))
  
  seuratObj$TNK_Type <- NA
  seuratObj$TNK_Type[seuratObj$IsNKCell] <- 'NK (CD3-/TCR-)'
  seuratObj$TNK_Type[seuratObj$IsGammaDelta] <- 'Gamma/Delta'
  seuratObj$TNK_Type[seuratObj$IsAlphaBeta] <- 'Alpha/Beta'
  seuratObj$TNK_Type[(seuratObj$IsNKCell + seuratObj$IsAlphaBeta + seuratObj$IsGammaDelta) > 1] <- 'Ambiguous'

  # This allows a cell with a gamma chain, but not evidence of A/B to be called as gamma/delta
  if (collapseGOnlyToGD) {
    seuratObj$TNK_Type[is.na(seuratObj$TNK_Type) & seuratObj$HasGammaChain] <-'Gamma/Delta'
  } else {
    seuratObj$TNK_Type[is.na(seuratObj$TNK_Type) & seuratObj$HasGammaChain] <-'Gamma Chain-Only'
  }

  seuratObj$TNK_Type[is.na(seuratObj$TNK_Type)] <- 'Unknown'
  seuratObj$TNK_Type[seuratObj$TNK_Type == 'Unknown' & !seuratObj$HasCDR3Data] <- 'Other (TCR-)'
  seuratObj$TNK_Type <- naturalsort::naturalfactor(seuratObj$TNK_Type)

  print(DimPlot(seuratObj, group.by = 'TNK_Type'))
  print(table(seuratObj$TNK_Type, useNA = 'always'))

  return(seuratObj)
}

#' @title MakeClonotypePlot
#' @description Generates a summary plot of clonotype data
#' @param seuratObj A Seurat object
#' @param outFile The output file path to which results will be written
#' @param xFacetField Passed to facet_grid
#' @param groupingFields The set of fields used for grouping data
#' @param activationFieldName The name of the field holding the score to be used to determine activation state
#' @param threshold The minimum value to consider a cell activated
#' @importFrom magrittr %>%
#' @return A modified Seurat object.
#' @export
SummarizeTNK_Activation <- function(seuratObj, outFile, xFacetField = 'Population', groupingFields = c('Stim', 'Population', 'SampleDate', 'AssayType'), activationFieldName = 'TandNK_Activation_UCell', threshold = 0.5) {
  seuratObj <- QueryAndApplyCdnaMetadata(seuratObj)

  if (!'HasCDR3Data' %in% names(seuratObj@meta.data)) {
    stop('This seurat object appears to be missing TCR data. See RDiscvr::DownloadAndAppendTcrClonotypes')
  }

  seuratObj$IsActive <- seuratObj[[activationFieldName]] >= threshold

  PA <- FeaturePlot(seuratObj, features = activationFieldName, min.cutoff = 'q02', max.cutoff = 'q98') +
    scale_color_gradientn(colors = c("navy", "cadetblue2", "gold", "red")) +
    ggtitle('TCR Signaling Score') +
    labs(x = 'UMAP_1', y = 'UMAP_2')
  print(PA)

  print('Summarizing T Cell activation by subject:')
  results <- NULL
  for (subjectId in sort(unique(seuratObj$SubjectId))) {
    for (chain in c('TRA', 'TRB')) {
      outFileTemp <- tempfile()
      print(MakeClonotypePlot(seuratObj, outFile = outFileTemp, subjectId = subjectId, chain = chain, threshold = threshold, activationFieldName = activationFieldName, groupingFields = groupingFields, xFacetField = xFacetField))
      dat <- read.table(outFileTemp, header = TRUE, sep = '\t')
      unlink(outFileTemp)

      if (all(is.null(results))) {
        results <- dat
      } else {
        results <- rbind(results, dat)
      }
    }

    write.table(results, file = outFile, quote = FALSE, sep = '\t')
  }
}

#' @title MakeClonotypePlot
#' @description Generates a summary plot of clonotype data
#' @param seuratObj A Seurat object
#' @param subjectId The subject Id to show
#' @param outFile An optional file where a TSV of the data will be written
#' @param chain The chain (i.e. TRA, TRB, TRG, TRD)
#' @param xFacetField Passed to facet_grid
#' @param groupingFields The set of fields used for grouping data
#' @param threshold The minimum value to consider a cell activated
#' @param activationFieldName The name of the field holding the score to be used to determine activation state
#' @param lowFreqThreshold Any clone not appearing above this threshold will be marked as low. frequency
#' @importFrom magrittr %>%
#' @return The plot object
#' @export
MakeClonotypePlot <- function(seuratObj, outFile = NULL, subjectId, chain, xFacetField = 'Population', groupingFields = c('Stim', 'Population', 'Tissue', 'SampleDate', 'AssayType'), threshold = 0.5, activationFieldName = 'TandNK_ActivationCore_UCell', lowFreqThreshold = 0.005) {
  dat <- seuratObj@meta.data %>%
    filter(SubjectId == subjectId) %>%
    mutate(IsActive = !!sym(activationFieldName) >= threshold)

  dat$CellBarcode <- rownames(dat)
  dat$CDR3 <- dat[[chain]]

  dat <- dat[!is.na(dat$CDR3),]
  dat <- dat %>%
    group_by(across(all_of(c(groupingFields, 'IsActive')))) %>%
    mutate(TotalCells = n_distinct(CellBarcode))

  dat <- dat %>% group_by(across(all_of(groupingFields))) %>%
    mutate(TotalForGroup = n_distinct(CellBarcode))

  dat <- dat %>% group_by(across(all_of(c(groupingFields, 'IsActive', 'TotalCells', 'TotalForGroup', 'CDR3')))) %>%
    summarise(Count = n())

  dat$Fraction <- dat$Count / dat$TotalForGroup
  dat$Label <- as.character(dat[['CDR3']])
  toKeep <- dat$Label[dat$Fraction > lowFreqThreshold]
  dat$Label[!dat$Label %in% toKeep] <- 'Low Freq'
  rm(toKeep)

  dat <- dat %>%
    group_by(across(all_of(c(groupingFields, 'CDR3')))) %>%
    mutate(TotalSubset = n_distinct(IsActive))

  dat$IsShared <- dat$TotalSubset > 1
  dat$IsShared[dat$Label == 'Low Freq'] <- FALSE
  dat$IsShared <- ifelse(dat$IsShared, yes = 'Yes', no = 'No')

  colorSteps <- max(min(length(unique(dat$Label[dat$Label != 'Low Freq'])), 9), 3)
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(colorSteps, 'Set1'))

  patternValues <- c('stripe', 'none')
  names(patternValues) <- c('Yes', 'No')

  dat$Label <- forcats::fct_reorder(dat$Label, desc(dat$Fraction))
  if ('Low Freq' %in% dat$Label) {
    dat$Label <- forcats::fct_relevel(dat$Label, 'Low Freq', after = 0)
  }

  cols <- getPalette(length(unique(dat$Label[dat$Label != 'Low Freq'])))
  if ('Low Freq' %in% dat$Label) {
    cols <- c('#ECECEC', cols)
  }
  names(cols) <- levels(dat$Label)

  # Prune fields, as possible:
  groupingFieldsOrig <- groupingFields
  for (idx in length(groupingFields):1) {
    if (length(unique(dat[[groupingFields[idx]]])) == 1) {
      print(paste0('Dropping group field with single value: ', groupingFields[idx]))
      groupingFields <- groupingFields[groupingFields != groupingFields[idx]]
    }
  }

  groupingFields <- groupingFields[groupingFields != xFacetField]

  # Restore fields if there is just one group:
  if (length(groupingFields) == 0) {
    groupingFields <- groupingFieldsOrig
  }

  dat <- dat %>% tidyr::unite(col = 'GroupField', {{groupingFields}}, remove = FALSE)
  dat$GroupField <- naturalsort::naturalfactor(dat$GroupField)

  dat$IsActiveLabel <- ifelse(dat$IsActive, yes = 'Activated', no = 'Not Activated')

  wrap_by <- function(xFacetField) {
    facet_grid(vars(IsActiveLabel), vars(!!sym(xFacetField)), scales = 'free_y')
  }

  PT <- ggplot(dat, aes(x = GroupField, y = Fraction, fill = Label, pattern = IsShared)) +
      ggpattern::geom_col_pattern(pattern_fill = "black", color = 'black',
                                  pattern_density = 0.2,
                                  pattern_spacing = 0.05,
                                  pattern_key_scale_factor = 0.6
      ) +
      ggpattern::scale_pattern_manual(values = patternValues) +
      wrap_by(xFacetField) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(labels = scales::percent) +
      labs(y = 'Fraction of Cells', x = '', fill = 'Clone') +
      theme_classic(base_size = 14) +
      theme(
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle(paste0(subjectId, ': ', chain))

  if (!is.null(outFile)) {
    dat$Chain <- chain
    dat$SubjectId <- subjectId

    write.table(dat, sep = '\t', quote = FALSE, row.names = FALSE, file = outFile)
  }

  PT
}

#' @title CalculateAndStoreTcrRepertoireStats
#' @description Calculates  a summary plot of clonotype data
#' @param seuratObj A Seurat object
#' @param outputFile An optional filepath, where the complete results will be written as a TSV
#' @importFrom magrittr %>%
#' @return The dataframe with results
#' @export
CalculateAndStoreTcrRepertoireStats <- function(seuratObj, outputFile = NULL) {
  if (! 'BarcodePrefix' %in% names(seuratObj@meta.data)) {
    stop(paste0('Column BarcodePrefix missing from the seuratObj'))
  }

  if (!'TRA' %in% names(seuratObj@meta.data) || !'TRB' %in% names(seuratObj@meta.data)) {
    prefixes <- unique(seuratObj$BarcodePrefix)
    if (length(prefixes) == 1) {
      hasTcr <- .HasTcrLibrary(prefixes[1])
      if (hasTcr == FALSE) {
        print('The seurat object does not include a VDJ library, skipping CalculateTcrRepertoireStatsByPopulation')
        return(NULL)
      }
    }
  }

  df <- CellMembrane::CalculateTcrRepertoireStatsByPopulation(seuratObj@meta.data, groupField = 'cDNA_ID')
  if (all(is.null(df))) {
    print('No results returned by CalculateTcrRepertoireStatsByPopulation')
    return(NULL)
  }

  if (! 'cDNA_ID' %in% names(df)) {
    stop(paste0('Expected cDNA_ID in column names. Found: ', paste0(names(df), collapse = ',')))
  }

  if (any(is.na(df$cDNA_ID))) {
    stop('Found NA values for cDNA_ID')
  }

  existingCDNA <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="singlecell",
    queryName="cdna_libraries",
    colSelect="rowid,container",
    colFilter=makeFilter(c("rowid", "IN", paste0(unique(df$cDNA_ID), collapse = ';'))),
    colNameOpt="rname"
  )

  for (rowid in unique(df$cDNA_ID)) {
    if (! rowid %in% existingCDNA$rowid) {
      stop(paste0('Unknown cDNA ID: ', rowid))
    }

    containerId <- existingCDNA$container[existingCDNA$rowid == rowid]

    existingRows <- labkey.selectRows(
      baseUrl=.getBaseUrl(),
      folderPath=.getLabKeyDefaultFolder(),
      schemaName="tcrdb",
      queryName="repertoire_stats",
      colSelect="rowid,container",
      colFilter=makeFilter(
        c("cdna_id", "EQUAL", rowid)
      ),
      colNameOpt="rname"
    )

    if (nrow(existingRows) > 0) {
      print(paste0('Deleting ', nrow(existingRows), ' existing rows'))
      deleted <- labkey.deleteRows(
        baseUrl=.getBaseUrl(),
        folderPath=.getLabKeyDefaultFolder(),
        schemaName="tcrdb",
        queryName="repertoire_stats",
        toDelete = existingRows
      )
    }

    toInsert <- df %>%
      filter(cDNA_ID == rowid) %>%
      mutate(container = containerId) %>%
      rename(
        cdna_id = 'cDNA_ID',
        metricName = 'MetricName',
        value = 'Value',
        cellType = 'population'
      )

    if (!is.null(outputFile)) {
      write.table(toInsert, file = outputFile, sep = '\t', quote = FALSE, row.names = FALSE)
    }

    # Added to limit the number of 'TopX' record we import:
    toInsert <- toInsert %>%
      filter(!grepl(metricName, pattern = '_Top') | grepl(metricName, pattern = '_Top5$') | grepl(metricName, pattern = '_Top10$'))

    inserted <- labkey.insertRows(
      baseUrl=.getBaseUrl(),
      folderPath=.getLabKeyDefaultFolder(),
      schemaName="tcrdb",
      queryName="repertoire_stats",
      toInsert = toInsert
    )
  }
}
