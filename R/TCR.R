#' @include LabKeySettings.R
#' @include Utils.R
#' @import utils

utils::globalVariables(
  names = c('sortOrder', 'SampleName', 'SubjectId', 'c_gene', 'cdna', 'count', 'd_gene', 'j_gene', 'population', 'raw_clonotype_id', 'raw_consensus_id', 'v_gene'),
  package = 'Rdiscvr',
  add = TRUE
)

#' @title DownloadAndAppendTcrClonotypes
#' @description Download And Append TCR Clonotypes data from Prime-Seq
#' @param seuratObject A Seurat object
#' @param outPath The output filepath
#' @param dropExisting If true, any existing clonotype data will be replaced
#' @param overwriteTcrTable If true, any existing table(s) of TCR clones will be overwritten and re-downloaded
#' @param allowMissing If true, samples missing data will be skipped. Otherwise, the function will fail.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendTcrClonotypes <- function(seuratObject, outPath = tempdir(), dropExisting = T, overwriteTcrTable = F, allowMissing = FALSE){
  if (all(is.null(seuratObject[['BarcodePrefix']]))){
    stop('Seurat object lacks BarcodePrefix column')
  }

  i <- 0
  allPrefixes <- unique(unlist(seuratObject[['BarcodePrefix']]))
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
        stop(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
      }
    }

    clonotypeFile <- file.path(outPath, paste0(barcodePrefix, '.', vloupeId, '.clonotypes.csv'))
    .DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = clonotypeFile, overwrite = overwriteTcrTable)
    if (!file.exists(clonotypeFile)){
      stop(paste0('Unable to download clonotype file for prefix: ', barcodePrefix))
    }

    doDropExisting <- i == 1 && dropExisting
    seuratObject <- .AppendTcrClonotypes(seuratObject, clonotypeFile, barcodePrefix = barcodePrefix, dropExisting = doDropExisting)
  }

  return(seuratObject)
}

.AppendTcrClonotypes <- function(seuratObject = NA, clonotypeFile = NA, barcodePrefix = NULL, dropExisting = F){
  tcr <- .ProcessAndAggregateTcrClonotypes(clonotypeFile)
  if (!is.null(barcodePrefix)){
    tcr$barcode <- as.character(tcr$barcode)
    tcr$barcode <- paste0(barcodePrefix, '_', tcr$barcode)
    tcr$barcode <- as.factor(tcr$barcode)
  }

  origRows <- nrow(tcr)

  datasetSelect <- seuratObject$BarcodePrefix == barcodePrefix
  gexBarcodes <- colnames(seuratObject)[datasetSelect]

  tcrIntersect <- tcr[tcr$barcode %in% gexBarcodes,]
  pct <- round(nrow(tcrIntersect) / origRows * 100, 2)
  pct2 <- round(nrow(tcrIntersect) / length(gexBarcodes) * 100, 2)

  print(paste0('Barcodes with clonotypes: ', origRows, ', intersecting with GEX data (total ', length(gexBarcodes),'): ', nrow(tcrIntersect), " (", pct, "% of TCR / ", pct2, "% of GEX)"))
  if (nrow(tcrIntersect) == 0) {
    print('no barcodes shared')
    print(paste0('first GEX barcodes:'))
    print(head(gexBarcodes))
    print(paste0('first TCR barcodes:'))
    print(head(tcr$barcode))
  }

  tcr <- tcrIntersect

  merged <- merge(data.frame(barcode = gexBarcodes, sortOrder = 1:length(gexBarcodes)), tcr, by = c('barcode'), all.x = T)
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

    if ((colName %in% names(seuratObject@meta.data)) && dropExisting) {
      seuratObject@meta.data[colName] <- NULL
    }

    # Handle legacy columns:
    if (grepl(pattern = '_', x = colName)) {
      legacyName <- gsub(colName, pattern = '_', replacement = '')
      if (legacyName %in% names(seuratObject@meta.data)) {
        print(paste0('Dropping legacy column: ', legacyName))
        seuratObject@meta.data[[legacyName]] <- NULL
      }
    }

    if (!(colName %in% names(seuratObject@meta.data))) {
      toUpdate <- rep(NA, ncol(seuratObject))
    } else {
      toUpdate <- unlist(seuratObject[[colName]])
    }

    # Convert to string in case levels do not match:
    if (is.factor(toUpdate)) {
      toUpdate <- as.character(toUpdate)
    }

    names(toUpdate) <- colnames(seuratObject)
    toUpdate[datasetSelect] <- toAdd
    seuratObject[[colName]] <- as.factor(toUpdate)
  }

  seuratObject$HasCDR3Data <- !is.na(seuratObject$CDR3s)
  seuratObject$HasTRAorB <- !is.na(seuratObject$TRA) | !is.na(seuratObject$TRB)
  seuratObject$HasTRGorD <- !is.na(seuratObject$TRG) | !is.na(seuratObject$TRD)

  return(seuratObject)
}

.FindMatchedVloupe <- function(loupeDataId) {
  rows <- labkey.selectRows(
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
  )

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
  )

  if (nrow(rows) > 1){
    print(paste0('more than one matching VLoupe file found, using most recent: ', tcrReadsetId))
    rows <- rows[1,,drop = F]
  } else if (nrow(rows) == 0) {
    print(paste0('Vloupe file not found for TCR readset: ', tcrReadsetId))
    return(NA)
  }

  return(rows$rowid[1])
}

.DownloadCellRangerClonotypes <- function(vLoupeId, outFile, overwrite = T, fileName = 'all_contig_annotations.csv') {
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

.ProcessTcrClonotypes <- function(clonotypeFile){
  tcr <- utils::read.table(clonotypeFile, header=T, sep = ',', fill = TRUE)
  tcr <- tcr[tcr$cdr3 != 'None' & tcr$cdr3 != '',]

  # drop cellranger '-1' suffix
  tcr$barcode <- gsub("-1", "", tcr$barcode)

  # Many TRDV genes can be used as either alpha or delta TCRs.  10x classifies and TRDV/TRAJ/TRAC clones as 'Multi'.  Re-classify these:
  tcr$chain[tcr$chain == 'Multi' & grepl(pattern = 'TRD', x = tcr$v_gene) & grepl(pattern = 'TRAJ', x = tcr$j_gene) & grepl(pattern = 'TRAC', x = tcr$c_gene)] <- c('TRA')

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

  # Add chain-specific columns:
  tcr$ChainCDR3s <- paste0(tcr$chain, ':', tcr$cdr3)
  for (l in c('TRA', 'TRB', 'TRD', 'TRG')){
    tcr[[l]] <- c(NA)
    tcr[[l]][tcr$chain == l] <- as.character(tcr$cdr3[tcr$chain == l])

    target <- paste0(l, 'V')
    tcr[[target]] <- c(NA)
    tcr[[target]][tcr$chain == l] <- as.character(tcr$v_gene[tcr$chain == l])

    target <- paste0(l, 'J')
    tcr[[target]] <- c(NA)
    tcr[[target]][tcr$chain == l] <- as.character(tcr$j_gene[tcr$chain == l])

    if (l %in% c('TRB', 'TRD')) {
      target <- paste0(l, 'D')
      tcr[[target]] <- c(NA)
      tcr[[target]][tcr$chain == l] <- as.character(tcr$d_gene[tcr$chain == l])
    }

    target <- paste0(l, 'C')
    tcr[[target]] <- c(NA)
    tcr[[target]][tcr$chain == l] <- as.character(tcr$c_gene[tcr$chain == l])
  }

  return(tcr)
}

#' @import Rlabkey
#' @importFrom dplyr %>% coalesce group_by summarise
#' @importFrom naturalsort naturalsort
.ProcessAndAggregateTcrClonotypes <- function(clonotypeFile){
  tcr <- .ProcessTcrClonotypes(clonotypeFile)

  # Summarise, grouping by barcode
  tcr <- tcr %>% group_by(barcode) %>% summarise(
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
