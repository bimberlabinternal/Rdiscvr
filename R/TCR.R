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
DownloadAndAppendTcrClonotypes <- function(seuratObject, outPath = '.', dropExisting = T, overwriteTcrTable = F, allowMissing = FALSE){
  if (all(is.null(seuratObject[['BarcodePrefix']]))){
    stop('Seurat object lacks BarcodePrefix column')
  }

  i <- 0
  for (barcodePrefix in unique(unique(unlist(seuratObject[['BarcodePrefix']])))) {
    i <- i + 1
    print(paste0('Adding TCR clonotypes for prefix: ', barcodePrefix))

    vloupeId <- .FindMatchedVloupe(barcodePrefix)
    if (is.na(vloupeId)){
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
  tcr <- tcr[tcr$cdr3 != 'None',]

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
		ChainCDR3s = paste0(sort(unique(ChainCDR3s)), collapse = ","),
		CDR3s = paste0(sort(unique(cdr3)), collapse = ","),
		TRA = paste0(sort(unique(as.character(TRA))), collapse = ","),
		TRB = paste0(sort(unique(as.character(TRB))), collapse = ","),
		TRD = paste0(sort(unique(as.character(TRD))), collapse = ","),
		TRG = paste0(sort(unique(as.character(TRG))), collapse = ","),
		TRAV = paste0(sort(unique(as.character(TRAV))), collapse = ","),
		TRBV = paste0(sort(unique(as.character(TRBV))), collapse = ","),
		TRDV = paste0(sort(unique(as.character(TRDV))), collapse = ","),
		TRGV = paste0(sort(unique(as.character(TRGV))), collapse = ","),
    TRAJ = paste0(sort(unique(as.character(TRAJ))), collapse = ","),
    TRBJ = paste0(sort(unique(as.character(TRBJ))), collapse = ","),
    TRDJ = paste0(sort(unique(as.character(TRDJ))), collapse = ","),
    TRGJ = paste0(sort(unique(as.character(TRGJ))), collapse = ","),
    TRAC = paste0(sort(unique(as.character(TRAC))), collapse = ","),
    TRBC = paste0(sort(unique(as.character(TRBC))), collapse = ","),
    TRDC = paste0(sort(unique(as.character(TRDC))), collapse = ","),
    TRGC = paste0(sort(unique(as.character(TRGC))), collapse = ","),
    TRBD = paste0(sort(unique(as.character(TRBD))), collapse = ","),
    TRDD = paste0(sort(unique(as.character(TRDD))), collapse = ","),
		CloneNames = paste0(sort(unique(CloneName)), collapse = ",")  #this is imprecise b/c we count a hit if we match any chain, but this is probably what we often want
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


# #' @title CalculateTCRFreqForActivatedCells
# #' @description For the supplied cDNA Rows, this will query their gene expression and TCR readsets, identifying a) existing seurat objects, b) vloupe files. If both are found,
# #  it will download them, build a whitelist of highly activated cells (using SGS), calculate TCR frequencies, and return a dataframe
# #' @return A dataframe of results
# #' @param cDndIds A vector of cDNA rowIDs
# #' @param geneSetName The gene set name to use for SGS
# #' @param positivityThreshold The threshold to use for calling cells as positive
# #' @param outPrefix A string that will be prepended to all saved files
# #' @param invert If TRUE, those cells NOT positive for the gene set will be summarized, instead of positive cells
# #' @param doCleanup If TRUE, any downloaded files will be deleted on completion
# #' @param reduction The reduction (i.e. tsne or umap) that will be used when plotting
# #' @export
# #' @import Seurat
# #' @importFrom Biostrings readDNAStringSet
# #' @importFrom dplyr %>% group_by select n summarize
# CalculateTCRFreqForActivatedCells <- function(cDndIds, geneSetName = 'HighlyActivated', positivityThreshold = 0.5, outPrefix = './', invert = FALSE, doCleanup = FALSE, reduction = NULL) {
# 	print(paste0('Total cDNA records: ', length(cDndIds)))
# 	rows <- labkey.selectRows(
# 		baseUrl=.getBaseUrl(),
# 		folderPath=.getLabKeyDefaultFolder(),
# 		schemaName="singlecell",
# 		queryName="cdna_libraries",
# 		viewName="",
# 		colSort="-rowid",
# 		colFilter = makeFilter(c("rowid", "IN", paste0(cDndIds, collapse = ";"))),
# 		colSelect="rowid,readsetid,tcrreadsetid,hashingreadsetid",
# 		containerFilter=NULL,
# 		colNameOpt="rname"
# 	)
#
# 	if (nrow(rows) != length(cDndIds)) {
# 		print(paste0('Not all requested cDNAs found.  Row IDs found: ', paste0(unique(rows$rowid), collapse = ',')))
# 		return(NA)
# 	}
#
# 	gexReadsets <- unique(rows$readsetid)
# 	print(paste0('total GEX readsets: ', length(gexReadsets)))
#
# 	# Identify, download seuratObj, created from the appropriate readsetId:
# 	seuratRows <- labkey.selectRows(
# 		baseUrl=.getBaseUrl(),
# 		folderPath=.getLabKeyDefaultFolder(),
# 		schemaName="sequenceanalysis",
# 		queryName="outputfiles",
# 		colSort="-rowid",
# 		colSelect="rowid,readset",
# 		colFilter=makeFilter(
# 			c("readset", "IN", paste0(gexReadsets, collapse = ';')),
# 			c("category", "EQUAL", "Seurat Data")
# 		),
# 		containerFilter=NULL,
# 		colNameOpt="rname"
# 	)
#
# 	# Possible to have a duplicate:
# 	seuratRows <- unique(seuratRows)
#
# 	if (length(unique(seuratRows$readset)) != length(gexReadsets)) {
# 		missing <- gexReadsets[!(gexReadsets %in% unique(seuratRows$readset))]
# 		print(paste0('Not all requested cDNAs have seurat objects.  Readsets missing: ', paste0(unique(missing), collapse = ',')))
#
# 		return(NA)
# 	}
#
# 	downloadedFiles <- c()
# 	ret <- NA
# 	for (gexReadset in gexReadsets) {
# 		print(paste0('processing readset: ', gexReadset))
# 		row <- seuratRows[seuratRows$readset == gexReadset,,drop = F]
# 		if (nrow(row) > 1) {
# 			print('More than one seurat row found, using the most recent')
# 			row <- row[1,,drop = F]
# 		}
#
# 		if (is.na(row[['rowid']])) {
# 			warning(paste0('Error: RowID was NA, skipping'))
# 			print(row)
# 			next
# 		}
#
# 		f <- paste0(outPrefix, row[['rowid']], '.seurat.rds')
# 		DownloadOutputFile(row[['rowid']], f, overwrite = F)
# 		downloadedFiles <- c(downloadedFiles, f)
#
# 		# For each, apply metadata, TCR clones
# 		seuratObj <- readRDS(f)
#
# 		seuratObj <- DownloadAndAppendCellHashing(seuratObject = seuratObj)
# 		seuratObj <- QueryAndApplyCdnaMetadata(seuratObj)
#
#
# 		# TODO: refactor this from OOSAP
# 		seuratObj <- OOSAP::ClassifySGSAndApply(seuratObj = seuratObj, geneSetName = 'Positive', geneList = OOSAP::Phenotyping_GeneList()[[geneSetName]], positivityThreshold = positivityThreshold, reduction = reduction)
# 		if (invert) {
# 			print('Selecting cells without the provided signature')
# 			barcodeWhitelist <- colnames(seuratObj)[!seuratObj$Positive.Call & !is.na(seuratObj$cDNA_ID)]
# 		} else {
# 			barcodeWhitelist <- colnames(seuratObj)[seuratObj$Positive.Call & !is.na(seuratObj$cDNA_ID)]
# 		}
#
# 		i <- 0
# 		for (barcodePrefix in unique(unlist(seuratObj[['BarcodePrefix']]))) {
# 			i <- i + 1
#
# 			vloupeId <- .FindMatchedVloupe(barcodePrefix)
# 			if (is.na(vloupeId)){
# 				stop(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
# 			}
#
# 			#TCR libraryId
# 			tcrLibRows <- labkey.selectRows(
# 				baseUrl=.getBaseUrl(),
# 				folderPath=.getLabKeyDefaultFolder(),
# 				schemaName="sequenceanalysis",
# 				queryName="outputfiles",
# 				colSort="-rowid",
# 				colSelect="library_id,analysis_id,readset",
# 				colFilter=makeFilter(c("rowid", "EQUALS", vloupeId)),
# 				containerFilter=NULL,
# 				colNameOpt="rname"
# 			)
# 			libraryId <- tcrLibRows$library_id[1]
# 			analysisId <- tcrLibRows$analysis_id[1]
# 			tcrReadset <- tcrLibRows$readset[1]
# 			print(paste0('TCR library ID: ', libraryId))
# 			print(paste0('TCR readset: ', tcrReadset))
#
# 			#All clonotypes
# 			clonotypeFile <- file.path(outPrefix, paste0(barcodePrefix, '_all_contig_annotations.csv'))
# 			.DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = clonotypeFile, overwrite = T)
# 			if (!file.exists(clonotypeFile)){
# 				stop(paste0('Unable to download clonotype file for prefix: ', barcodePrefix))
# 			}
# 			downloadedFiles <- c(downloadedFiles, clonotypeFile)
#
# 			#FASTA:
# 			fastaFile <- file.path(outPrefix, paste0(barcodePrefix, '_consensus.fasta'))
# 			.DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = fastaFile, overwrite = T, fileName = 'consensus.fasta')
# 			if (!file.exists(fastaFile)){
# 				stop(paste0('Unable to download clonotype FASTA for prefix: ', barcodePrefix))
# 			}
# 			downloadedFiles <- c(downloadedFiles, fastaFile)
#
# 			tcrData <- .ProcessTcrClonotypes(clonotypeFile)
# 			if (!is.null(barcodePrefix)){
# 				tcrData$barcode <- as.character(tcrData$barcode)
# 				tcrData$barcode <- paste0(barcodePrefix, '_', tcrData$barcode)
# 				tcrData$barcode <- as.factor(tcrData$barcode)
# 			}
#
# 			retain <- intersect(barcodeWhitelist, tcrData$barcode)
#
# 			print(paste0('initial barcodes with TCR call: ', nrow(tcrData)))
# 			pct1 <- round(length(retain) / length(unique(tcrData$barcode)), 2)
# 			pct2 <- round(length(retain) / length(barcodeWhitelist), 2)
#
# 			print(paste0('overlapping with activated cells: ', length(retain), ' (', pct1, ' of TCR calls, ',pct2,' of positive cells)'))
# 			tcrData <- tcrData[tcrData$barcode %in% retain,]
# 			if (nrow(tcrData) == 0) {
# 				print(paste0('no rows for prefix: ', barcodePrefix))
# 				next
# 			}
#
# 			#summarize metadata
# 			meta <- data.frame(
# 				barcode = colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix],
# 				SubjectId = as.character(seuratObj$SubjectId[seuratObj$BarcodePrefix == barcodePrefix]),
# 				Stim = as.character(seuratObj$Stim[seuratObj$BarcodePrefix == barcodePrefix]),
# 				population = as.character(seuratObj$Population[seuratObj$BarcodePrefix == barcodePrefix]),
# 				date = as.character(seuratObj$SampleDate[seuratObj$BarcodePrefix == barcodePrefix]),
# 				cdna = as.character(seuratObj$cDNA_ID[seuratObj$BarcodePrefix == barcodePrefix]),
# 				libraryId = c(libraryId),
# 				analysisId = c(analysisId)
# 			)
# 			meta$SampleName <- paste0(meta$SubjectId, '_', meta$Stim)
# 			tcrData <- merge(tcrData, meta, by = c('barcode'), all.x = T)
#
# 			if (sum(is.na(tcrData$SubjectId)) > 0) {
# 				f <- paste0(outPrefix, 'temp.txt')
# 				write.table(tcrData, file = f, sep = '\t', quote = F, row.names = F)
# 				stop(paste0('Missing subject IDs!  See ', f, ' for table of results'))
# 			}
#
# 			# Group
# 			tcrData <- tcrData %>% group_by(SampleName, SubjectId, population, date, cdna, libraryId, analysisId, chain, cdr3, v_gene, d_gene, j_gene, c_gene, raw_clonotype_id, raw_consensus_id) %>% summarize(count = dplyr::n())
# 			names(tcrData) <- c('SampleName', 'SubjectId', 'population', 'date', 'cdna', 'libraryId', 'analysisId', 'locus', 'cdr3', 'vHit', 'dHit', 'jHit', 'cHit', 'cloneId', 'consensus_id', 'count')
#
# 			tcrData <- tcrData %>% group_by(cdna) %>% mutate(totalCells = dplyr::n())
# 			tcrData$fraction = tcrData$count / tcrData$totalCells
#
# 			# Merge sequence:
# 			fastaData <- Biostrings::readDNAStringSet(fastaFile)
# 			seqDf <- data.frame(consensus_id = names(fastaData), sequence = paste(fastaData))
#
# 			tcrData <- merge(tcrData, seqDf, by = c('consensus_id'), all.x = T)
# 			tcrData <- tcrData[!(names(tcrData) %in% c('consensus_id', 'totalCells', 'Stim'))]
# 			tcrData[tcrData == 'None'] <- NA
#
# 			if (all(is.na(ret))) {
# 				ret <- tcrData
# 			} else {
# 				ret <- rbind(ret, tcrData)
# 			}
# 		}
# 	}
#
# 	if (doCleanup) {
# 		print('Cleaning up downloaded files')
# 		for (f in downloadedFiles) {
# 			unlink(f)
# 		}
# 	}
#
# 	return(ret)
# }

#' @title CalculateTCRFreqForActivatedCellsAndImport
#' @description For the supplied cDNA Rows, this will call CalculateTCRFreqForActivatedCells(), and also import the results into the provided assay
#' @return Returns the object representation of the experiment batch (from Rlabkey::labkey.experiment.saveBatch)
#' @param cDndIds A vector of cDNA rowIDs
#' @param geneSetName The gene set name to use for SGS
#' @param positivityThreshold The threshold to use for calling cells as positive
#' @param outPrefix A string that will be prepended to all saved files
#' @param workbook The target workbook to save results
#' @param assayName The name of the target assay
#' @param populationNameSuffix This string will be added to the end of the population field on the saved assay rows, appended to the associated cDNA population
#' @param invert If TRUE, those cells NOT positive for the gene set will be summarized, instead of positive cells
#' @param doCleanup If TRUE, any downloaded files will be deleted on completion
#' @param minCells If provided, data will only be imported if at least this many cells exist for the cDNA library
#' @export
#' @import Seurat
#' @importFrom dplyr %>% group_by select n summarize
CalculateTCRFreqForActivatedCellsAndImport <- function(cDndIds, workbook = NULL, geneSetName = 'HighlyActivated', positivityThreshold = 0.5, outPrefix = './', assayName = 'TCRdb', populationNameSuffix = '-HA', invert = FALSE, doCleanup = FALSE, minCells = NULL) {
	folder <- .getLabKeyDefaultFolder()
	if (!is.null(workbook)) {
		folder <- paste0(folder, workbook)
	}

	resultDataFrame <- CalculateTCRFreqForActivatedCells(cDndIds, geneSetName = geneSetName, positivityThreshold = positivityThreshold, outPrefix = outPrefix, invert = invert)
	if (all(is.na(resultDataFrame))) {
		print('no results, skipping')
		return(NULL)
	}

	resultDataFrame$calculatedPopulation <- paste0(as.character(resultDataFrame$population), populationNameSuffix)

	analysisIds <- unique(resultDataFrame$analysisId)
	print(paste0('Total analyses: ', length(analysisIds)))

	runList <- list()
	for (analysisId in analysisIds) {
		df <- resultDataFrame[resultDataFrame$analysisId == analysisId,]
		print(paste0('Preparing run for analysis Id: ',analysisId,', total rows: ', nrow(df), ', total cells: ', sum(df$count)))
		if (!is.null(minCells)) {
			totals <- df %>% group_by(cdna) %>% summarize(total = sum(count))
			toKeep <- unique(totals$cdna[totals$total >= minCells])
			if (length(toKeep) != length(unique(df$cdna))) {
				print('The following cDNA libraries will be skipped due to low cells:')
				print(paste0(totals$cdna[totals$total < minCells], ': ', totals$total[totals$total < minCells]))

				df <- df[df$cdna %in% toKeep,]
			}
		}

		if (sum(is.na(df$SubjectId)) > 0) {
			f <- paste0(outPrefix, 'temp.txt')
			write.table(df, file = f, sep = '\t', quote = F, row.names = F)
			stop(paste0('Missing subject IDs!  See ', f, ' for table of results'))
		}

		if (nrow(df) > 0) {
			print('run passed, adding')
			run <- labkey.experiment.createRun(list(name = paste0('AnalysisId: ', analysisId, populationNameSuffix), properties = list(assayName = '10x', analysisId = analysisId)), dataRows = df)
			runList <- append(runList, list(run))
		} else {
			print('no rows passed, skipping')
		}
	}

	print(paste0('total runs: ', length(runList)))
	if (length(runList) > 0) {
		labkey.experiment.saveBatch(
			baseUrl=.getBaseUrl(),
			folderPath=folder,
			assayConfig=list(assayName=assayName, providerName='TCRdb'),
			runList = runList
		)
	} else {
		print('no runs to save')
	}
}