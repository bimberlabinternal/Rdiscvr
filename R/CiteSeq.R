#' @include LabKeySettings.R
#' @include Utils.R

#' @title DownloadAndAppendCiteSeq
#'
#' @description Downloads matching Cite-seq counts using barcodePrefix on the seurat object
#' @param seuratObj, A Seurat object.
#' @param featureLabelTable An optional TSV file listing the markers and metadata.  It must contain a header row with at least the columns (lowercase): tagname, sequence. If it contains the column markername, this will be used to replace the rownames of the matrix.
#' @param minRowSum  If provided, any ADTs with a rowSum less than this value will be dropped.
#' @param adtWhitelist An optional character vector of tag names, which must exactly match the rownames of the ADT count matrix.  If provided, data will be limited to these tags.
#' @param assayName The name of the assay to store the ADT data.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendCiteSeq <- function(seuratObj, outPath = '.', assayName = 'ADT', featureLabelTable = NULL, minRowSum = 5, adtWhitelist = NULL){
	if (is.null(seuratObj[['BarcodePrefix']])){
		stop('Seurat object lacks BarcodePrefix column')
	}

	for (barcodePrefix in unique(unique(unlist(seuratObj[['BarcodePrefix']])))) {
		print(paste0('Possibly adding cell hashing data for prefix: ', barcodePrefix))

		citeseqId <- .FindMatchedCiteSeq(barcodePrefix)
		if (is.null(citeseqId)){
			print(paste0('CITE-seq not used for prefix: ', barcodePrefix, ', skipping'))
			next
		} else if (is.na(citeseqId)){
			stop(paste0('Unable to find CITE-seq matrix file for prefix: ', barcodePrefix))
		}

		countDir <- file.path(outPath, paste0(barcodePrefix, '_citeseqCounts'))
		if (!dir.exists(countDir)) {
			dir.create(countDir)
		}

		countDir <- .DownloadCiteSeqDir(outputFileId = citeseqId, localBaseDir = countDir)
		if (!dir.exists(countDir)){
			stop(paste0('Unable to download calls table for prefix: ', barcodePrefix, ', expected file: ', countDir))
		}

		seuratObj <- AppendCiteSeq(seuratObj = seuratObj, countMatrixDir = countDir, barcodePrefix = barcodePrefix, assayName = assayName, featureLabelTable = featureLabelTable, adtWhitelist = adtWhitelist, minRowSum = minRowSum, skipNormalize = T)

		seuratObj <- NormalizeData(seuratObj, assay = assayName, normalization.method = "CLR")
		seuratObj <- ScaleData(seuratObj, assay = assayName)
	}

	cellhashR::PlotCiteSeqCountData(seuratObj, assayName = assayName)

	return(seuratObj)
}

.DownloadCiteSeqDir <- function(outputFileId, localBaseDir = './', overwriteFiles = T, mergeFolders = F) {
	rows <- labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		viewName="",
		colSort="-rowid",
		colSelect="rowid,workbook/workbookid,dataid/webdavurlrelative",
		colFilter=makeFilter(c("rowid", "EQUAL", outputFileId)),
		containerFilter=NULL,
		colNameOpt="rname"
	)

	if (nrow(rows) != 1) {
		stop(paste0('More than one matching file found, this should not occur.  RowId: ', outputFileId))
	}

	wb <- rows[['workbook_workbookid']]
	if (is.na(wb) || is.null(wb)){
		wb <- ''
	}

	#The database tracks the matrix, but we want the whole folder
	remotePath <- rows[['dataid_webdavurlrelative']]
	remotePath <- gsub(x = remotePath, pattern = 'matrix.mtx.gz', replacement = '')

	success <- labkey.webdav.downloadFolder(
		baseUrl=lkBaseUrl,
		folderPath=paste0(.getLabKeyDefaultFolder(),wb),
		remoteFilePath = remotePath,
		overwriteFiles = overwriteFiles,
		mergeFolders = mergeFolders,
		localBaseDir = localBaseDir
	)

	if (!success || !dir.exists(localBaseDir) || !file.exists(paste0(localBaseDir, '/umi_count/matrix.mtx.gz'))) {
		stop(paste0('labkey.webdav.downloadFolder failed, expected: ', localBaseDir))
	}

	return(paste0(localBaseDir, '/umi_count'))
}

#' @import Rlabkey
.FindMatchedCiteSeq <- function(loupeDataId){
	#Note: the seurat object gets associated with the GEX readset, so look based on this:
	rows <- labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		viewName="",
		colSort="-rowid",
		colSelect="readset,library_id",
		colFilter=makeFilter(c("rowid", "EQUAL", loupeDataId)),
		containerFilter=NULL,
		colNameOpt="rname"
	)

	if (nrow(rows) == 0) {
		print(paste0("Loupe File ID: ", loupeDataId, " not found"))
		return(NA)
	}

	readset <- unique(rows[['readset']])

	if (is.na(readset) || is.null(readset)) {
		print("readset is NA/NULL")
		return(NA)
	}

	libraryId <- unique(rows[['library_id']])

	#determine whether we expect cell hashing to be used:
	cDNAs <- labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="tcrdb",
		queryName="cdnas",
		viewName="",
		colSort="-rowid",
		colFilter = makeFilter(c("readsetId", "EQUALS", readset)),
		colSelect="rowid,readsetid,citeseqreadsetid",
		containerFilter=NULL,
		colNameOpt="rname"
	)

	if (nrow(cDNAs) == 0) {
		stop(paste0('No cDNA records found for GEX readset: ', readset))
	} else if (sum(!is.na(cDNAs$citeseqreadsetid)) == 0) {
		print(paste0('The cDNA library does not use citeseq, aborting'))
		return(NULL)
	}

	rows <- suppressWarnings(labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		colSort="-rowid",
		colSelect="rowid",
		colFilter=makeFilter(c("readset", "EQUAL", readset),
		c("category", "EQUAL", "Seurat CITE-Seq Count Matrix"),
		c("library_id", "EQUAL", libraryId)),
		containerFilter=NULL,
		colNameOpt="rname"
	))

	ret <- NULL
	if (nrow(rows) == 0){
		print(paste0("Output of type 'Seurat CITE-Seq Count Matrix' not found.  Readset: ", readset, ", libraryId: ", libraryId))
	} else {
		ret <- rows[1]
	}

	if (all(is.null(ret))){
		print("Trying to find output of type: 'CITE-Seq Count Matrix' using GEX readset")

		rows <- suppressWarnings(labkey.selectRows(
			baseUrl=lkBaseUrl,
			folderPath=.getLabKeyDefaultFolder(),
			schemaName="sequenceanalysis",
			queryName="outputfiles",
			colSort="-rowid",
			colSelect="rowid,",
			colFilter=makeFilter(c("readset", "EQUAL", readset),
			c("category", "EQUAL", "CITE-Seq Count Matrix"),
			c("library_ld", "EQUAL", libraryId)
			),
			containerFilter=NULL,
			colNameOpt="rname"
		))

		if (nrow(rows) == 0){
			print("Not found")
		} else {
			ret <- rows[1]
			print("Found output")
		}
	}

	if (all(is.null(ret))){
		print("Trying to find output of type: 'CITE-Seq Count Matrix', using cite-seq readset")
		citeseqReadset <- unique(cDNAs$citeseqreadsetid)

		rows <- suppressWarnings(labkey.selectRows(
			baseUrl=lkBaseUrl,
			folderPath=.getLabKeyDefaultFolder(),
			schemaName="sequenceanalysis",
			queryName="outputfiles",
			colSort="-rowid",
			colSelect="rowid,",
			colFilter=makeFilter(c("readset", "EQUAL", citeseqReadset),
			c("category", "EQUAL", "CITE-Seq Count Matrix"),
			c("library_ld", "EQUAL", libraryId)
			),
			containerFilter=NULL,
			colNameOpt="rname"
		))

		if (nrow(rows) == 0){
			print("Not found")
		} else {
			ret <- rows[1]
			print("Found output")
		}
	}

	if (all(is.null(ret))) {
		return(NA)
	} else {
		if (length(ret) > 0){
			print('More than one matching file found, using most recent')
		}

		return(ret$rowid[1])
	}
}