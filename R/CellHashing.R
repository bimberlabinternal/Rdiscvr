#' @include LabKeySettings.R
#' @include Utils.R


#' @title DownloadAndAppendCellHashing
#'
#' @description Downloads matching Cell Hashing data using barcodePrefix on the seurat object and appends it to metadata
#' @param seuratObject, A Seurat object.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendCellHashing <- function(seuratObject, outPath = '.'){
	if (is.null(seuratObject[['BarcodePrefix']])){
		stop('Seurat object lacks BarcodePrefix column')
	}

	for (barcodePrefix in unique(unique(unlist(seuratObject[['BarcodePrefix']])))) {
		print(paste0('Possibly adding cell hashing data for prefix: ', barcodePrefix))

		cellHashingId <- .FindMatchedCellHashing(barcodePrefix)
		if (is.null(cellHashingId)){
			print(paste0('Cell hashing not used for prefix: ', barcodePrefix, ', skipping'))
			next
		} else if (is.na(cellHashingId)){
			stop(paste0('Unable to find cellHashing calls table file for prefix: ', barcodePrefix))
		}

		callsFile <- file.path(outPath, paste0(barcodePrefix, '_cellHashingCalls.csv'))
		DownloadOutputFile(outputFileId = cellHashingId, outFile = callsFile, overwrite = T)
		if (!file.exists(callsFile)){
			stop(paste0('Unable to download calls table for prefix: ', barcodePrefix))
		}

		seuratObject <- cellhashR::AppendCellHashing(seuratObj = seuratObject, barcodeCallFile = callsFile, barcodePrefix = barcodePrefix)
	}

	return(seuratObject)
}

#' @import Rlabkey
.FindMatchedCellHashing <- function(loupeDataId){
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
		colSelect="rowid,readsetid,hashingreadsetid,sortid/hto",
		containerFilter=NULL,
		colNameOpt="rname"
	)

	if (nrow(cDNAs) == 0) {
		stop(paste0('No cDNA records found for GEX readset: ', readset))
	} else if (sum(!is.na(cDNAs$sortid_hto)) == 0 || sum(!is.na(cDNAs$hashingreadsetid)) == 0) {
		print(paste0('The cDNA library does not use cell hashing, aborting'))
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
		c("category", "EQUAL", "Seurat Cell Hashing Calls"),
		c("library_id", "EQUAL", libraryId)),
		containerFilter=NULL,
		colNameOpt="rname"
	))

	ret <- NULL
	if (nrow(rows) == 0){
		print(paste0("Output of type 'Seurat Cell Hashing Calls' not found.  Readset: ", readset, ", libraryId: ", libraryId))
	} else {
		ret <- rows[1]
	}

	if (all(is.null(ret))){
		print("Trying to find output of type: '10x GEX Cell Hashing Calls'")
		rows <- suppressWarnings(labkey.selectRows(
			baseUrl=lkBaseUrl,
			folderPath=.getLabKeyDefaultFolder(),
			schemaName="sequenceanalysis",
			queryName="outputfiles",
			colSort="-rowid",
			colSelect="rowid,",
			colFilter=makeFilter(c("readset", "EQUAL", readset),
			c("category", "EQUAL", "10x GEX Cell Hashing Calls"),
			c("library_ld", "EQUAL", libraryId)),
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