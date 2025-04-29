#' @import utils

#' @title DownloadOutputFile
#' @description Downloads an output file tracked in LabKey to the local filesystem.
#' @param outputFileId The rowid of the sequence outputfiles on the webserver
#' @param outFile The local path to write this file
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param pathTranslator There are instances where the goal is to download some file located relative to the actual outputfile. For example, 10x data has an outputfile for the Loupe file, but it might be useful to download the raw counts. To accomplish this, a function can be provided that accepts the primary file remote path, and returns a modified path to download.
#' @param showProgressBar If true, a progress bar will be shown
#' @export
#'
#' @import Rlabkey
DownloadOutputFile <- function(outputFileId, outFile, overwrite = T, pathTranslator = NULL, showProgressBar = FALSE) {
	return(.DownloadOutputFileOrDirectory(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = pathTranslator, asDirectory = FALSE, showProgressBar = showProgressBar))
}

#' @title DownloadOutputFile
#' @description Downloads an output file tracked in LabKey to the local filesystem.
#' @param outputFileId The rowid of the sequence outputfiles on the webserver
#' @param outFile The local path to write this file
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param pathTranslator There are instances where the goal is to download some file located relative to the actual outputfile. For example, 10x data has an outputfile for the Loupe file, but it might be useful to download the raw counts. To accomplish this, a function can be provided that accepts the primary file remote path, and returns a modified path to download.
#' @param showProgressBar If true, a progress bar will be shown
#' @export
#'
#' @import Rlabkey
DownloadOutputDirectoryFromOutputFile <- function(outputFileId, outFile, overwrite = T, pathTranslator = NULL, showProgressBar = FALSE) {
	if (is.null(pathTranslator)) {
		stop('When attempting to download a folder relative to an outputfile, you must provide a pathTranslator function')
	}
	return(.DownloadOutputFileOrDirectory(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = pathTranslator, asDirectory = TRUE, showProgressBar = showProgressBar))
}

.DownloadOutputFileOrDirectory <- function(outputFileId, outFile, asDirectory, overwrite = T, pathTranslator = NULL, showProgressBar = FALSE) {
	if (is.na(outputFileId)) {
		stop('Output file ID cannot be NA')
	}

	# NOTE: labkey.webdav.downloadFolder expects the base folder to exist, so only perform this check for files.
	# Otherwise let labkey.webdav.downloadFolder handle overwrite behaviors
	if (!asDirectory) {
		if (file.exists(outFile)) {
			if (!overwrite) {
				print(paste0("File exists, will not overwrite: ", outFile))
				return(outFile)
			} else {
				unlink(outFile)
			}
		}
	}

	rows <- labkey.selectRows(
		baseUrl=.getBaseUrl(),
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

	if (nrow(rows) > 1) {
		stop(paste0('More than one matching file found, this should not occur.  RowId: ', outputFileId))
	} else if (nrow(rows) == 0) {
			stop(paste0('File not found. RowId: ', outputFileId))
	}

	wb <- rows[['workbook_workbookid']]
	if (is.na(wb) || is.null(wb)){
		wb <- ''
	}

	remotePath <- rows[['dataid_webdavurlrelative']]
	if (!is.null(pathTranslator)) {
		remotePath <- pathTranslator(remotePath)
	}

	if (asDirectory) {
		success <- labkey.webdav.downloadFolder(
			baseUrl=.getBaseUrl(),
			folderPath=paste0(.getLabKeyDefaultFolder(),wb),
			remoteFilePath = remotePath,
			overwriteFiles = overwrite,
			localBaseDir = outFile,
			showProgressBar = showProgressBar
		)

		if (!success) {
			stop(paste0('labkey.webdav.downloadFolder failed for file: ', remotePath))
		}
	} else {
		outFileTmp <- paste0(outFile, ".tmp")
		if (file.exists(outFileTmp)) {
			unlink(outFileTmp)
		}

		success <- labkey.webdav.get(
			baseUrl=.getBaseUrl(),
			folderPath=paste0(.getLabKeyDefaultFolder(),wb),
			remoteFilePath = remotePath,
			overwrite = overwrite,
			localFilePath = outFileTmp,
			showProgressBar = showProgressBar
		)

		if (success) {
			file.rename(outFileTmp, outFile)
		}

		if (!success || !file.exists(outFile)) {
			stop(paste0('labkey.webdav.get failed for file: ', remotePath))
		}
	}

	return(outFile)
}

#' @title UploadOutputFile
#' @description Uploads a local file and saves as an output file tracked in LabKey.
#' @param localPath The local filepath to the file
#' @param workbook The numeric ID of the target workbook
#' @param name The name for this output
#' @param category The category for this output
#' @param description A description for this output
#' @param genomeId The integer genome ID associated with this output
#' @param readsetId The integer readsetId associated with this output
#' @param remoteFilename If null, the local name will be used

#' @export
#'
#' @import Rlabkey
UploadOutputFile <- function(localPath, workbook, name, category, description, readsetId = NULL, genomeId, remoteFilename = NULL) {
	if (!file.exists(localPath)) {
		stop(paste0('Unable to find file: ', localPath))
	}

	if (is.null(remoteFilename)) {
		remoteFilename <- basename(localPath)
	}

	remotePath <- paste0('sequenceOutputs/', remoteFilename)
	folderExists <- labkey.webdav.pathExists(
		baseUrl=.getBaseUrl(),
		folderPath=paste0(.getLabKeyDefaultFolder(), workbook),
		remoteFilePath = 'sequenceOutputs'
	)

	if (!folderExists) {
		success <- labkey.webdav.mkDirs(
			baseUrl=.getBaseUrl(),
			folderPath=paste0(.getLabKeyDefaultFolder(), workbook),
			remoteFilePath = 'sequenceOutputs'
		)

		if (!success) {
			stop('Unable to create folder: sequenceOutputs')
		}
	}

	remoteExists <- labkey.webdav.pathExists(
		baseUrl=.getBaseUrl(),
		folderPath=paste0(.getLabKeyDefaultFolder(), workbook),
		remoteFilePath = remotePath
	)

	if (remoteExists) {
		stop(paste0('Remote path already exists: ', remotePath))
	}

	print('Uploading file')
	success <- labkey.webdav.put(
		baseUrl=.getBaseUrl(),
		folderPath=paste0(.getLabKeyDefaultFolder(), workbook),
		remoteFilePath = remotePath,
		localFile = localPath
	)
	print('Done')

	if (!success) {
		stop('Failed to upload file')
	}

	#TODO: need to create ExpData object from this file
	expDataId <- NA

	toInsert <- data.frame(
		fileid = expDataId,
		category = category,
		description = description,
		name = name,
		readset = readsetId,
		library_id = genomeId
	)

	inserted <- labkey.insertRows(
		baseUrl=.getBaseUrl(),
		folderPath=paste0(.getLabKeyDefaultFolder(), workbook),
		schemaName = 'sequenceanalysis',
		queryName = 'outputfiles',
		toInsert = toInsert
	)

	return(inserted$rows[[1]]$rowid)
}

#' @title UploadOutputFile
#' @description Uploads a local file and saves as an output file tracked in LabKey.
#' @param localPath The local filepath to the file
#' @param workbook The numeric ID of the target workbook
#' @param name The name for this output
#' @param category The category for this output
#' @param description A description for this output
#' @param readsetId The integer readsetId associated with this output
#' @param genomeId The integer genome ID associated with this output
#' @param remoteFilename If null, the local name will be used

#' @export
#'
#' @import Rlabkey
UploadSeuratObject <- function(localPath, workbook, name, category, description, readsetId = NULL, genomeId, remoteFilename = NULL) {
	return(UploadOutputFile(localPath, workbook = workbook, name = name, category = 'Seurat Object', description = description, readsetId = readsetId, genomeId = genomeId, remoteFilename = remoteFilename))
}

#' @title DownloadMetadataForSeuratObject
#' @description Downloads an output file tracked in LabKey to the local filesystem.
#' @param outputFileId The rowid of the sequence outputfiles on the webserver
#' @param outFile The local path to write this file
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param returnDataFrame If true, the metadata will be read to a data.frame and returned
#' @param deleteFile If true, the outFile will be deleted. This only makes sense when used with returnDataFrame=TRUE
#' @param showProgressBar If true, a progress bar will be shown
#' @export
DownloadMetadataForSeuratObject <- function(outputFileId, outFile, overwrite = TRUE, returnDataFrame = FALSE, deleteFile = returnDataFrame, showProgressBar = FALSE) {
	DownloadOutputFile(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = function(x){
		return(gsub(x, pattern = 'seurat.rds', replacement = 'seurat.meta.txt.gz'))
	}, showProgressBar = showProgressBar)

	if (returnDataFrame) {
		df <- read.table(outFile, header = T, sep = ',')
		if (deleteFile) {
			unlink(outFile)
		}

		return(df)
	}
}

.ResolveLoupeIdFromDeleted <- function(loupeId, throwOnError = FALSE) {
	print(paste0('The following loupeId was not found, looking for deleted records: ', loupeId))
	translatedRows <- suppressWarnings(labkey.selectRows(
		baseUrl=.getBaseUrl(),
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="singlecell",
		queryName="singlecellDatasets",
		colFilter = makeFilter(c("loupeFileId", "EQUALS", loupeId)),
		colSelect="loupeFileId,readsetId",
		containerFilter=NULL,
		colNameOpt="rname"
	))

	if (nrow(translatedRows) > 0) {
		print(paste0('The following ID was matched to a deleted object: ', paste0(translatedRows$loupefileid, collapse = ';')))
		return(translatedRows)
	}

	if (throwOnError) {
		stop('The following loupeId was not found, even after looking for deleted records: ', loupeId)
	}

	return(NULL)
}


#' @title GenerateSRATable
#' @description Returns a table with SRA-related information for a given set of cDNA_IDs
#' @param cDNA_IDs A vector of integer cDNA_IDs
#' @export
#'
#' @import Rlabkey
GenerateSRATable <- function(cDNA_IDs) {
	if (any(is.na(cDNA_IDs))) {
		stop('Cannot provide NA cDNA IDs')
	}

	# Avoid factors:
	cDNA_IDs <- as.integer(as.character(cDNA_IDs))

	if (any(is.na(cDNA_IDs))) {
		stop('Non-numeric cDNA IDs provided')
	}

	dat <- labkey.selectRows(
		baseUrl=.getBaseUrl(),
		folderPath=.getLabKeyDefaultFolder(),
		schemaName="singlecell",
		queryName="cdna_libraries",
		colSelect="rowid,sortId/sampleId/subjectId,sortId/sampleId/sampledate,sortId/sampleId/stim,sortId/sampleId/assayType,sortId/sampleId/tissue,sortId/population,sortId/hto,sortId/hto/adaptersequence,readsetId/sraRuns,tcrReadsetId/sraRuns,hashingReadsetId/sraRuns,citeseqReadsetId/sraRuns",
		colFilter=makeFilter(c("rowid", "IN", paste0(cDNA_IDs, collapse = ';'))),
		colNameOpt="rname"
	)
	names(dat) <- c('cDNA_ID', 'SubjectId', 'SampleDate', 'Stim', 'AssayType', 'Tissue', 'Population', 'HTO', 'HTO Sequence', 'GEX SRA', 'VDJ SRA', 'Hashing SRA', 'CITE-seq SRA')

	return(dat)
}

.GetAssayMetaSlotName <- function(assayObj) {
	slotName <- ifelse('meta.features' %in% methods::slotNames(assayObj), yes = 'meta.features', no = 'meta.data')
	if (! slotName %in% methods::slotNames(assayObj)) {
		stop(paste0('Assay object lacks slot: ', slotName))
	}

	return(slotName)
}

.GetAssayMeta <- function(assayObj) {
	return(methods::slot(assayObj, .GetAssayMetaSlotName(assayObj)))
}
