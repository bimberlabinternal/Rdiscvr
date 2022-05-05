#' @import utils

#' @title DownloadOutputFile
#' @description Downloads an output file tracked in LabKey to the local filesystem.
#' @param outputFileId The rowid of the sequence outputfiles on the webserver
#' @param outFile The local path to write this file
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param pathTranslator There are instances where the goal is to download some file located relative to the actual outputfile. For example, 10x data has an outputfile for the Loupe file, but it might be useful to download the raw counts. To accomplish this, a function can be provided that accepts the primary file remote path, and returns a modified path to download.
#' @export
#'
#' @import Rlabkey
DownloadOutputFile <- function(outputFileId, outFile, overwrite = T, pathTranslator = NULL) {
	return(.DownloadOutputFileOrDirectory(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = pathTranslator, asDirectory = F))
}

#' @title DownloadOutputFile
#' @description Downloads an output file tracked in LabKey to the local filesystem.
#' @param outputFileId The rowid of the sequence outputfiles on the webserver
#' @param outFile The local path to write this file
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param pathTranslator There are instances where the goal is to download some file located relative to the actual outputfile. For example, 10x data has an outputfile for the Loupe file, but it might be useful to download the raw counts. To accomplish this, a function can be provided that accepts the primary file remote path, and returns a modified path to download.
#' @export
#'
#' @import Rlabkey
DownloadOutputDirectoryFromOutputFile <- function(outputFileId, outFile, overwrite = T, pathTranslator = NULL) {
	if (is.null(pathTranslator)) {
		stop('When attempting to download a folder relative to an outputfile, you must provide a pathTranslator function')
	}
	return(.DownloadOutputFileOrDirectory(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = pathTranslator, asDirectory = T))
}

.DownloadOutputFileOrDirectory <- function(outputFileId, outFile, asDirectory, overwrite = T, pathTranslator = NULL) {
	if (is.na(outputFileId)) {
		stop('Output file ID cannot be NA')
	}

	# NOTE: labkey.webdav.downloadFolder expects the base folder to exist, so only perform this check for files.
	# Otherwise let labkey.webdav.downloadFolder handle overwrite behaviors
	if (!asDirectory && file.exists(outFile) & !overwrite) {
		print(paste0("File exists, will not overwrite: ", outFile))
		return(outFile)
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

	if (nrow(rows) != 1) {
		stop(paste0('More than one matching file found, this should not occur.  RowId: ', outputFileId))
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
			localBaseDir = outFile
		)

		if (!success) {
			stop(paste0('labkey.webdav.downloadFolder failed for file: ', remotePath))
		}
	} else {
		success <- labkey.webdav.get(
			baseUrl=.getBaseUrl(),
			folderPath=paste0(.getLabKeyDefaultFolder(),wb),
			remoteFilePath = remotePath,
			overwrite = overwrite,
			localFilePath = outFile
		)

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