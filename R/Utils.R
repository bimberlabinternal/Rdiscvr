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