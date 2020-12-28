
.lkdefaults <- new.env(parent=emptyenv());

#' @title SetLabKeyDefaults
#' @description Can be used to set the default LabKey baseUrl and folderPath for all queries of this module.
#'
#' @param baseUrl The baseUrl of the LabKey server.
#' @param defaultFolder The default folder path (which should generally be the parent folder and not a workbook) for queries.
#' @export
SetLabKeyDefaults <- function(baseUrl="", defaultFolder = ""){
	if (baseUrl != "") {
		.lkdefaults$baseUrl = .normalizeSlash(baseUrl, trailing = TRUE, leading = FALSE);
	}

	if (defaultFolder != "") {
		.lkdefaults$defaultFolder = .normalizeSlash(defaultFolder, trailing = TRUE);
	}
}

.getLabKeyDefaultFolder <- function(defaultFolder = NULL) {
	folder <- .lkdefaults$defaultFolder
	if (is.null(folder)) {
		stop(paste("defaultFolder is null or has not been set yet. See SetLabKeyDefaults"))
	}

	return (folder)
}

.getBaseUrl <- function(baseUrl = NULL) {
	url <- .lkdefaults$baseUrl
	if (is.null(url)) {
		stop(paste("baseUrl is null or has not been set yet. See SetLabKeyDefaults"))
	}

	return (url)
}

.normalizeSlash <- function(folderPath, leading = T, trailing = T) {
	## Formatting
	folderPath <- gsub("[\\]", "/", folderPath)

	if (trailing) {
		if(substr(folderPath, nchar(folderPath), nchar(folderPath))!="/")
			folderPath <- paste(folderPath,"/",sep="")
	} else {
		if(substr(folderPath, nchar(folderPath), nchar(folderPath))=="/")
			folderPath <- substr(folderPath,1, nchar(folderPath)-1)
	}

	if (leading) {
		if(substr(folderPath, 1, 1)!="/")
			folderPath <- paste("/",folderPath,sep="")
	} else {
		if(substr(folderPath, 1, 1)=="/")
			folderPath <- substr(folderPath,2, nchar(folderPath))
	}

	return(folderPath)
}