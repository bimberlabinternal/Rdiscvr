#' @include LabKeySettings.R
#' @include Utils.R

utils::globalVariables(
  names = c('SortOrder'),
  package = 'Rdiscvr',
  add = TRUE
)
#' @title QueryAndApplyCdnaMetadata
#' @description This will read the barcodePrefix from the seurat object, query the LK server and apply the appropriate metadata from the cDNA table
#'
#' @param seuratObj, A Seurat object.
#' @param fieldSelect The set of fields to query
#' @param fieldNames The labels to use for the fields
#' @param overwriteExisting If true and the output exists, it will be overwritten
#' @return A modified Seurat object.
#' @export
#' @importFrom dplyr %>% group_by_at summarise_at arrange
QueryAndApplyCdnaMetadata <- function(seuratObj,
                                      fieldSelect = c('rowid', 'sortid/population', 'sortid/stimid/animalId', 'sortid/stimid/date', 'sortid/stimid/stim', 'sortid/stimid/tissue'),
                                      fieldNames = c('cDNA_ID', 'Population', 'SubjectId', 'SampleDate', 'Stim', 'Tissue'), overwriteExisting = F) {
  if (length(fieldSelect) != length(fieldNames)) {
    stop('The length of fields must equal the length of fieldNames')
  }

  if (sum(duplicated(fieldSelect)) > 0) {
    stop('All values for fieldSelect must be unique')
  }

  #spike in readset, since we need this to join dataframes
  fieldSelect <- tolower(fieldSelect)
  readsetAdded <- F
  if (!('readsetid' %in% fieldSelect)) {
    readsetAdded <- T
    fieldSelect <- c('readsetid', fieldSelect)
    fieldNames <- c('Readset', fieldNames)
  }
  readsetIdx <- which('readsetid' %in% fieldSelect)
  readsetLabel <- fieldNames[readsetIdx]

  #also HTO:
  htoAdded <- F
  if (!('sortid/hto' %in% fieldSelect)) {
    htoAdded <- T
    fieldSelect <- c('sortid/hto', fieldSelect)
    fieldNames <- c('HTO', fieldNames)
  }
  htoIdx <- which('sortid/hto' %in% fieldSelect)
  htoLabel <- fieldNames[htoIdx]

  #Download info, based on BarcodePrefix:
  outputFiles <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colFilter=makeFilter(c('rowId', "IN", paste0(unique(seuratObj$BarcodePrefix), collapse = ';'))),
    colSelect="rowid,readset",
    containerFilter=NULL,
    colNameOpt="rname"
  )
  names(outputFiles) <- c('BarcodePrefix', readsetLabel)
  outputFiles$BarcodePrefix <- as.character(outputFiles$BarcodePrefix)
  print(paste0('total outputfile rows: ', nrow(outputFiles)))
  if (nrow(outputFiles) != length(unique(seuratObj$BarcodePrefix))) {
    missing <- sort(c(setdiff(unique(seuratObj$BarcodePrefix), unique(outputFiles$BarcodePrefix)), setdiff(unique(outputFiles$BarcodePrefix), unique(seuratObj$BarcodePrefix))))
    warning(paste0('Did not find output file record for all prefixes!  Missing: ', paste0(missing, collapse = ',')))
  }

  rows <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="tcrdb",
    queryName="cdnas",
    viewName="",
    colSort="-rowid",
    colFilter = makeFilter(c("readsetId", "IN", paste0(unique(outputFiles$Readset), collapse = ";"))),
    colSelect=paste0(fieldSelect, collapse = ','),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  names(rows) <- fieldNames
  print(paste0('total cDNA rows: ', nrow(rows)))

  rows <- merge(rows, outputFiles, by = c(readsetLabel), all.x = T)
  rows <- unique(rows)

  #Force unique values
  groupVars <- c(htoLabel, 'BarcodePrefix')
  colSummary <- fieldNames[!(fieldNames %in% groupVars)]
  rows <- rows %>% group_by_at(groupVars) %>% summarise_at(colSummary, function(x){
    paste(sort(unique(x)), collapse=',')
  })

  origBarcodes <- colnames(seuratObj)
  hasHTO <- !all(is.na(rows[htoLabel]))
  if (hasHTO) {
    if (!('HTO' %in% names(seuratObj@meta.data))) {
      print('Adding HTO columns to seurat object')
      seuratObj$HTO <- c(NA)
      seuratObj$HTO_Classification <- c(NA)
    }

    df <- data.frame(HTO = as.character(seuratObj$HTO), BarcodePrefix = as.character(seuratObj$BarcodePrefix), Barcode = origBarcodes, SortOrder = 1:length(origBarcodes))

    #Allow for libraries that have a non-null HTO, but have only a single per library (which is effectively the same as not being hashed):
    rows2 <- rows %>% group_by(BarcodePrefix) %>% summarise(Total = dplyr::n_distinct(HTO)) %>% filter(Total == 1)
    if (nrow(rows2) > 0) {
      df$HTO <- as.character(df$HTO)
      for (bc in unique(rows2$BarcodePrefix)) {
        r <- rows$HTO[rows$BarcodePrefix == bc]
        df$HTO[df$BarcodePrefix == bc] <- r
      }
      df$HTO <- naturalsort::naturalfactor(df$HTO)
    }

    names(df) <- c(htoLabel, 'BarcodePrefix', 'Barcode', 'SortOrder')
    df <- merge(df, rows, by = c(htoLabel, 'BarcodePrefix'), all.x = T)
  } else {
    df <- data.frame(HTO = NA, BarcodePrefix = as.character(seuratObj$BarcodePrefix), Barcode = origBarcodes, SortOrder = 1:length(origBarcodes))
    names(df) <- c(htoLabel, 'BarcodePrefix', 'Barcode', 'SortOrder')
    df <- merge(df, rows, by = c('BarcodePrefix'), all.x = T)
  }

  df <- dplyr::arrange(df, SortOrder)
  df <- df[colnames(df) != 'SortOrder']
  rownames(df) <- df$Barcode
  df <- df[colnames(df) != 'Barcode']

  if (nrow(df) != ncol(seuratObj)) {
    stop('Length of original seurat object and metadata not equal. Something went wrong merging')
  }

  if (sum(origBarcodes != rownames(df)) > 0) {
    stop('BarcodePrefix does not match for all cells.  Something went wrong merging')
  }

  #strings to factors:
  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.factor)

  #now apply the result:
  for (idx in 1:length(fieldSelect)) {
    fieldName <- fieldNames[idx]
    fieldKey <- fieldSelect[idx]
    if (readsetAdded && fieldName == readsetLabel) {
      next
    }

    if (htoAdded && fieldName == htoLabel) {
      next
    }

    if (!overwriteExisting && (fieldName %in% names(seuratObj@meta.data))) {
      print(paste0('column already exists, skipping: ', fieldName))
      next
    }

    #print(paste0('Adding column: ', fieldName, ' (', fieldKey, ')'))
    v <- df[[fieldName]]
    names(v) <- names(df)
    seuratObj[[fieldName]] <- v
  }

  return(seuratObj)
}




