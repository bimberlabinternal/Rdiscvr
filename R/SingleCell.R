#' @include LabKeySettings.R
#' @include Utils.R

utils::globalVariables(
  names = c('SortOrder', 'BarcodePrefix', 'Total'),
  package = 'Rdiscvr',
  add = TRUE
)
#' @title QueryAndApplyCdnaMetadata
#' @description This will read the barcodePrefix from the seurat object, query the LK server and apply the appropriate metadata from the cDNA table
#'
#' @param seuratObj A Seurat object.
#' @param fieldSelect The set of fields to query
#' @param fieldNames The labels to use for the fields
#' @param overwriteExisting If true and the output exists, it will be overwritten
#' @return A modified Seurat object.
#' @export
#' @importFrom dplyr %>% group_by_at summarise_at arrange
QueryAndApplyCdnaMetadata <- function(seuratObj,
                                      fieldSelect = c('rowid', 'sortid/population', 'sortid/sampleid/subjectid', 'sortid/sampleid/sampledate', 'sortid/sampleid/stim', 'sortid/sampleid/tissue', 'plateid', 'workbook/workbookid', 'sortid/sampleid/assaytype'),
                                      fieldNames = c('cDNA_ID', 'Population', 'SubjectId', 'SampleDate', 'Stim', 'Tissue', 'PlateId', 'WorkbookId', 'AssayType'), overwriteExisting = TRUE) {
  if (length(fieldSelect) != length(fieldNames)) {
    stop('The length of fields must equal the length of fieldNames')
  }

  if (sum(duplicated(fieldSelect)) > 0) {
    stop('All values for fieldSelect must be unique')
  }

  httr::set_config(httr::timeout(60000))

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
  prefixes <- unique(seuratObj$BarcodePrefix)
  outputFiles <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    colSort="-rowid",
    colFilter=makeFilter(c('rowId', "IN", paste0(prefixes, collapse = ';'))),
    colSelect="rowid,readset",
    containerFilter=NULL,
    colNameOpt="rname"
  )

  # resolve using deleted
  prefixNotFound <- prefixes[!prefixes %in% outputFiles$rowid]
  if (length(prefixNotFound) > 0) {
    print(paste0('The following barcodePrefix values were not found, looking for deleted records: ', paste0(prefixNotFound, collapse = ';')))
    translatedRows <- suppressWarnings(labkey.selectRows(
      baseUrl=.getBaseUrl(),
      folderPath=.getLabKeyDefaultFolder(),
      schemaName="singlecell",
      queryName="singlecellDatasets",
      colFilter = makeFilter(c("loupeFileId", "IN", paste0(prefixNotFound, collapse = ";"))),
      colSelect="loupeFileId,readsetId",
      containerFilter=NULL,
      colNameOpt="rname"
    ))

    if (nrow(translatedRows) > 0) {
      names(translatedRows) <- c('rowid', 'readset')
      print(paste0('The following IDs were matched to deleted objects: ', paste0(translatedRows$rowid, collapse = ';')))
      outputFiles <- rbind(outputFiles, translatedRows)
      prefixNotFound <- prefixes[!prefixes %in% outputFiles$rowid]
    }

    if (length(prefixNotFound) > 0) {
      stop('The following barcodePrefix values were not found, even after looking for deleted records: ', paste0(prefixNotFound, collapse = ','))
    }
  }

  names(outputFiles) <- c('BarcodePrefix', readsetLabel)
  outputFiles$BarcodePrefix <- as.character(outputFiles$BarcodePrefix)
  print(paste0('total outputfile rows: ', nrow(outputFiles)))
  if (nrow(outputFiles) != length(unique(seuratObj$BarcodePrefix))) {
    missing <- sort(c(setdiff(unique(seuratObj$BarcodePrefix), unique(outputFiles$BarcodePrefix)), setdiff(unique(outputFiles$BarcodePrefix), unique(seuratObj$BarcodePrefix))))
    warning(paste0('Did not find output file record for all prefixes!  Missing: ', paste0(missing, collapse = ',')))
  }

  rows <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="singlecell",
    queryName="cdna_libraries",
    viewName="",
    colSort="-rowid",
    colFilter = makeFilter(c("readsetId", "IN", paste0(unique(outputFiles[[readsetLabel]]), collapse = ";"))),
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

  if ('cDNA_ID' %in% names(rows) && any(grepl(rows$cDNA_ID, pattern = ','))){
    dups <- sort(unique(grep(rows$cDNA_ID, pattern = ',', value = TRUE)))
    stop(paste0('The data contained rows matching multiple cDNA_IDs: ', paste0(dups, collapse = '; ')))
  }

  origBarcodes <- colnames(seuratObj)
  hasHTO <- !all(is.na(rows[htoLabel]))
  if (hasHTO) {
    if (!('HTO' %in% names(seuratObj@meta.data))) {
      print('Adding HTO column to seurat object')
      seuratObj$HTO <- NA
    }

    if (!('HTO.Classification' %in% names(seuratObj@meta.data))) {
      print('Adding HTO.Classification column to seurat object')
      seuratObj$HTO.Classification <- NA
    }

    df <- data.frame(HTO = as.character(seuratObj$HTO), BarcodePrefix = as.character(seuratObj$BarcodePrefix), Barcode = origBarcodes, SortOrder = seq_along(origBarcodes))

    #Allow for libraries that have a non-null HTO, but have only a single per library (which is effectively the same as not being hashed):
    rows2 <- rows %>% dplyr::group_by(BarcodePrefix) %>% dplyr::summarise(Total = dplyr::n_distinct(HTO)) %>% dplyr::filter(Total == 1)
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
    df <- data.frame(HTO = NA, BarcodePrefix = as.character(seuratObj$BarcodePrefix), Barcode = origBarcodes, SortOrder = seq_along(origBarcodes))
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
  for (idx in seq_along(fieldSelect)) {
    fieldName <- fieldNames[idx]
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

    # Dont allow empty strings:
    v <- df[[fieldName]]
    v[v == ''] <- NA
    names(v) <- rownames(df)
    seuratObj[[fieldName]] <- v
  }

  return(seuratObj)
}

#' @title QueryAndApplyMetadataUsingCDNA
#' @description This will read the cDNA_ID column from the seurat object, query the LK server and apply the appropriate metadata from the cDNA table
#'
#' @param seuratObj A Seurat object.
#' @param fieldSelect The set of fields to query
#' @param fieldNames The labels to use for the fields
#' @param overwriteExisting If true and the output exists, it will be overwritten
#' @return A modified Seurat object.
#' @export
#' @importFrom dplyr %>% group_by_at summarise_at arrange
QueryAndApplyMetadataUsingCDNA <- function(seuratObj,
                                      fieldSelect = c('rowid', 'sortid/population', 'sortid/sampleid/subjectid', 'sortid/sampleid/sampledate', 'sortid/sampleid/stim', 'sortid/sampleid/tissue', 'plateid', 'workbook/workbookid', 'sortid/sampleid/assayType'),
                                      fieldNames = c('cDNA_ID', 'Population', 'SubjectId', 'SampleDate', 'Stim', 'Tissue', 'PlateId', 'WorkbookId', 'AssayType'), overwriteExisting = TRUE) {
  if (length(fieldSelect) != length(fieldNames)) {
    stop('The length of fields must equal the length of fieldNames')
  }

  if (sum(duplicated(fieldSelect)) > 0) {
    stop('All values for fieldSelect must be unique')
  }

  if (!'cDNA_ID' %in% names(seuratObj@meta.data)) {
    stop('Missing required metadata field: cDNA_ID')
  }

  #spike in readset, since we need this to join dataframes
  fieldSelect <- tolower(fieldSelect)
  if (!('rowid' %in% fieldSelect)) {
    print('Missing rowid/cDNA_ID to the field select')
    fieldSelect <- c('rowid', fieldSelect)
    fieldNames <- c('cDNA_ID', fieldNames)
  }

  rows <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="singlecell",
    queryName="cdna_libraries",
    viewName="",
    colSort="-rowid",
    colFilter = makeFilter(c("rowid", "IN", paste0(unique(seuratObj$cDNA_ID), collapse = ";"))),
    colSelect=paste0(fieldSelect, collapse = ','),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  names(rows) <- fieldNames
  print(paste0('Total cDNA rows: ', nrow(rows)))

  df <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj), SortOrder = 1:ncol(seuratObj))
  df <- merge(df, rows, by = 'cDNA_ID', all.x = TRUE)
  rownames(df) <- df$CellBarcode
  df <- df %>% arrange(SortOrder)
  df <- df[,!names(df) %in% c('SortOrder', 'CellBarcode', 'cDNA_ID'), drop = FALSE]

  if (nrow(df) != ncol(seuratObj)) {
    stop('Length of original seurat object and metadata not equal. Something went wrong merging')
  }

  if (sum(colnames(seuratObj) != rownames(df)) > 0) {
    stop('CellBarcode does not match for all cells.  Something went wrong merging')
  }

  #strings to factors:
  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.factor)

  #now apply the result:
  for (fieldName in names(df)) {
    if (!overwriteExisting && (fieldName %in% names(seuratObj@meta.data))) {
      print(paste0('column already exists, skipping: ', fieldName))
      next
    }

    # Dont allow empty strings:
    v <- df[[fieldName]]
    v[v == ''] <- NA
    names(v) <- rownames(df)
    seuratObj[[fieldName]] <- v
  }

  return(seuratObj)
}

