#' @include LabKeySettings.R
#' @include Utils.R

library(dplyr)

#' @title DownloadAndAppendNimble
#' @description This read a given seurat object and download/append all associated files to the object
#'
#' @param seuratObj A Seurat object.
#' @param outPath The path to download nimble files to
#' @param enforceUniqueFeatureNames Whether or not to fail if we discover multiple nimble files with the same feature name
#' @return A modified Seurat object.
#' @export
DownloadAndAppendNimble <- function(seuratObject, outPath=tempdir(), enforceUniqueFeatureNames=TRUE) {#, genomes=NULL) {
  
  # Ensure we have a DatasetId column
  if (is.null(seuratObject@meta.data[['DatasetId']])) {
    stop('Seurat object lacks DatasetId column')
  }
  
  # Produce a nimble file id/DatasetId vector for each DatasetId
  nimbleFileComponents <- list()
  for (datasetID in unique(seuratObject@meta.data[['DatasetId']])) {
    print(paste0('Possibly adding nimble data for dataset: ', datasetID))
    
    nimbleIds <- .queryNimble(loupeDataId=datasetID)#, genomeFilter=genomes)
    
    # TODO Make this configurable whether or not we skip or return
    if (is.null(nimbleIds)) {
      print(paste0('Nimble file(s) not found for dataset: ', datasetID, ', skipping'))
      next
    }
    
    for(id in nimbleIds) {
      nimbleFile <- file.path(outPath, paste0(id, '_nimbleCounts.tsv'))
      
      DownloadOutputFile(outputFileId = id, outFile = nimbleFile, overwrite = T)
      
      if (!file.exists(nimbleFile)) {
        stop(paste0('Unable to download calls table for nimbleId: ', nimbleId, ' datasetId: ', datasetId))
      }
      
      nimbleFileComponents <- append(nimbleFileComponents, list(c(nimbleFile, datasetID)))
    }
  }
  
  print(nimbleFileComponents)
  df <- .mergeNimbleFiles(fileComponents=nimbleFileComponents, enforceUniqueFeatureNames)
  
  outFile <- file.path(outPath, paste0("mergedNimbleCounts.tsv"))
  
  print(paste0("Writing combined nimble file to disk at: ", outFile))
  write.table(df, outFile, sep="\t", col.names=F, row.names=F, quote=F)
  
  seuratObject <- CellMembrane::AppendNimbleCounts(seuratObject=seuratObject, file=outFile)
  return(seuratObject)
}

#' @title queryNimble
#' @description Given a loupe data Id, get the associated nimble files
#'
#' @param loupeDataId A Loupe data id
.queryNimble <- function(loupeDataId) {
  rows <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="readset",
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
    
  rows <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    colSort="-rowid",
    colSelect="rowid,library_id",
    colFilter=makeFilter(c("readset", "EQUAL", readset),
                         c("category", "EQUAL", "Nimble Alignment")),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  
  ret <- NULL
  if (nrow(rows) == 0) {
    print(paste0("No nimble objects found for readset: ", readset))
  } else {
    df <- data.frame(rows$rowid, rows$library_id)
    groups <- df %>% group_by(rows.library_id) %>% filter(row_number()==1)
    ret <- groups$rows.rowid
  }
  
  return(ret)
}

#' @title .mergeNimbleFiles
#' @description Takes a set of nimble file/dataID tuples and appends them into one file
#' @param enforceUniqueFeatureNames Whether or not to fail if we discover multiple nimble files with the same feature name
#' @param fileComponents A vector of nimble file/dataID tuples
.mergeNimbleFiles <- function(fileComponents=c(), enforceUniqueFeatureNames=TRUE) {
  df <- NULL
  
  for(component in fileComponents) {
    file <- component[1]
    datasetId <- component[2]
    
    if (!file.exists(file)) {
      stop(paste0('Unable to open local file for nimbleId: ', file, ' datasetId: ', datasetId))
    }
    
    nimbleTable <- read.table(file, sep="\t", header=FALSE)
    nimbleTable$V3 <- paste0(datasetId, "_", nimbleTable$V3)
    
    if(is.null(df)) {
      df <- nimbleTable
    } else {
      dfFeatures <- unique(df$'V1')
      tableFeatures <- unique(nimbleTable$'V1')
      sharedFeatures <- tableFeatures %in% dfFeatures
      
      if(any(sharedFeatures) && enforceUniqueFeatureNames) {
        stop(paste0("Cannot merge nimble files: features shared between libraries: ", dfFeatures$V1[sharedFeatures]))  
      }
      
      df <- rbind(df, nimbleTable)
    }
  }
  return(df)
}