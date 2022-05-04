#' @include LabKeySettings.R
#' @include Utils.R

#' @import dplyr

#' @title DownloadAndAppendNimble
#' @description This read a given seurat object and download/append all associated files to the object
#'
#' @param seuratObj A Seurat object.
#' @param outPath The path to which nimble files will be downloaded and saved
#' @param enforceUniqueFeatureNames Whether or not to fail if we discover multiple nimble files with the same feature name
#' @param ensureSamplesShareAllGenomes If true, the function will fail unless all samples have data from the same set of genomes
#' @return A modified Seurat object.
#' @export
DownloadAndAppendNimble <- function(seuratObject, outPath=tempdir(), enforceUniqueFeatureNames=TRUE, allowableGenomes=NULL, ensureSamplesShareAllGenomes = TRUE) {
  # Ensure we have a DatasetId column
  if (is.null(seuratObject@meta.data[['DatasetId']])) {
    stop('Seurat object lacks DatasetId column')
  }
  
  # Produce a nimble file id/DatasetId vector for each DatasetId
  nimbleFileComponents <- list()
  genomeToDataset <- list()
  for (datasetID in unique(seuratObject@meta.data[['DatasetId']])) {
    print(paste0('Possibly adding nimble data for dataset: ', datasetID))
    
    nimbleToGenome <- .queryNimble(loupeDataId=datasetID, allowableGenomes=allowableGenomes)

    # TODO Make this configurable whether or not we skip or return
    if (is.null(nimbleIds)) {
      print(paste0('Nimble file(s) not found for dataset: ', datasetID, ', skipping'))
      genomeToDataset[[datasetID]] <- integer()

      next
    }

    genomeToDataset[[datasetID]] <- names(nimbleToGenome)

    nimbleFiles <- character()
    for (genomeId in names(nimbleToGenome)) {
      outputFileId <- nimbleToGenome[[genomeId]]
      nimbleFile <- file.path(outPath, paste0('nimbleCounts.', datasetId, '.', genomeId, '.', id, '.tsv'))
      DownloadOutputFile(outputFileId = id, outFile = nimbleFile, overwrite = T)
      if (!file.exists(nimbleFile)) {
        stop(paste0('Unable to download calls table for nimbleId: ', nimbleId, ' datasetId: ', datasetId))
      }

      nimbleFiles <- c(nimbleFiles, nimbleFile)
    }

    nimbleFileComponents[[datasetID]] <- nimbleFiles
  }

  # Ensure genomes identical across samples:
  if (ensureSamplesShareAllGenomes) {
    consensusGenomes <- NULL
    for (datasetId in names(genomeToDataset)) {
      if (all(is.null(consensusGenomes))) {
        consensusGenomes <- genomeToDataset[[datasetId]]
      } else {
        if (length(consensusGenomes) != length(genomeToDataset[[datasetId]]) || length(setdiff(genomeToDataset[[datasetId]], consensusGenomes)) > 0) {
           stop(paste0('Unequal genomes found between datasets. Current set: ', paste0(consensusGenomes, collapse = ','), ', but ', datasetId, ' uses: ', paste0(genomeToDataset[[datasetId]], collapse = ',')))
        }
      }
    }
  }

  df <- .mergeNimbleFiles(fileComponents=nimbleFileComponents, enforceUniqueFeatureNames)
  
  outFile <- file.path(outPath, paste0("mergedNimbleCounts.tsv"))
  write.table(df, outFile, sep="\t", col.names=F, row.names=F, quote=F)
  
  seuratObject <- CellMembrane::AppendNimbleCounts(seuratObject=seuratObject, file=outFile)
  unlink(outFile)

  return(seuratObject)
}

.queryNimble <- function(loupeDataId, allowableGenomes=NULL) {
  rows <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
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

  filterClauses <- list(c("readset", "EQUAL", readset),
                        c("category", "EQUAL", "Nimble Alignment")
  )
  if (!is.null(allowableGenomes)) {
    filterClauses <- append(filterClauses, list(c("library_id", "IN", paste0(allowableGenomes, collapse = ';'))))
  }

  filter <- do.call(makeFilter, filterClauses)

  rows <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    colSort="-rowid",
    colSelect="rowid,library_id",
    colFilter=filter,
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) == 0) {
    print(paste0("No nimble objects found for readset: ", readset))
    return(NULL)
  } else {
    df <- data.frame(rowid = rows$rowid, library_id = rows$library_id)

    # TODO: consider failing if expected library_ids are not present

    # Always take the most recent output by genome:
    groups <- df %>% group_by(library_id) %>% summarise(rowid = max(rowid))
    ret <- groups$rowid
    names(ret) <- groups$library_id

    return(ret)
  }
}

.mergeNimbleFiles <- function(fileComponents, enforceUniqueFeatureNames=TRUE) {
  df <- NULL
  
  for (datasetId in names(fileComponents)) {
    files <- fileComponents[[datasetId]]
    for (fn in files) {
      if (!file.exists(fn)) {
        stop(paste0('Unable to open local file for nimbleId: ', fn, ' datasetId: ', datasetId))
      }

      nimbleTable <- read.table(fn, sep="\t", header=FALSE)
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
  }

  return(df)
}