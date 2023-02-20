#' @include LabKeySettings.R
#' @include Utils.R

utils::globalVariables(
  names = c('V1','V2','V3','library_id','rowid'),
  package = 'Rdiscvr',
  add = TRUE
)

#' @import dplyr

#' @title DownloadAndAppendNimble
#' @description This read a given seurat object and download/append all associated files to the object
#'
#' @param seuratObject A Seurat object.
#' @param targetAssayName The target assay. If this assay exists, features will be appended (and an error thrown if there are duplicates). Otherwise a new assay will be created.
#' @param outPath The path to which nimble files will be downloaded and saved
#' @param enforceUniqueFeatureNames Whether or not to fail if we discover multiple nimble files with the same feature name
#' @param allowableGenomes An optional vector of genomeIds. If provided, nimble results from these genomes will be appended. If any dataset lacks a nimble file for a genome, it will fail
#' @param ensureSamplesShareAllGenomes If true, the function will fail unless all samples have data from the same set of genomes
#' @param dropAmbiguousFeatures If true, any ambiguous feature (identified as containing a comma)
#' @param reuseExistingDownloads If true, any pre-existing downloaded nimble TSVs will be re-used
#' @param performDietSeurat If true, DietSeurat will be run, which removes existing reductions. This may or may not be required based on your usage, but the default is true out of caution.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendNimble <- function(seuratObject, targetAssayName, outPath=tempdir(), enforceUniqueFeatureNames=TRUE, allowableGenomes=NULL, ensureSamplesShareAllGenomes = TRUE, dropAmbiguousFeatures = TRUE, reuseExistingDownloads = FALSE, performDietSeurat = TRUE) {
  # Ensure we have a DatasetId column
  if (is.null(seuratObject@meta.data[['DatasetId']])) {
    stop('Seurat object lacks DatasetId column')
  }

  if (performDietSeurat) {
    seuratObject <- Seurat::DietSeurat(seuratObject)
  }

  # Produce a nimble file id/DatasetId vector for each DatasetId
  nimbleFileComponents <- list()
  genomeToDataset <- list()
  print(paste0('Total datasets: ', length(unique(seuratObject@meta.data[['DatasetId']]))))
  for (datasetId in unique(seuratObject@meta.data[['DatasetId']])) {
    print(paste0('Possibly adding nimble data for dataset: ', datasetId))
    
    nimbleToGenome <- .queryNimble(loupeDataId=datasetId, allowableGenomes=allowableGenomes)

    if (is.null(nimbleToGenome)) {
      if (ensureSamplesShareAllGenomes) {
        stop(paste0('Nimble file(s) not found for dataset: ', datasetId, ', for genome(s): ', paste0(allowableGenomes, collapse = ';')))
      } else {
        print(paste0('Nimble file(s) not found for dataset: ', datasetId, ', skipping'))
        genomeToDataset[[datasetId]] <- integer()
      }

      next
    }

    genomeToDataset[[datasetId]] <- names(nimbleToGenome)

    nimbleFiles <- character()
    for (genomeId in names(nimbleToGenome)) {
      outputFileId <- nimbleToGenome[[genomeId]]
      nimbleFile <- file.path(outPath, paste0('nimbleCounts.', datasetId, '.', genomeId, '.tsv'))
      DownloadOutputFile(outputFileId = outputFileId, outFile = nimbleFile, overwrite = !reuseExistingDownloads)
      if (!file.exists(nimbleFile)) {
        stop(paste0('Unable to download calls table for genome: ', genomeId, ' datasetId: ', datasetId))
      }

      nimbleFiles <- c(nimbleFiles, nimbleFile)
    }

    nimbleFileComponents[[datasetId]] <- nimbleFiles
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

  print('Merging into single matrix')
  df <- .mergeNimbleFiles(fileComponents=nimbleFileComponents, enforceUniqueFeatureNames)
  print(paste0('Total features: ', length(unique(df$V1))))
  
  outFile <- file.path(outPath, paste0("mergedNimbleCounts.tsv"))
  write.table(df, outFile, sep="\t", col.names=F, row.names=F, quote=F)

  print(paste0('Appending counts to ', targetAssayName))
  seuratObject <- AppendNimbleCounts(seuratObject=seuratObject, targetAssayName = targetAssayName, nimbleFile=outFile, dropAmbiguousFeatures = dropAmbiguousFeatures, performDietSeurat = FALSE)
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
    message(paste0("Loupe File ID: ", loupeDataId, " not found"))
    return(NA)
  }
  
  readset <- unique(rows[['readset']])
  
  if (is.na(readset) || is.null(readset)) {
    message("Readset is NA/NULL for loupe file: ", loupeDataId)
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
    message(paste0("No nimble objects found for readset: ", readset))
    return(NULL)
  } else {
    df <- data.frame(rowid = rows$rowid, library_id = rows$library_id)

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

    # Features shared across datasets is fine. The problem arises if two genomes from the same dataset do.
    featuresForDataset <- NULL
    for (fn in files) {
      print(paste0('Reading nimble TSV for: ', datasetId, ': ', fn))

      if (!file.exists(fn)) {
        stop(paste0('Unable to open local file for nimbleId: ', fn, ' datasetId: ', datasetId))
      }

      # detect/handle gzip:
      f <- file(fn)
      filetype <- summary(f)$class
      close.connection(f)
      if (filetype == 'gzfile') {
        message('Opening gzipped nimble output')
        f <- gzfile(fn)

        # Check for an empty file:
        numLines <- length(readLines(f))
        if (numLines == 0) {
          print(paste0('No lines in file, skipping: ', fn))
          close(f)
          nimbleTable <- data.frame(V1 = character(), V2 = integer(), V3 = character())
        } else {
          nimbleTable <- read.table(f, sep="\t", header=FALSE)
        }

        # Note: read.table seems to automatically close the gzfile connection, but check just in case:
        dat <- as.data.frame(showConnections(all = TRUE))
        if (fn %in% dat$description) {
          print('gzip connection was not closed!')
          close.connection(f)
        }
      } else {
        # NOTE: nimble could run on a dataset and produce a zero-line output if no passing alignments are found
        if (file.info(fn)$size == 0) {
          print(paste0('No lines in file, skipping: ', fn))
          nimbleTable <- data.frame(V1 = character(), V2 = integer(), V3 = character())
        } else {
          nimbleTable <- read.table(fn, sep="\t", header=FALSE)
        }
      }

      if (nrow(nimbleTable) > 0) {
        nimbleTable$V3 <- paste0(datasetId, "_", nimbleTable$V3)
      }

      # NOTE: this could occur if a job was restarted after a failure. Prior versions of nimble used append instead of overwrite for the output.
      # TODO: remove this eventually
      nimbleTableGrouped <- nimbleTable %>% group_by(V1, V3) %>% summarize(V2 = sum(V2), InputRows = n())
      if (sum(nimbleTableGrouped$InputRows > 1) > 0) {
        warning(paste0('There were duplicate rows from the same cell-barcode/feature in the input for dataset. This could occur if one job overwrote an existing input: ', datasetId))
        print(paste0('Rows: ', nrow(nimbleTable)))
        print(paste0('Unique Rows: ', nrow(unique(nimbleTable))))

        # The rationale is that if the full dataframe is an even multiple of the unique rows, this was likely re-writting over the input file
        if (nrow(nimbleTable) %% nrow(unique(nimbleTable)) == 0) {
          print('Keeping unique rows')
          nimbleTable <- unique(nimbleTable)
        }
      }

      # TODO: this should also be removed eventually, since this nimble bug was fixed. This code exists to allow legacy files to be processed:
      if (sum(nimbleTable$V1 == "") > 0) {
        warning("The nimble data contains blank feature names. This should not occur. They will be removed")
        nimbleTable <- nimbleTable[nimbleTable$V1 != "",]
        nimbleTable <- nimbleTable %>% group_by(V1, V3) %>% summarize(V2 = sum(V2))
        nimbleTable <- nimbleTable[c('V1', 'V2', 'V3')]
      }

      if (sum(grepl(nimbleTable$V1, pattern = "^,")) > 0) {
        warning("The nimble data contains features with leading commas. This should not occur. They will be removed")
        nimbleTable$V1 <- as.character(nimbleTable$V1)
        nimbleTable$V1 <- gsub(nimbleTable$V1, pattern = "^,", replacement = "")
        nimbleTable <- nimbleTable %>% group_by(V1, V3) %>% summarize(V2 = sum(V2))
        nimbleTable <- nimbleTable[c('V1', 'V2', 'V3')]
      }

      if (sum(grepl(nimbleTable$V1, pattern = ",$")) > 0) {
        warning("The nimble data contains features with trailing commas. This should not occur. They will be removed")
        nimbleTable$V1 <- as.character(nimbleTable$V1)
        nimbleTable$V1 <- gsub(nimbleTable$V1, pattern = ",$", replacement = "")
        nimbleTable <- nimbleTable %>% group_by(V1, V3) %>% summarize(V2 = sum(V2))
        nimbleTable <- nimbleTable[c('V1', 'V2', 'V3')]
      }

      nimbleTableGrouped <- nimbleTable %>% group_by(V1, V3) %>% summarize(V2 = sum(V2), InputRows = n())
      if (sum(nimbleTableGrouped$InputRows > 1) > 0) {
        stop(paste0('There were duplicate rows from the same cell-barcode/feature in the input for dataset. This could occur if one job appended to the input file: ', datasetId))
      }
      nimbleTable <- nimbleTable[c('V1', 'V2', 'V3')]

      if (is.null(df)) {
        df <- nimbleTable
        featuresForDataset <- sort(unique(nimbleTable$V1))
      } else {
        incomingFeatures <- sort(unique(nimbleTable$V1))
        sharedFeatures <- incomingFeatures %in% featuresForDataset

        if (sum(sharedFeatures) > 0 && enforceUniqueFeatureNames) {
            stop(paste0("Cannot merge nimble files: features shared between libraries: ", paste0(incomingFeatures[sharedFeatures], collapse = '; ')))
        }

        featuresForDataset <- sort(c(featuresForDataset, incomingFeatures))

        df <- rbind(df, nimbleTable)
      }
    }
  }

  return(df)
}