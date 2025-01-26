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
#' @param maxAmbiguityAllowed If provided, any features representing more than ths value will be discarded. For example, 'Feat1,Feat2,Feat3' represents 3 features. maxAmbiguityAllowed=1 results in removal of all ambiguous features.
#' @param reuseExistingDownloads If true, any pre-existing downloaded nimble TSVs will be re-used
#' @param performDietSeurat If true, DietSeurat will be run, which removes existing reductions. This may or may not be required based on your usage.
#' @param normalizeData Passed directly to AppendNimbleCounts()
#' @param assayForLibrarySize Passed directly to AppendNimbleCounts()
#' @param maxLibrarySizeRatio Passed directly to AppendNimbleCounts()
#' @param queryDatabaseForLineageUpdates If true, after downloading the raw nimble output, the code will query any feature not ending with 'g' against the database and replace that name with the current value of lineage.
#' @param replaceExistingAssayData If true, any existing data in the targetAssay will be deleted
#' @return A modified Seurat object.
#' @export
DownloadAndAppendNimble <- function(seuratObject, targetAssayName, outPath=tempdir(), enforceUniqueFeatureNames=TRUE, allowableGenomes=NULL, ensureSamplesShareAllGenomes = TRUE, maxAmbiguityAllowed = 1, reuseExistingDownloads = FALSE, performDietSeurat = FALSE, normalizeData = TRUE, assayForLibrarySize = 'RNA', maxLibrarySizeRatio = 0.05, queryDatabaseForLineageUpdates = FALSE, replaceExistingAssayData = TRUE) {
  # Ensure we have a DatasetId column
  if (is.null(seuratObject@meta.data[['DatasetId']])) {
    stop('Seurat object lacks DatasetId column')
  }

  if (performDietSeurat) {
    print('Running DietSeurat')
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
        stop(paste0('Nimble file(s) not found for datasetId: ', datasetId, ', for genome(s): ', paste0(allowableGenomes, collapse = ';')))
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
      nimbleFile <- file.path(outPath, paste0('nimbleCounts.', datasetId, '.', outputFileId, '.', genomeId, '.tsv'))
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
  df <- .mergeNimbleFiles(fileComponents=nimbleFileComponents, enforceUniqueFeatureNames, queryDatabaseForLineageUpdates = queryDatabaseForLineageUpdates)
  print(paste0('Total features: ', length(unique(df$V1))))
  
  outFile <- file.path(outPath, paste0("mergedNimbleCounts.tsv"))
  write.table(df, outFile, sep="\t", col.names=F, row.names=F, quote=F)

  print(paste0('Appending counts to ', targetAssayName))
  seuratObject <- AppendNimbleCounts(seuratObject=seuratObject, targetAssayName = targetAssayName, nimbleFile=outFile, maxAmbiguityAllowed = maxAmbiguityAllowed, performDietSeurat = FALSE, normalizeData = normalizeData, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = maxLibrarySizeRatio, replaceExistingAssayData = replaceExistingAssayData)
  unlink(outFile)

  return(seuratObject)
}

.queryNimble <- function(loupeDataId, allowableGenomes=NULL) {
  rows <- suppressWarnings(labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    colSort="-rowid",
    colSelect="readset",
    colFilter=makeFilter(c("rowid", "EQUAL", loupeDataId)),
    containerFilter=NULL,
    colNameOpt="rname"
  ))
  
  if (nrow(rows) == 0) {
    translated <- .ResolveLoupeIdFromDeleted(loupeDataId, throwOnError = FALSE)
    if (all(is.null(translated))) {
      print(paste0("Loupe File ID: ", loupeDataId, " not found"))
      return(NA)
    }

    rows <- translated[,'readsetid',drop = FALSE]
    names(rows) <- c('readset')
  }
  
  readset <- unique(rows[['readset']])
  
  if (is.na(readset) || is.null(readset)) {
    message("Readset is NA/NULL for loupe file: ", loupeDataId)
    return(NA)
  }

  filterClauses <- list(c("readset", "EQUAL", readset),
                        c("category", "EQUAL", "Nimble Results")
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

.mergeNimbleFiles <- function(fileComponents, enforceUniqueFeatureNames = TRUE, enforceCellBarcodeFormat = TRUE, queryDatabaseForLineageUpdates = FALSE) {
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

      # strip hyphens
      nimbleTable$V3 <- gsub(nimbleTable$V3, pattern = '-[0-9]$', replacement = '')

      if (enforceCellBarcodeFormat) {
        suspiciousRows <- !grepl(nimbleTable$V3, pattern = '^[ATGCN]+$')
        if (sum(suspiciousRows) > 0) {
          fn <- paste0('suspiciousCellBarcodes.', datasetId, '.txt')
          nimbleTable$RowNumber <- 1:nrow(nimbleTable)
          badRows <- nimbleTable[suspiciousRows,]
          write.table(badRows, file = fn, sep = '\t', row.names = FALSE, quote = FALSE)
          stop(paste0('The dataset', datasetId, ' had suspicious cell barcodes. Rows written to: ', fn))
        }
      }

      if (nrow(nimbleTable) > 0) {
        nimbleTable$V3 <- paste0(datasetId, "_", nimbleTable$V3)
      }

      # Note: this was an old nimble bug that has been fixed, but retain this check:
      if (sum(nimbleTable$V1 == "") > 0) {
        stop(paste0("The nimble data contains blank feature names. This should not occur. Dataset Id: ", datasetId))
      }

      if (sum(grepl(nimbleTable$V1, pattern = "^,")) > 0) {
        stop(paste0("The nimble data contains features with leading commas. This should not occur. Dataset Id: ", datasetId))
      }

      if (sum(grepl(nimbleTable$V1, pattern = ",$")) > 0) {
        stop(paste0("The nimble data contains features with trailing commas. This should not occur. Dataset Id: ", datasetId))
      }

      if (queryDatabaseForLineageUpdates) {
        genomeId <- unlist(strsplit(basename(fn), split = '\\.'))[4]
        df <- .UpdateLineages(df, libraryId = genomeId)
      }

      nimbleTableGrouped <- nimbleTable %>% group_by(V1, V3) %>% summarize(V2 = sum(V2), InputRows = n())
      if (sum(nimbleTableGrouped$InputRows > 1) > 0) {
        stop(paste0('There were duplicate rows from the same cell-barcode/feature in the input for dataset. This could occur if one job appended to the input file: ', datasetId))
      }
      nimbleTable <- nimbleTable[c('V1', 'V2', 'V3')]

      if (enforceCellBarcodeFormat) {
        suspiciousRows <- !grepl(nimbleTable$V3, pattern = '^[0-9]+_[ATGCN]+$')
        if (sum(suspiciousRows) > 0) {
          fn <- paste0('suspiciousCellBarcodes.', datasetId, '.txt')
          badRows <- nimbleTable[suspiciousRows,]
          badRows$RowNumber <- 1:nrow(nimbleTable)[suspiciousRows]
          write.table(badRows, file = fn, sep = '\t', row.names = FALSE, quote = FALSE)
          stop(paste0('The dataset', datasetId, ' had suspicious cell barcodes. Rows written to: ', fn))
        }
      }

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

.UpdateLineages <- function(df, libraryId) {
  print(paste0('Querying DB to update lineages for: ', libraryId))

  feats <- unique(df$V1)
  feats <- feats[!grepl(feats, pattern = 'g$')]
  if (length(feats) == 0) {
    return(df)
  }

  dat <- suppressWarnings(labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber",
    schemaName="sequenceanalysis",
    queryName="reference_library_members",
    colSelect="ref_nt_id/name,ref_nt_id/lineage",
    colSort="ref_nt_id/name",
    colFilter=makeFilter(
      c("library_id", "EQUALS", libraryId),
      c("ref_nt_id/lineage", "NOT_MISSING", ''),
      c("ref_nt_id/name", "IN", paste0(feats, collapse = ';'))
    ),
    colNameOpt="rname"
  ))
  names(dat) <- c('name', 'lineage')

  if (nrow(dat) == 0){
    return(df)
  }

  print(paste0('Total features with updated lineages: ', nrow(dat)))
  df$V1 <- as.character(df$V1)
  for (idx in seq_along(dat$name)) {
    df$V1[df$V1 == dat$name[idx]] <- dat$lineage[idx]
  }

  df <- df %>%
    group_by(V1, V3) %>%
    summarize(V2 = sum(V2))

  print(paste0('Distinct features after re-grouping: ', length(unique(df$V1))))

  return(df)
}