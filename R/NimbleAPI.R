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
#' @param seuratObj A Seurat object.
#' @param targetAssayName The target assay. If this assay exists, features will be appended (and an error thrown if there are duplicates). Otherwise a new assay will be created.
#' @param outPath The path to which nimble files will be downloaded and saved
#' @param enforceUniqueFeatureNames Whether or not to fail if we discover multiple nimble files with the same feature name
#' @param allowableGenomes An optional vector of genomeIds. If provided, nimble results from these genomes will be appended. If any dataset lacks a nimble file for a genome, it will fail
#' @param ensureSamplesShareAllGenomes If true, the function will fail unless all samples have data from the same set of genomes
#' @param maxAmbiguityAllowed If provided, any features representing more than ths value will be discarded. For example, 'Feat1,Feat2,Feat3' represents 3 features. maxAmbiguityAllowed=1 results in removal of all ambiguous features.
#' @param reuseExistingDownloads If true, any pre-existing downloaded nimble TSVs will be re-used
#' @param performDietSeurat If true, DietSeurat will be run, which removes existing reductions. This may or may not be required based on your usage.
#' @param normalizeData Passed directly to nimbleR::AppendNimbleCounts()
#' @param assayForLibrarySize Passed directly to nimbleR::AppendNimbleCounts()
#' @param maxLibrarySizeRatio Passed directly to nimbleR::AppendNimbleCounts()
#' @param queryDatabaseForLineageUpdates If true, after downloading the raw nimble output, the code will query any feature not ending with 'g' against the database and replace that name with the current value of lineage.
#' @param replaceExistingAssayData If true, any existing data in the targetAssay will be deleted
#' @param featureRenameList An optional named list in the format <OLD_NAME> = <NEW_NAME>. If any <OLD_NAME> are present, the will be renamed to <NEW_NAME>. The intention of this is to recover specific ambiguous classes.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendNimble <- function(seuratObj, targetAssayName, outPath=tempdir(), enforceUniqueFeatureNames=TRUE, allowableGenomes=NULL, ensureSamplesShareAllGenomes = TRUE, maxAmbiguityAllowed = 1, reuseExistingDownloads = FALSE, performDietSeurat = FALSE, normalizeData = TRUE, assayForLibrarySize = 'RNA', maxLibrarySizeRatio = 0.05, queryDatabaseForLineageUpdates = FALSE, replaceExistingAssayData = TRUE, featureRenameList = NULL) {
  # Ensure we have a DatasetId column
  if (is.null(seuratObj@meta.data[['DatasetId']])) {
    stop('Seurat object lacks DatasetId column')
  }

  if (performDietSeurat) {
    print('Running DietSeurat')
    seuratObj <- Seurat::DietSeurat(seuratObj)
  }

  # Produce a nimble file id/DatasetId vector for each DatasetId
  nimbleFileComponents <- list()
  genomeToDataset <- list()
  print(paste0('Total datasets: ', length(unique(seuratObj@meta.data[['DatasetId']]))))
  for (datasetId in unique(seuratObj@meta.data[['DatasetId']])) {
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
  df <- .mergeNimbleFiles(seuratObj, fileComponents=nimbleFileComponents, enforceUniqueFeatureNames, queryDatabaseForLineageUpdates = queryDatabaseForLineageUpdates)
  print(paste0('Total features: ', length(unique(df$V1))))
  
  outFile <- file.path(outPath, paste0("mergedNimbleCounts.tsv"))
  df <- df[c('V1', 'V2', 'V3')] # grouping could re-order this
  write.table(df, outFile, sep="\t", col.names=F, row.names=F, quote=F)

  print(paste0('Appending counts to ', targetAssayName))
  seuratObj <- nimbleR::AppendNimbleCounts(seuratObj = seuratObj, targetAssayName = targetAssayName, nimbleFile=outFile, maxAmbiguityAllowed = maxAmbiguityAllowed, performDietSeurat = FALSE, normalizeData = normalizeData, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = maxLibrarySizeRatio, replaceExistingAssayData = replaceExistingAssayData, featureRenameList = featureRenameList, doPlot = TRUE)
  unlink(outFile)

  return(seuratObj)
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
    names(rows) <- 'readset'
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

.mergeNimbleFiles <- function(seuratObj, fileComponents, enforceUniqueFeatureNames = TRUE, enforceCellBarcodeFormat = TRUE, queryDatabaseForLineageUpdates = FALSE) {
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

      d <- as.integer(as.character(nimbleTable$V2))
      if (any(is.na(d))){
        stop(paste0('Non-integer count values found, were: ', paste0(head(unique(df$nimbleTable[is.na(d)])), collapse = ',')))
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

      origRows <- nrow(nimbleTable)
      nimbleTable <- nimbleTable[nimbleTable$V3 %in% colnames(seuratObj),]
      print(paste0('Original rows: ', origRows, '. After intersecting with seurat object: ', nrow(nimbleTable)))

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
        nimbleTable <- .UpdateLineages(nimbleTable, libraryId = genomeId)
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

  df <- df[c('V1', 'V2', 'V3')]

  print(paste0('Distinct features after re-grouping: ', length(unique(df$V1))))

  return(df)
}

.FindLibraryByName <- function(libraryName) {
  dat <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="reference_libraries",
    colSelect="rowid",
    colSort="name",
    colFilter=makeFilter(c("name", "EQUAL", libraryName)),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(dat) == 0) {
    stop(paste0('Unable to find library: ', libraryName))
  }

  if (nrow(dat) > 1) {
    stop(paste0('More than one library matched: ', libraryName))
  }

  return(dat$rowid[1])
}

#' @title PerformDefaultNimbleAppend
#' @description This is designed to wrap a series of nimble download commands into one, along with some domain-specific logic for each type of data
#' @param seuratObj A Seurat object.
#' @param isotypeFilterThreshold When calculating isotype, any isotype representing below this fraction of reads in the cell is discarded. If this value is 0.1, then a cell with 5 percent of isotype reads for IHGM and 95 percent IGHA would be labeled IGHA.
#' @param maxLibrarySizeRatio Passed directly to nimbleR::AppendNimbleCounts()
#' @param assayForLibrarySize Passed directly to nimbleR::AppendNimbleCounts()
#' @param maxAmbiguityAllowedForKIR Passed to maxAmbiguityAllowed in nimbleR::AppendNimbleCounts() for KIR data specifically
#' @param appendMHC If true, MHC data will be appended
#' @param appendKIR = If true, KIR data will be appended
#' @param appendNKG = If true, NKG2 data will be appended
#' @param appendIG = If true, Ig data will be appended
#' @param appendViral = If true, viral data will be appended
#' @return A modified Seurat object.
#' @export
PerformDefaultNimbleAppend <- function(seuratObj, isotypeFilterThreshold = 0.1, maxLibrarySizeRatio = 100, assayForLibrarySize = 'RNA', maxAmbiguityAllowedForKIR = 2, appendMHC = TRUE, appendKIR = TRUE, appendNKG = TRUE, appendIG = TRUE, appendViral = TRUE) {
  # MHC:
  if (appendMHC) {
    seuratObj <- DownloadAndAppendNimble(seuratObj,
                                       allowableGenomes = .FindLibraryByName('Rhesus Macaque MHC V2'),
                                       targetAssayName = 'MHC',
                                       assayForLibrarySize = assayForLibrarySize,
                                       normalizeData = TRUE,
                                       maxLibrarySizeRatio = maxLibrarySizeRatio,
                                       replaceExistingAssayData = TRUE
    )
    seuratObj <- .GroupMhcData(seuratObj, targetAssay = 'MHC_Grouped')
  }

  # KIR:
  if (appendKIR) {
    seuratObj <- DownloadAndAppendNimble(seuratObj,
                                       allowableGenomes = .FindLibraryByName('Rhesus_KIR'),
                                       targetAssayName = 'KIR',
                                       assayForLibrarySize = assayForLibrarySize,
                                       normalizeData = TRUE,
                                       maxAmbiguityAllowed = maxAmbiguityAllowedForKIR,
                                       maxLibrarySizeRatio = maxLibrarySizeRatio,
                                       replaceExistingAssayData = TRUE,
                                       featureRenameList = NULL
    )
    seuratObj <- .GroupKirData(seuratObj)
  }

  # NKG:
  if (appendNKG) {
    seuratObj <- DownloadAndAppendNimble(seuratObj,
                                       allowableGenomes = .FindLibraryByName('RhesusSupplementalFeatures'),
                                       targetAssayName = 'Nimble',
                                       assayForLibrarySize = assayForLibrarySize,
                                       normalizeData = TRUE,
                                       maxLibrarySizeRatio = maxLibrarySizeRatio,
                                       replaceExistingAssayData = TRUE,
                                       featureRenameList = list(
                                         'NKG2C-KLRC2,NKG2E-KLRC3' = 'NKG2C/E'
                                       )
    )

    seuratObj <- .GroupNkgData(seuratObj)
    seuratObj$NKG_Status <- .IterativeFeatureFiltering(seuratObj, features = c("NKG2A", "NKG2C/E",  "NKG2D"), threshold = 0, maxAllowedClasses = 1, assayName = 'NKG')
    print(sort(table(seuratObj$NKG_Status)))
    print(DimPlot(seuratObj, group.by = 'NKG_Status'))
  }

  # Ig
  if (appendIG) {
    seuratObj <- DownloadAndAppendNimble(seuratObj,
                                       allowableGenomes = .FindLibraryByName('Rhesus_Ig'),
                                       targetAssayName = 'IG',
                                       assayForLibrarySize = assayForLibrarySize,
                                       normalizeData = FALSE,
                                       maxLibrarySizeRatio = maxLibrarySizeRatio,
                                       replaceExistingAssayData = TRUE,
                                       featureRenameList = NULL
    )

    seuratObj <- .MergeNimbleAndRnaIg(seuratObj, assayForLibrarySize = assayForLibrarySize)
    seuratObj <- CalculateIsotype(seuratObj, assayName = 'IG', isotypeFilterThreshold = isotypeFilterThreshold)
  }

  # Viruses:
  if (appendViral) {
    seuratObj <- DownloadAndAppendNimble(seuratObj,
                                       allowableGenomes = .FindLibraryByName('Viral_Genomes'),
                                       targetAssayName = 'Virus',
                                       assayForLibrarySize = assayForLibrarySize,
                                       normalizeData = TRUE,
                                       maxLibrarySizeRatio = NULL,
                                       replaceExistingAssayData = TRUE,
                                       featureRenameList = NULL
    )
  }

  return(seuratObj)
}

.MergeNimbleAndRnaIg <- function(seuratObj, igAssay = 'IG', assayForLibrarySize = 'RNA'){
  ad <- GetAssayData(seuratObj, assay = igAssay, layer = 'counts')
  ighgMap <- c(
    "IGHG1-like" = "LOC708891",
    "IGHG2-like" = "LOC114679691",
    "IGHG4-like" = c("LOC114679690", "LOC710905")
  )

  for (feat in names(ighgMap)) {
    dat <- suppressWarnings(Seurat::FetchData(seuratObj, vars = ighgMap[[feat]]))
    if (ncol(dat) != length(ighgMap[[feat]])) {
      stop(paste0('Did not find all features: ', paste0(ighgMap[[feat]], collapse = ',')))
    }

    dat <- rowSums(dat)
    if (any(names(dat) != colnames(ad))) {
      stop('Cell barcodes not equal!')
    }

    if (feat %in% rownames(ad)) {
      ad[feat] <- dat
    } else {
      ad <- rbind(ad, dat)
      rownames(ad)[length(rownames(ad))] <- feat
    }
  }

  newAssay <- Seurat::CreateAssayObject(counts = ad)
  seuratObj[[igAssay]] <- NULL
  seuratObj[[igAssay]] <- newAssay
  seuratObj <- nimbleR::LogNormalizeUsingAlternateAssay(seuratObj, assay = igAssay, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = NULL)

  return(seuratObj)
}

.GroupNkgData <- function(seuratObj, targetAssay = 'NKG', sourceAssay = 'Nimble', assayForLibrarySize = 'RNA') {
  groupedData <- .RegroupCountMatrix(Seurat::GetAssayData(seuratObj, assay = sourceAssay, layer = 'counts'), featureTransform = function(x){
    return(dplyr::case_when(
      x == 'NKG2A-KLRC1' ~ 'NKG2A',
      x == 'NKG2C-KLRC2' ~ 'NKG2C/E',
      x == 'NKG2E-KLRC3' ~ 'NKG2C/E',
      x == 'NKG2C-KLRC2,NKG2E-KLRC3' ~ 'NKG2C/E',
      .default = x
    ))
  })

  groupedData <- groupedData[c('NKG2A', 'NKG2C/E'),]
  dat <- Matrix::t(Seurat::as.sparse(suppressWarnings(Seurat::FetchData(seuratObj, vars = 'NKG2D'))))
  groupedData <- rbind(groupedData, dat)
  seuratObj[[targetAssay]] <- Seurat::CreateAssayObject(counts = groupedData)
  seuratObj <- nimbleR::LogNormalizeUsingAlternateAssay(seuratObj, assay = targetAssay, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = NULL)

  for (feat in rownames(seuratObj@assays[[targetAssay]])){
    print(FeaturePlot(seuratObj, features = paste0(seuratObj@assays[[targetAssay]]@key, feat)))
  }

  return(seuratObj)
}

.GroupMhcData <- function(seuratObj, targetAssay, sourceAssay = 'MHC', prefix = 'Mamu-', assayForLibrarySize = 'RNA') {
  seuratObj <- .AssignLocusToMhcFeatures(seuratObj, sourceAssayName = sourceAssay)

  ad <- Seurat::GetAssay(seuratObj, assay = sourceAssay)
  groupedMHC <- .RegroupCountMatrix(Seurat::GetAssayData(seuratObj, assay = sourceAssay, layer = 'counts'), featureTransform = function(x){
    return(ad@meta.features$locus[rownames(ad) == x])
  })
  rownames(groupedMHC) <- paste0(prefix, rownames(groupedMHC))

  seuratObj[[targetAssay]] <- Seurat::CreateAssayObject(counts = groupedMHC)
  seuratObj <- nimbleR::LogNormalizeUsingAlternateAssay(seuratObj, assay = targetAssay, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = NULL)

  for (feat in rownames(seuratObj@assays[[targetAssay]])){
    print(FeaturePlot(seuratObj, features = feat))
  }

  return(seuratObj)
}

.GroupKirData <- function(seuratObj, targetAssay = 'KIR_Grouped', sourceAssay = 'KIR', assayForLibrarySize = 'RNA') {
  groupedMat <- .RegroupCountMatrix(Seurat::GetAssayData(seuratObj, assay = sourceAssay, layer = 'counts'), featureTransform = function(x){
    x <- unlist(strsplit(x, split = ','))
    feats <- sapply(x, function(y) {
      return(dplyr::case_when(
        grepl(y, pattern = 'KIR[0-9]DL') ~ 'Inhibitory-KIR',
        grepl(y, pattern = 'KIR[0-9]DS') ~ 'Activating-KIR',
        .default = y
      ))
    })

    feats <- sort(unique(feats))
    if (length(feats) > 1) {
      return('Multiple')
    }

    return(paste0(feats, collapse = ','))
  })

  seuratObj[[targetAssay]] <- Seurat::CreateAssayObject(counts = groupedMat)
  seuratObj <- nimbleR::LogNormalizeUsingAlternateAssay(seuratObj, assay = targetAssay, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = NULL)

  for (feat in rownames(seuratObj@assays[[targetAssay]])){
    print(FeaturePlot(seuratObj, features = paste0(seuratObj@assays[[targetAssay]]@key, feat)))
  }

  seuratObj$KIR_Status <- .IterativeFeatureFiltering(seuratObj, assayName = 'KIR_Grouped', features = rownames(seuratObj@assays$KIR_Grouped), threshold = 0, maxAllowedClasses = 1)

  print(sort(table(seuratObj$KIR_Status)))
  print(DimPlot(seuratObj, group.by = 'KIR_Status'))

  return(seuratObj)
}

#' @title CalculateIsotype
#' @description Uses nimble data to score isotype of each cell
#' @param seuratObj A Seurat object.
#' @param assayName The name of the source assay
#' @param isotypeFilterThreshold When calculating isotype, any isotype representing below this fraction of reads in the cell is discarded. If this value is 0.1, then a cell with 5 percent of isotype reads for IHGM and 95 percent IGHA would be labeled IGHA.
#' @return A modified Seurat object.
#' @export
CalculateIsotype <- function(seuratObj, assayName = 'IG', isotypeFilterThreshold = 0.1) {
  seuratObj$Isotype <- .IterativeFeatureFiltering(seuratObj, assayName = assayName, features = c("IGHA", "IGHG1-like", "IGHG2-like", "IGHG4-like", "IGHM", "IGHD"), maxAllowedClasses = 1, threshold = isotypeFilterThreshold, featureTransform = function(x){
    return(substr(x = x, start = 1, stop = 4))
  })
  print(sort(table(seuratObj$Isotype)))
  print(DimPlot(seuratObj, group.by = 'Isotype'))
  
  seuratObj$IsotypeDetailed <- .IterativeFeatureFiltering(seuratObj, assayName = assayName, features = c("IGHA", "IGHG1-like", "IGHG2-like", "IGHG4-like", "IGHM", "IGHD"), maxAllowedClasses = 2)
  print(sort(table(seuratObj$IsotypeDetailed)))
  print(DimPlot(seuratObj, group.by = 'IsotypeDetailed'))
  
  seuratObj$ClassSwitchedStatus <- .IterativeFeatureFiltering(seuratObj, assayName = assayName, features = c("IGHA", "IGHG1-like", "IGHG2-like", "IGHG4-like", "IGHM", "IGHD"), maxAllowedClasses = 1, threshold = isotypeFilterThreshold, featureTransform = function(x){
    return(dplyr::case_when(
      grepl(x, pattern = 'IGHA|IGHG') ~ 'Class-switched',
      grepl(x, pattern = 'IGHM|IGHD') ~ 'Not Class-switched',
      .default = NA
    ))
  })
  print(sort(table(seuratObj$ClassSwitchedStatus)))
  print(DimPlot(seuratObj, group.by = 'ClassSwitchedStatus'))
  
  return(seuratObj)
}

.RegroupCountMatrix <- function(mat, featureTransform) {
  if (is.null(featureTransform)) {
    return(mat)
  }
  
  updatedRows <- sapply(rownames(mat), featureTransform)
  if (any(is.na(updatedRows))) {
    stop('NAs were produced by featureTransform')
  }
  
  if (any(rownames(mat) != updatedRows)) {
    newMat <- NULL
    for (val in unique(updatedRows)) {
      toSelect <- names(updatedRows)[updatedRows == val]
      m <- matrix(data = Matrix::colSums(mat[toSelect,,drop = FALSE]), nrow = 1)
      rownames(m) <- val
      newMat <- rbind(newMat, Seurat::as.sparse(m))
    }

    colnames(newMat) <- colnames(mat)
    
    return(newMat)
  } else {
    return(mat)
  }
}

.IterativeFeatureFiltering <- function(seuratObj, assayName, features, threshold = 0.05, maxAllowedClasses = 3, featureTransform = NULL) {
  countsMatrix <- Matrix::t(Seurat::FetchData(seuratObj, vars = paste0(seuratObj[[assayName]]@key, features)))
  rownames(countsMatrix) <- gsub(rownames(countsMatrix), pattern = seuratObj[[assayName]]@key, replacement = '')

  cellLabels <- rep(NA, ncol(countsMatrix))
  names(cellLabels) <- colnames(countsMatrix)
  
  countsMatrix <- .RegroupCountMatrix(countsMatrix, featureTransform = featureTransform)

  for (i in 1:ncol(countsMatrix)) {
    v <- countsMatrix[, i]
    repeat {
      total_counts <- sum(v)
      if (total_counts == 0)
        break

      threshold_counts <- total_counts * threshold
      features_to_drop <- which(v < threshold_counts & v > 0)

      if (length(features_to_drop) == 0)
        break

      v[features_to_drop] <- 0
    }

    feats <- sort(rownames(countsMatrix)[v > 0])
    if (length(feats) == 0) {
      next
    }
    
    if (length(feats) > maxAllowedClasses) {
      feats <- 'Multiple'
    }

    cellLabels[i] <- paste0(feats, collapse = ',')
  }

  return(cellLabels)
}
