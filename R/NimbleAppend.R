#' @import Seurat
#' @importFrom tidyr pivot_wider

#' @title AppendNimbleCounts
#' @description Reads a given seurat object and a nimble file, and appends the nimble data to the object.
#'
#' @param seuratObject A Seurat object.
#' @param nimbleFile A nimble file, which is a TSV of feature counts created by nimble
#' @param dropAmbiguousFeatures If true, any ambiguous features (defined as containing a comma) will be discarded
#' @param targetAssayName The target assay. If this assay exists, features will be appended (and an error thrown if there are duplicates). Otherwise a new assay will be created.
#' @param renameConflictingFeatures If true, when appending to an existing assay, any conflicting feature names will be renamed, appending the value of duplicateFeatureSuffix
#' @param duplicateFeatureSuffix If renameConflictingFeatures is true, this string will be appended to duplicated feature names
#' @param normalizeData If true, data will be normalized after appending/creating the assay. This will default to CellMembrane::LogNormalizeUsingAlternateAssay; however, if assayForLibrarySize equals targetAssayName then Seurat::NormalizeData is used.
#' @param performDietSeurat If true, DietSeurat will be run, which removes existing reductions. This may or may not be required based on your usage, but the default is true out of caution.
#' @param assayForLibrarySize If normalizeData is true, then this is the assay used for librarySize when normalizing. If assayForLibrarySize equals targetAssayName, Seurat::NormalizeData is used.
#' @param maxLibrarySizeRatio If normalizeData is true, then this is passed to CellMembrane::LogNormalizeUsingAlternateAssay
#' @return A modified Seurat object.
#' @export
AppendNimbleCounts <- function(seuratObject, nimbleFile, targetAssayName, dropAmbiguousFeatures = TRUE, renameConflictingFeatures = TRUE, duplicateFeatureSuffix = ".Nimble", normalizeData = TRUE, performDietSeurat = TRUE, assayForLibrarySize = 'RNA', maxLibrarySizeRatio = 0.05) {
  if (!file.exists(nimbleFile)) {
    stop(paste0("Nimble file not found: ", nimbleFile))
  }
  
  # Read file and construct df
  df <- NULL
  tryCatch({
    df <- read.table(nimbleFile, sep="\t", header=FALSE)
  }, error = function(e){
    if (conditionMessage(e) != 'no lines available in input') {
      stop(e)
    } else {
      print(paste0('No lines in nimble file: ', nimbleFile))
      df <- data.frame(V1 = character(), V2 = character(), V3 = character())
    }
  })

  # Indicates no data in TSV
  if (all(is.null(df))){
    return(seuratObject)
  }

  if (sum(df$V1 == "") > 0) {
    stop("The nimble data contains blank feature names. This should not occur.")
  }

  if (sum(grepl(df$V1, pattern = "^,")) > 0) {
    stop("The nimble data contains features with leading commas. This should not occur.")
  }

  if (sum(grepl(df$V1, pattern = ",$")) > 0) {
    stop("The nimble data contains features with trailing commas. This should not occur.")
  }

  #Remove ambiguous features
  ambigFeatRows <- grepl(",", df$V1)
  if (sum(ambigFeatRows) > 0) {
    if (dropAmbiguousFeatures) {
      print(paste0('Dropping ', sum(ambigFeatRows), ' rows with ambiguous features. (', sum(ambigFeatRows),' of ', nrow(df), ')'))
      x <- df$V1[ambigFeatRows]

      # For the purposes of reporting only, collapse highly ambiguous results
      totalHitsByRow <- sapply(x, function(y){
        return(length(unlist(strsplit(y, split = ','))))
      })
      x[totalHitsByRow > 3] <- 'ManyHits'

      x <- sort(table(x), decreasing = T)
      x <- data.frame(Feature = names(x), Total = as.numeric(unname(x)))
      print(x)
      df <- df[!ambigFeatRows, , drop = F]
    }
  }

  tryCatch({
    df <- tidyr::pivot_wider(df, names_from=V3, values_from=V2, values_fill=0)
  }, error = function(e){
    write.table(df, file = 'debug.nimble.txt', sep = '\t', quote = F, row.names = F)

    print('Error pivoting input data, results saved to: debug.nimble.txt')
    print(conditionMessage(e))
    traceback()
    e$message <- paste0('Error pivoting nimble data. target assay: ', targetAssayName)
    stop(e)
  })

  appendToExistingAssay <- targetAssayName %in% names(seuratObject@assays)

  # Remove barcodes from nimble that aren't in seurat
  seuratBarcodes <- colnames(seuratObject@assays[[Seurat::DefaultAssay(seuratObject)]])
  barcodeDiff <- colnames(df) %in% seuratBarcodes
  barcodeDiff[1] <- TRUE # retain feature name
  numColsToDrop <- length(barcodeDiff) - sum(barcodeDiff)
  if (numColsToDrop > 0) {
    print(paste0('Dropping ', numColsToDrop, ' cell barcodes not in the seurat object (out of ', (ncol(df)-1), ')'))
  }
  df <- df[barcodeDiff]
  
  # Fill zeroed barcodes that are in seurat but not in nimble
  zeroedBarcodes <- setdiff(seuratBarcodes, colnames(df)[-1])
  print(paste0('Total cells lacking nimble data: ', length(zeroedBarcodes), ' of ', length(seuratBarcodes), ' cells'))
  for (barcode in zeroedBarcodes) {
    df[barcode] <- 0
  }
  
  # Cast nimble df to matrix
  featureNames <- df$V1
  df <- subset(df, select=-(V1))
  m <- Reduce(methods::cbind2, lapply(df, Matrix::Matrix, sparse = TRUE))
  dimnames(m) <- list(featureNames, colnames(df))
  if (is.null(colnames(m))) {
    stop(paste0('Error: no column names in nimble count matrix, size: ', paste0(dim(m), collapse = ' by ')))
  }

  m <- m[,seuratBarcodes, drop=FALSE] # Ensure column order matches
  if (appendToExistingAssay && ncol(m) != ncol(seuratObject@assays[[targetAssayName]])) {
    stop(paste0('Error parsing nimble data, ncol not equal after subset, was ', ncol(m)))
  }

  if (is.null(colnames(m))) {
    stop(paste0('Error: no column names in matrix after subset, size: ', paste0(dim(m), collapse = ' by ')))
  }

  if (appendToExistingAssay) {
    if (any(rownames(m) %in% rownames(seuratObject@assays[[targetAssayName]]))) {
      conflicting <- rownames(m)[rownames(m) %in% rownames(seuratObject@assays[[targetAssayName]])]

      if (renameConflictingFeatures) {
        print(paste0('The following nimble features have conflicts in the existing assay and will be renamed: ', paste0(conflicting, collapse = ',')))
        newNames <- rownames(m)
        names(newNames) <- newNames
        newNames[conflicting] <- paste0(conflicting, duplicateFeatureSuffix)
        newNames <- unname(newNames)
        rownames(m) <- newNames
      } else {
        stop(paste0('The following nimble features conflict with features in the seuratObj: ', paste0(conflicting, collapse = ',')))
      }
    }

    # NOTE: always perform DietSeurat() to drop products:
    if (performDietSeurat) {
      seuratObject <- Seurat::DietSeurat(seuratObject)
    }

    # Append nimble matrix to seurat count matrix
    existingBarcodes <- colnames(Seurat::GetAssayData(seuratObject, assay = targetAssayName, slot = 'counts'))
    if (sum(colnames(m) != existingBarcodes) > 0) {
      stop('cellbarcodes do not match on matrices')
    }

    # If feature source exists, retain it. Otherwise assume these are from cellranger:
    if ('FeatureSource' %in% names(seuratObject@assays[[targetAssayName]]@meta.features)) {
      fs <- seuratObject@assays[[targetAssayName]]@meta.features$FeatureSource
    } else {
      fs <- rep('CellRanger', nrow(seuratObject@assays[[targetAssayName]]))
    }

    fs <- c(fs, rep('Nimble', nrow(m)))

    seuratObject[[targetAssayName]] <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(rbind(Seurat::GetAssayData(seuratObject, assay = targetAssayName, slot = 'counts'), m)))

    names(fs) <- rownames(seuratObject@assays[[targetAssayName]])
    seuratObject@assays[[targetAssayName]] <- Seurat::AddMetaData(seuratObject@assays[[targetAssayName]], metadata = fs, col.name = 'FeatureSource')

    if (sum(colnames(Seurat::GetAssayData(seuratObject, assay = targetAssayName, slot = 'counts')) != existingBarcodes) > 0) {
      stop('cellbarcodes do not match on matrices after assay replacement')
    }
  } else {
    # Add nimble as separate assay
    seuratObject[[targetAssayName]] <- Seurat::CreateAssayObject(counts = m, min.features = 0, min.cells = 0)

    fs <- rep('Nimble', nrow(seuratObject@assays[[targetAssayName]]))
    names(fs) <- rownames(seuratObject@assays[[targetAssayName]])
    seuratObject@assays[[targetAssayName]] <- Seurat::AddMetaData(seuratObject@assays[[targetAssayName]], metadata = fs, col.name = 'FeatureSource')
  }

  if (normalizeData) {
    if (targetAssayName == assayForLibrarySize) {
      seuratObject <- Seurat::NormalizeData(seuratObject, assay = targetAssayName, verbose = FALSE)
    } else {
      seuratObject <- CellMembrane::LogNormalizeUsingAlternateAssay(seuratObject, assay = targetAssayName, assayForLibrarySize = assayForLibrarySize, maxLibrarySizeRatio = maxLibrarySizeRatio)
    }
  }
  
  return(seuratObject)
}

.AssignLocusToMhcFeatures <- function(seuratObj, sourceAssayName = 'MHC', featurePrefix = 'Mamu-', delimiter = '*', ambiguousFeatureDelim = ',', stripNumbersFromLocus = TRUE) {
  seuratObj[[sourceAssayName]]@meta.features$locus <- NA
  for (featName in rownames(seuratObj[[sourceAssayName]])) {
    feats <- unlist(strsplit(x = featName, split = ambiguousFeatureDelim))
    loci <- c()
    for (feat in feats){
      qs <- paste0('^', featurePrefix)
      if (!grepl(x = feat, pattern = qs)) {
        warning(paste0('Feature lacks prefix: ', feat))
      }

      locus <- gsub(x = feat, pattern = featurePrefix, replacement = '')
      if (!grepl(x = locus, pattern = delimiter, fixed = TRUE)) {
        warning(paste0('Feature lacks delimiter: ', feat))
      }

      locus <- unlist(strsplit(x = locus, split = delimiter, fixed = TRUE))[1]
      if (stripNumbersFromLocus) {
        locus <- gsub(x = locus, pattern = '[0-9]+.*$', replacement = '')
      }

      loci <- unique(c(loci, locus))
    }

    if (length(loci) > 1) {
      warning(paste0('Feature matched multiple loci: ', featName, ', ', paste0(loci, collapse = ',')))
    }

    seuratObj[[sourceAssayName]]@meta.features$locus[rownames(seuratObj[[sourceAssayName]]) == feat] <- paste0(loci, collapse = ',')
  }

  return(seuratObj)
}

#' @title PerformMhcNormalization
#' @description This is a fairly specific normalization step for MHC data. It will divide the raw counts for each feature by the sum of counts in that cell from that locus (e.g., MHC-A, MHC-B, MHC-E, MHC-I, DPA, DPB)
#'
#' @param seuratObj A Seurat object
#' @param sourceAssayName The assay to normalize
#' @param featurePrefix This prefix is stripped from the start of all feature names
#' @param delimiter Used to split the locus from allele designation
#' @param ambiguousFeatureDelim This character is used to split feature names in the case of ambiguous features. If a feature is ambiguous, the locus is assigned as the unique loci of the feature set.
#' @param perCell If true, the feature counts are scaled based on the library size of features from that locus in that cell. If false, it is scaled based on the library size of features in that locus from all cells matching cellGroupingVariable
#' @param cellGroupingVariable If perCell is FALSE, the library size is calculated by taking the sum of features from that locus across all cells where this metadata variable matches the current cell
#' @param stripNumbersFromLocus If true, numeric values will be stripped from all locus strings 
#' @return A modified Seurat object.
#' @import ggplot2
#' @export
PerformMhcNormalization <- function(seuratObj, sourceAssayName = 'MHC', featurePrefix = 'Mamu-', delimiter = '*', ambiguousFeatureDelim = ',', perCell = TRUE, cellGroupingVariable = 'DatasetId', stripNumbersFromLocus = TRUE) {
  seuratObj <- .AssignLocusToMhcFeatures(seuratObj, sourceAssayName = sourceAssayName, featurePrefix = featurePrefix, delimiter = delimiter, ambiguousFeatureDelim = ambiguousFeatureDelim, stripNumbersFromLocus = stripNumbersFromLocus)

  dat <- Seurat::GetAssayData(seuratObj, assay = sourceAssayName, slot = 'counts')
  margin <- 2

  if (!perCell) {
    if (!cellGroupingVariable %in% names(seuratObj@meta.data)) {
      stop(paste0('The cellGroupingVariable of ', cellGroupingVariable, ' was not present seuratObj@meta.data'))
    }
    
    if (any(is.na(seuratObj[[cellGroupingVariable]]))) {
      stop(paste0('The cellGroupingVariable cannot have NAs in it: ', cellGroupingVariable))
    }
    
    librarySizeData <- NULL
    for (locus in sort(unique(seuratObj[[sourceAssayName]]@meta.features$locus))) {
      groupNames <- unique(seuratObj@meta.data[[cellGroupingVariable]])
      for (gn in groupNames) {
        scale.factor <- sum(seuratObj@meta.data[[cellGroupingVariable]] == gn)  # the number of cells in the group
        librarySize <- sum(dat[seuratObj[[sourceAssayName]]@meta.features$locus == locus, colnames(seuratObj)[seuratObj@meta.data[[cellGroupingVariable]] == gn], drop = TRUE])
        
        toAdd <- data.frame(locus = locus, groupName = gn, librarySize = librarySize, scale.factor = scale.factor)
        if (all(is.null(librarySizeData))) {
          librarySizeData <- toAdd
        } else {
          librarySizeData <- rbind(librarySizeData, toAdd)
        }
      }
    }
  }

  for (locus in sort(unique(seuratObj[[sourceAssayName]]@meta.features$locus))) {
    print(paste0('Normalizing locus: ', locus))
    librarySizes <- c()
    toNormalize <- dat[seuratObj[[sourceAssayName]]@meta.features$locus == locus,,drop = FALSE]
    ncells <- dim(x = toNormalize)[margin]

    for (i in seq_len(length.out = ncells)) {
      x <- toNormalize[, i]
      if (perCell) {
        librarySize <- sum(x)
        scale.factor <- 1
      } else {
        groupName <- seuratObj@meta.data[[cellGroupingVariable]][i]
        scale.factor <- librarySizeData$scale.factor[librarySizeData$locus == locus & librarySizeData$groupName == groupName]
        librarySize <- librarySizeData$librarySize[librarySizeData$locus == locus & librarySizeData$groupName == groupName]
      }
      
      librarySizes <- c(librarySizes, librarySize)

      if (librarySize == 0) {
        xnorm <- 0
      } else {
        xnorm <- x / librarySize * scale.factor
      }

      toNormalize[, i] <- xnorm
    }
    
    print(paste0('Total features: ', nrow(toNormalize)))
    print(paste0('Mean  library size: ', mean(librarySizes), ', min: ', min(librarySizes), ', max: ', max(librarySizes)))
    
    print(ggplot(data.frame(x = as.numeric(toNormalize)), aes(x = x)) +
      geom_density() +
      egg::theme_presentation(base_size = 12) +
      labs(x = 'Normalized Value', y = 'Density') +
      ggtitle(paste0('Normalized Data: ', locus)))

    dat[seuratObj[[sourceAssayName]]@meta.features$locus == locus] <- toNormalize
  }

  seuratObj <- Seurat::SetAssayData(seuratObj, assay = sourceAssayName, slot = 'data', new.data = toNormalize)

  return(seuratObj)
}
