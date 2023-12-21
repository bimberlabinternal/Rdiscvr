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
#' @return A modified Seurat object.
#' @export
AppendNimbleCounts <- function(seuratObject, nimbleFile, targetAssayName, dropAmbiguousFeatures = TRUE, renameConflictingFeatures = TRUE, duplicateFeatureSuffix = ".Nimble", normalizeData = TRUE, performDietSeurat = TRUE, assayForLibrarySize = 'RNA') {
  if (!file.exists(nimbleFile)) {
    stop(paste0("Nimble file not found: ", nimbleFile))
  }
  
  # Read file and construct df
  df <- read.table(nimbleFile, sep="\t", header=FALSE)
  
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
      print(paste0('Dropping ', sum(ambigFeatRows), ' ambiguous features. (', sum(ambigFeatRows),' of ', nrow(df), ')'))
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
    stop('Error preparing nimble data')
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
    existingBarcodes <- colnames(seuratObject@assays[[targetAssayName]]@counts)
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

    seuratObject[[targetAssayName]] <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(rbind(seuratObject@assays[[targetAssayName]]@counts, m)))

    print('adding 1')
    names(fs) <- rownames(seuratObject@assays[[targetAssayName]])
    seuratObject@assays[[targetAssayName]] <- Seurat::AddMetaData(seuratObject@assays[[targetAssayName]], metadata = fs, col.name = 'FeatureSource')

    if (sum(colnames(seuratObject@assays[[targetAssayName]]@counts) != existingBarcodes) > 0) {
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
      seuratObject <- CellMembrane::LogNormalizeUsingAlternateAssay(seuratObject, assay = targetAssayName, assayForLibrarySize = assayForLibrarySize)
    }
  }
  
  return(seuratObject)
}
