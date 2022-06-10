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
#' @param normalizeData If true, Seurat::NormalizeData will be run after appending/creating the assay
#' @return A modified Seurat object.
#' @export
AppendNimbleCounts <- function(seuratObject, nimbleFile, targetAssayName, dropAmbiguousFeatures = TRUE, renameConflictingFeatures = TRUE, duplicateFeatureSuffix = ".Nimble", normalizeData = TRUE) {
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
      x <- sort(table(df$V1[ambigFeatRows]), decreasing = T)
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
  dimnames(m) <- list(featureNames, seuratBarcodes)
  if (is.null(colnames(m))) {
    stop(paste0('Error: no column names in nimble count matrix, size: ', paste0(dim(m), collapse = ' by ')))
  }

  m <- m[,seuratBarcodes] # Ensure column order matches
  if (appendToExistingAssay && ncol(m) != ncol(seuratObject@assays[[targetAssayName]])) {
    stop(paste0('Error parsing nimble data, ncol not equal after subset, was ', ncol(m)))
  }

  if (is.null(colnames(m))) {
    stop('Error: no column names in matrix after subset')
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
    seuratObject <- DietSeurat(seuratObject)

    # Append nimble matrix to seurat count matrix
    existingBarcodes <- colnames(seuratObject@assays[[targetAssayName]]@counts)
    if (sum(colnames(m) != existingBarcodes) > 0) {
      stop('cellbarcodes do not match on matrices')
    }

    seuratObject[[targetAssayName]] <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(rbind(seuratObject@assays[[targetAssayName]]@counts, m)))

    if (sum(colnames(seuratObject@assays[[targetAssayName]]@counts) != existingBarcodes) > 0) {
      stop('cellbarcodes do not match on matrices after assay replacement')
    }
  } else {
    # Add nimble as separate assay
    seuratObject[[targetAssayName]] <- Seurat::CreateAssayObject(counts = m, min.features = 0, min.cells = 0)
  }

  if (normalizeData) {
    seuratObject <- Seurat::NormalizeData(seuratObject, assay = targetAssayName, verbose = FALSE)
  }
  
  return(seuratObject)
}
