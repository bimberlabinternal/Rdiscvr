#' @include Utils.R

#' @title PerformMhcDimRedux
#' @description Performs PCA/UMAP on MHC typing data, often for the purpose of
#'
#' @param seuratObj A Seurat object.
#' @param sourceAssay The assay holding MHC data
#' @param groupField An optional field holding the Subject/Sample variable. This is used for plotting only.
#' @param resolution The cluster resolution
#'
#' @export
PerformMhcDimRedux <- function(seuratObj, sourceAssay = 'MHC', groupField = 'SubjectId', resolution = 0.05) {
  origDefaultAssay <- DefaultAssay(seuratObj)
  origDefaultVariableFeatures <- VariableFeatures(seuratObj)

  DefaultAssay(seuratObj) <- sourceAssay

  fn <- paste0("MHC_snn_res.", resolution)
  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE, assay = sourceAssay)
  seuratObj <- Seurat::ScaleData(seuratObj, verbose = FALSE, assay = sourceAssay)
  seuratObj <- FindVariableFeatures(seuratObj, verbose = FALSE, assay = sourceAssay)

  seuratObj <- RunPCA(seuratObj, verbose = FALSE, assay = sourceAssay, reduction.name = 'mhc.pca', reduction.key = 'mhcPCA_')
  seuratObj <- FindNeighbors(seuratObj, verbose = FALSE, assay = sourceAssay, reduction = 'mhc.pca', graph.name = 'MHC_snn')
  seuratObj <- FindClusters(object = seuratObj,
                               resolution = resolution,
                               verbose = FALSE,
                               random.seed = CellMembrane::GetSeed(),
                               method = 'igraph',
                               algorithm = 4,
                               graph.name = 'MHC_snn',
                               cluster.name = fn
  )
  seuratObj <- RunUMAP(seuratObj, dims = 1:10, verbose = FALSE, assay = sourceAssay, reduction = 'mhc.pca', reduction.name = 'mhc.umap', reduction.key = "mhcUMAP_")

  print(DimPlot(seuratObj, group.by = fn, reduction = 'mhc.umap'))

  if (! is.null(groupField)) {
    print(CellMembrane::PlotSeuratVariables(seuratObj, xvar = fn, yvar = groupField, labelDimplot = TRUE, reduction = 'mhc.umap'))
  }

  DefaultAssay(seuratObj) <- origDefaultAssay
  VariableFeatures(seuratObj) <- origDefaultVariableFeatures

  return(seuratObj)
}

#' @title GroupByMhcSimilarity
#' @description Uses MHC typing data to pseudobulk data by the supplied variable and then reports a matrix of distances
#'
#' @param seuratObj A Seurat object.
#' @param sourceAssay The assay holding MHC data
#' @param groupField The field on which to group the data
#'
#' @export
GroupByMhcSimilarity <- function(seuratObj, groupField, sourceAssay = 'MHC', dist.method = 'euclidean') {
  origDefaultAssay <- DefaultAssay(seuratObj)
  origDefaultVariableFeatures <- VariableFeatures(seuratObj)

  seuratObjGrouped <- Seurat::AggregateExpression(seuratObj, group.by = groupField, assays = sourceAssay, verbose = FALSE, return.seurat = TRUE)
  seuratObjGrouped <- PerformMhcNormalization(seuratObjGrouped, sourceAssayName = sourceAssay)
  mhcData <- Seurat::GetAssayData(seuratObjGrouped, assay = sourceAssay, layer = 'data')
  dm <- stats::dist(Matrix::t(mhcData), method = dist.method)

  ComplexHeatmap::Heatmap(as.matrix(dm))

  return(dm)
}

.AssignLocusToMhcFeatures <- function(seuratObj, sourceAssayName = 'MHC', featurePrefix = 'Mamu-', delimiter = '*', ambiguousFeatureDelim = ',', stripNumbersFromLocus = TRUE) {
  slotName <- .GetAssayMetaSlotName(seuratObj[[sourceAssayName]])
  methods::slot(seuratObj[[sourceAssayName]], slotName)$loci <- NA

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

    methods::slot(seuratObj[[sourceAssayName]], slotName)$locus[rownames(seuratObj[[sourceAssayName]]) == feat] <- paste0(loci, collapse = ',')
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
  assayMeta <- .GetAssayMeta(seuratObj[[sourceAssayName]])
  margin <- 2

  if (!perCell) {
    if (!cellGroupingVariable %in% names(seuratObj@meta.data)) {
      stop(paste0('The cellGroupingVariable of ', cellGroupingVariable, ' was not present seuratObj@meta.data'))
    }

    if (any(is.na(seuratObj[[cellGroupingVariable]]))) {
      stop(paste0('The cellGroupingVariable cannot have NAs in it: ', cellGroupingVariable))
    }

    librarySizeData <- NULL
    for (locus in sort(unique(assayMeta$locus))) {
      groupNames <- unique(seuratObj@meta.data[[cellGroupingVariable]])
      for (gn in groupNames) {
        scale.factor <- sum(seuratObj@meta.data[[cellGroupingVariable]] == gn)  # the number of cells in the group
        librarySize <- sum(dat[assayMeta$locus == locus, colnames(seuratObj)[seuratObj@meta.data[[cellGroupingVariable]] == gn], drop = TRUE])

        toAdd <- data.frame(locus = locus, groupName = gn, librarySize = librarySize, scale.factor = scale.factor)
        if (all(is.null(librarySizeData))) {
          librarySizeData <- toAdd
        } else {
          librarySizeData <- rbind(librarySizeData, toAdd)
        }
      }
    }
  }

  for (locus in sort(unique(assayMeta$locus))) {
    print(paste0('Normalizing locus: ', locus))
    librarySizes <- c()
    toNormalize <- dat[assayMeta$locus == locus,,drop = FALSE]
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

    dat[assayMeta$locus == locus] <- toNormalize
  }

  seuratObj <- Seurat::SetAssayData(seuratObj, assay = sourceAssayName, slot = 'data', new.data = dat)

  return(seuratObj)
}