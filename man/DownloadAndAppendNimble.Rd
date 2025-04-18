% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NimbleAPI.R
\name{DownloadAndAppendNimble}
\alias{DownloadAndAppendNimble}
\title{DownloadAndAppendNimble}
\usage{
DownloadAndAppendNimble(
  seuratObj,
  targetAssayName,
  outPath = tempdir(),
  enforceUniqueFeatureNames = TRUE,
  allowableGenomes = NULL,
  ensureSamplesShareAllGenomes = TRUE,
  maxAmbiguityAllowed = 1,
  reuseExistingDownloads = FALSE,
  performDietSeurat = FALSE,
  normalizeData = TRUE,
  assayForLibrarySize = "RNA",
  maxLibrarySizeRatio = 0.05,
  queryDatabaseForLineageUpdates = FALSE,
  replaceExistingAssayData = TRUE,
  featureRenameList = NULL
)
}
\arguments{
\item{seuratObj}{A Seurat object.}

\item{targetAssayName}{The target assay. If this assay exists, features will be appended (and an error thrown if there are duplicates). Otherwise a new assay will be created.}

\item{outPath}{The path to which nimble files will be downloaded and saved}

\item{enforceUniqueFeatureNames}{Whether or not to fail if we discover multiple nimble files with the same feature name}

\item{allowableGenomes}{An optional vector of genomeIds. If provided, nimble results from these genomes will be appended. If any dataset lacks a nimble file for a genome, it will fail}

\item{ensureSamplesShareAllGenomes}{If true, the function will fail unless all samples have data from the same set of genomes}

\item{maxAmbiguityAllowed}{If provided, any features representing more than ths value will be discarded. For example, 'Feat1,Feat2,Feat3' represents 3 features. maxAmbiguityAllowed=1 results in removal of all ambiguous features.}

\item{reuseExistingDownloads}{If true, any pre-existing downloaded nimble TSVs will be re-used}

\item{performDietSeurat}{If true, DietSeurat will be run, which removes existing reductions. This may or may not be required based on your usage.}

\item{normalizeData}{Passed directly to nimbleR::AppendNimbleCounts()}

\item{assayForLibrarySize}{Passed directly to nimbleR::AppendNimbleCounts()}

\item{maxLibrarySizeRatio}{Passed directly to nimbleR::AppendNimbleCounts()}

\item{queryDatabaseForLineageUpdates}{If true, after downloading the raw nimble output, the code will query any feature not ending with 'g' against the database and replace that name with the current value of lineage.}

\item{replaceExistingAssayData}{If true, any existing data in the targetAssay will be deleted}

\item{featureRenameList}{An optional named list in the format <OLD_NAME> = <NEW_NAME>. If any <OLD_NAME> are present, the will be renamed to <NEW_NAME>. The intention of this is to recover specific ambiguous classes.}
}
\value{
A modified Seurat object.
}
\description{
This read a given seurat object and download/append all associated files to the object
}
