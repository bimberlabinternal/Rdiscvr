% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCR.R
\name{CreateMergedTcrClonotypeFile}
\alias{CreateMergedTcrClonotypeFile}
\title{CreateMergedTcrClonotypeFile}
\usage{
CreateMergedTcrClonotypeFile(
  seuratObj,
  outputFile,
  overwriteTcrTable = F,
  downloadPath = tempdir(),
  allowMissing = FALSE,
  cellRangerType = "filtered_contig_annotations.csv",
  dropConflictingVJSegments = TRUE
)
}
\arguments{
\item{seuratObj}{A Seurat object}

\item{outputFile}{The path where the merged CSV will be written}

\item{overwriteTcrTable}{If true, any existing table(s) of TCR clones will be overwritten and re-downloaded}

\item{downloadPath}{The output filepath for per-dataset files}

\item{allowMissing}{If true, samples missing data will be skipped. Otherwise, the function will fail.}

\item{cellRangerType}{The type of cellranger data to download. Either all_contig_annotations.csv or filtered_contig_annotations.csv}

\item{dropConflictingVJSegments}{If true, any TRB rows with a TRA/D V/J segments will as dropped, as are TRA rows with TRB/G segments}
}
\description{
Download TCR Clonotypes for all datasets in a seuratObj, update their cellbarcodes with barcodePrefix, and create a merged table
}
