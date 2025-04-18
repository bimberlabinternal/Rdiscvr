% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCR.R
\name{DownloadAndAppendTcrClonotypes}
\alias{DownloadAndAppendTcrClonotypes}
\title{DownloadAndAppendTcrClonotypes}
\usage{
DownloadAndAppendTcrClonotypes(
  seuratObj,
  outPath = tempdir(),
  dropExisting = T,
  overwriteTcrTable = F,
  allowMissing = FALSE,
  dropConflictingVJSegments = TRUE
)
}
\arguments{
\item{seuratObj}{A Seurat object}

\item{outPath}{The output filepath}

\item{dropExisting}{If true, any existing clonotype data will be replaced}

\item{overwriteTcrTable}{If true, any existing table(s) of TCR clones will be overwritten and re-downloaded}

\item{allowMissing}{If true, samples missing data will be skipped. Otherwise, the function will fail.}

\item{dropConflictingVJSegments}{If true, any TRB rows with a TRA/D V/J segments will as dropped, as are TRA rows with TRB/G segments}
}
\value{
A modified Seurat object.
}
\description{
Download And Append TCR Clonotypes data from Prime-Seq
}
