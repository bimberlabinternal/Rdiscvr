% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MHC.R
\name{PerformMhcDimRedux}
\alias{PerformMhcDimRedux}
\title{PerformMhcDimRedux}
\usage{
PerformMhcDimRedux(
  seuratObj,
  sourceAssay = "MHC",
  groupField = "SubjectId",
  resolution = 0.05,
  nVariableFeatures = 200,
  dimsForUmap = 1:10
)
}
\arguments{
\item{seuratObj}{A Seurat object.}

\item{sourceAssay}{The assay holding MHC data}

\item{groupField}{An optional field holding the Subject/Sample variable. This is used for plotting only.}

\item{resolution}{The cluster resolution}

\item{nVariableFeatures}{Passed directly to Seurat::FindVariableFeatures}

\item{dimsForUmap}{Passed to the dims argument of Seurat::RunUMAP}
}
\description{
Performs PCA/UMAP on MHC typing data, often for the purpose of
}
