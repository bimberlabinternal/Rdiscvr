% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NimbleAPI.R
\name{PerformDefaultNimbleAppend}
\alias{PerformDefaultNimbleAppend}
\title{PerformDefaultNimbleAppend}
\usage{
PerformDefaultNimbleAppend(
  seuratObj,
  isotypeFilterThreshold = 0.1,
  maxLibrarySizeRatio = 100,
  assayForLibrarySize = "RNA",
  maxAmbiguityAllowedForKIR = 2
)
}
\arguments{
\item{seuratObj}{A Seurat object.}

\item{isotypeFilterThreshold}{When calculating isotype, any isotype representing below this fraction of reads in the cell is discarded. If this value is 0.1, then a cell with 5 percent of isotype reads for IHGM and 95 percent IGHA would be labeled IGHA.}

\item{maxLibrarySizeRatio}{Passed directly to AppendNimbleCounts()}

\item{assayForLibrarySize}{Passed directly to AppendNimbleCounts()}

\item{maxAmbiguityAllowedForKIR}{Passed to maxAmbiguityAllowed in AppendNimbleCounts() for KIR data specifically}
}
\value{
A modified Seurat object.
}
\description{
This is designed to wrap a series of nimble download commands into one, along with some domain-specific logic for each type of data
}
