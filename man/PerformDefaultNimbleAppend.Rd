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
  maxAmbiguityAllowedForKIR = 2,
  appendMHC = TRUE,
  appendKIR = TRUE,
  appendNKG = TRUE,
  appendIG = TRUE,
  appendViral = TRUE
)
}
\arguments{
\item{seuratObj}{A Seurat object.}

\item{isotypeFilterThreshold}{When calculating isotype, any isotype representing below this fraction of reads in the cell is discarded. If this value is 0.1, then a cell with 5 percent of isotype reads for IHGM and 95 percent IGHA would be labeled IGHA.}

\item{maxLibrarySizeRatio}{Passed directly to nimbleR::AppendNimbleCounts()}

\item{assayForLibrarySize}{Passed directly to nimbleR::AppendNimbleCounts()}

\item{maxAmbiguityAllowedForKIR}{Passed to maxAmbiguityAllowed in nimbleR::AppendNimbleCounts() for KIR data specifically}

\item{appendMHC}{If true, MHC data will be appended}

\item{appendKIR}{= If true, KIR data will be appended}

\item{appendNKG}{= If true, NKG2 data will be appended}

\item{appendIG}{= If true, Ig data will be appended}

\item{appendViral}{= If true, viral data will be appended}
}
\value{
A modified Seurat object.
}
\description{
This is designed to wrap a series of nimble download commands into one, along with some domain-specific logic for each type of data
}
