% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCell.R
\name{QueryAndApplyMetadataUsingCDNA}
\alias{QueryAndApplyMetadataUsingCDNA}
\title{QueryAndApplyMetadataUsingCDNA}
\usage{
QueryAndApplyMetadataUsingCDNA(
  seuratObj,
  fieldSelect = c("rowid", "sortid/population", "sortid/sampleid/subjectid",
    "sortid/sampleid/sampledate", "sortid/sampleid/stim", "sortid/sampleid/tissue",
    "plateid", "workbook/workbookid"),
  fieldNames = c("cDNA_ID", "Population", "SubjectId", "SampleDate", "Stim", "Tissue",
    "PlateId", "WorkbookId"),
  overwriteExisting = TRUE
)
}
\arguments{
\item{seuratObj}{A Seurat object.}

\item{fieldSelect}{The set of fields to query}

\item{fieldNames}{The labels to use for the fields}

\item{overwriteExisting}{If true and the output exists, it will be overwritten}
}
\value{
A modified Seurat object.
}
\description{
This will read the cDNA_ID column from the seurat object, query the LK server and apply the appropriate metadata from the cDNA table
}
