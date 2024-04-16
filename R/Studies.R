
#' @title ApplyPC475Metadata
#' @description Applies standard metadata related to PC475 / BMGF HIV Reservoirs project
#'
#' @param seuratObj A Seurat object.
#' @param errorIfUnknownIdsFound If true, the function will fail if the seurat object contains unknown IDs
#' @param reApplyMetadata If true, QueryAndApplyCdnaMetadata will be re-run
#' @return A modified Seurat object.
#' @export
ApplyPC475Metadata <- function(seuratObj, errorIfUnknownIdsFound = TRUE, reApplyMetadata = TRUE) {
  if (reApplyMetadata) {
    seuratObj <- .ApplyMetadata(seuratObj)
  }

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/587",
    schemaName="lists",
    queryName="PC475",
    colNameOpt="rname",
    colSelect = 'libraryid,label,timepoint,wpi,LibraryId/sortId/sampleId/subjectId,Nx_pVL,SIV_RNA,SIV_DNA,CAVL_SampleName',
  )
  names(metadata) <- c('cDNA_ID', 'TimepointLabel', 'Timepoint', 'WPI', 'SubjectId', 'Nx_pVL', 'SIV_RNA', 'SIV_DNA', 'CAVL_SampleName')

  metadata$TimepointLabel <- naturalsort::naturalfactor(metadata$TimepointLabel)
  metadata$Timepoint <- naturalsort::naturalfactor(metadata$Timepoint)

  metadata2 <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber",
    schemaName="laboratory",
    queryName="project_usage",
    colNameOpt="rname",
    colSelect = 'subjectid,groupname',
    colFilter = makeFilter(c('project', 'EQUALS', 'PC475 / HIV Cure'))
  )
  names(metadata2) <- c('SubjectId', 'PC475_Group')

  metadata <- merge(metadata, metadata2, by = 'SubjectId', all.x = T)
  metadata <- metadata[names(metadata) != 'SubjectId']

  if (errorIfUnknownIdsFound && !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID)) {
    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }
  
  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- 1:nrow(toAdd)
  toAdd <- merge(toAdd, metadata, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))

  return(seuratObj)
}

#' @title ApplyTBMetadata
#' @description Applies standard metadata related to TB projects
#'
#' @param seuratObj A Seurat object.
#' @param errorIfUnknownIdsFound If true, the function will fail if the seurat object contains unknown IDs
#' @param reApplyMetadata If true, QueryAndApplyCdnaMetadata will be re-run
#' @return A modified Seurat object.
#' @export
ApplyTBMetadata <-function(seuratObj, errorIfUnknownIdsFound = TRUE, reApplyMetadata = TRUE) {
  if (reApplyMetadata) {
    seuratObj <- .ApplyMetadata(seuratObj)
  }

  metadata <- .GetTbMetadata()

  if (errorIfUnknownIdsFound && !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID)) {
    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }

  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- 1:nrow(toAdd)
  toAdd <- merge(toAdd, metadata, by = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))

  seuratObj$Tissue <- naturalsort::naturalfactor(seuratObj$Tissue)

  return(seuratObj)
}

.GetTbMetadata <- function() {
  cDNA <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/627",
    schemaName="lists",
    queryName="TB_cDNA_Libraries",
    colSelect="cDNA_ID,SampleType,TimepointLabel,CFU_Homogenate,CFU_Tissue,cDNA_ID/sortid/sampleid/subjectid",
    colNameOpt="rname"
  )
  names(cDNA) <- c('cDNA_ID', 'SampleType', 'TimepointLabel', 'CFU_Homogenate', 'CFU_Tissue', 'SubjectId')

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/627",
    schemaName="lists",
    queryName="TB_Studies",
    colSelect="subjectid,study,vaccine,challenge,ChallengeType,NecropsyDate,PID,PathScore,LungPathScore,ChallengeDate,R_Caudal_Lung_1_Disease",
    colNameOpt="rname"
  )
  names(metadata) <- c('SubjectId', 'TB_Study', 'Vaccine', 'Challenge', 'ChallengeType', 'NecropsyDate', 'PID', 'PathScore', 'LungPathScore', 'ChallengeDate', 'R_Caudal_Lung_1_Disease')

  # CFU_Tissue_Rescaled
  cDNA$CFU_Homogenate[is.na(cDNA$CFU_Homogenate)] <- 0
  cDNA$CFU_Homogenate <- as.numeric(cDNA$CFU_Homogenate)
  cDNA$CFU_Homogenate_Rescaled <- scales::rescale(asinh(cDNA$CFU_Homogenate + 1), to = c(0,1))

  cDNA$CFU_Tissue[is.na(cDNA$CFU_Tissue)] <- 0
  cDNA$CFU_Tissue <- as.numeric(cDNA$CFU_Tissue)
  cDNA$CFU_Tissue_Rescaled <- scales::rescale(asinh(cDNA$CFU_Tissue + 1), to = c(0,1))

  metadata$IsMockChallenged <- metadata$Challenge == 'Mock-challenged'

  #Round to week:
  metadata$Timepoint <- metadata$PID
  metadata$Timepoint[!is.na(metadata$Timepoint)] <- round(metadata$Timepoint[!is.na(metadata$Timepoint)]/7, 0)*7
  metadata$Timepoint[!is.na(metadata$Timepoint)] <- paste0('Day ', metadata$Timepoint[!is.na(metadata$Timepoint)])

  # NOTE: make a simpler timepoint/mock field that lumps Mock-challenged samples into a different timepoint level
  metadata$Timepoint[metadata$IsMockChallenged] <- paste0('Mock / ', metadata$Vaccine[metadata$IsMockChallenged])
  if ('Mock / Unvaccinated' %in% metadata$Timepoint) {
    metadata$Timepoint[metadata$Timepoint == 'Mock / Unvaccinated'] <- 'Mock'
  }

  # And then place the Mocks first:
  metadata$Timepoint <- naturalsort::naturalfactor(metadata$Timepoint)
  mocks <- grep(metadata$Timepoint, pattern = '^Mock', value = TRUE)
  if (length(mocks) > 0) {
    for (l in rev(mocks)) {
      if (l %in% unique(metadata$Timepoint)) {
        metadata$Timepoint <- forcats::fct_relevel(metadata$Timepoint, l, after = 0)
      }
    }
  }

  metadata$Group <- as.character(metadata$Vaccine)
  if ('Mock-challenged' %in% metadata$Challenge) {
    metadata$Group[metadata$Challenge == 'Mock-challenged'] <- paste0(metadata$Vaccine[metadata$Challenge == 'Mock-challenged'], '-Mock')
  }
  metadata$Group <- naturalsort::naturalfactor(metadata$Group)

  # Establish order:
  expectedOrder <- c(
    'Unvaccinated-Mock',
    'IV-BCG-Mock',
    'RhCMV-TB/9Ag-Mock',
    
    'RhCMV/Gag',
    'RhCMV-Malaria',
    
    'Unvaccinated',
    
    'BCG (adult)',
    'BCG (at birth)',
    'BCG (at birth) / RhCMV-TB/6Ag',
    'IV-BCG',
    
    '∆pp71 RhCMV-TB/6Ag',
    'RhCMV-TB/6Ag',
    'RhCMV-TB/9Ag'
  )
  
  for (l in rev(expectedOrder)) {
    if (l %in% unique(metadata$Group)) {
      metadata$Group <- forcats::fct_relevel(metadata$Group, l, after = 0)
    }
  }
  
  cDNA <- merge(cDNA, metadata, by = 'SubjectId', all.x = T)
  
  return(cDNA)
}



#' @title ApplyMalariaMetadata
#' @description Applies standard metadata related to the Wilder Malaria study
#'
#' @param seuratObj A Seurat object.
#' @param errorIfUnknownIdsFound If true, the function will fail if the seurat object contains unknown IDs
#' @param reApplyMetadata If true, QueryAndApplyCdnaMetadata will be re-run
#' @return A modified Seurat object.
#' @export
ApplyMalariaMetadata <- function(seuratObj, errorIfUnknownIdsFound = TRUE, reApplyMetadata = TRUE) {
  if (reApplyMetadata) {
    seuratObj <- .ApplyMetadata(seuratObj)
  }
  
  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1172",
    schemaName="lists",
    queryName="MalariaAnimals",
    colNameOpt="rname",
    colSelect = 'subjectid,sex,groupname,mmr,d0,mmr_date,cvac1,cvac2,cvac3,challengedate,PreExposureDate,Protection',
  )
  names(metadata) <- c('SubjectId', 'Sex', 'GroupName', 'MMR', 'D0', 'MMR_Date', 'CVac1', 'CVac2', 'CVac3', 'ChallengeDate', 'PreExposureDate', 'Protection')
  
  metadata2 <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1172",
    schemaName="lists",
    queryName="MalariaSamples",
    colNameOpt="rname",
    colSelect = 'libraryid,timepointlabel,LibraryId/sortId/sampleId/subjectId'
  )
  names(metadata2) <- c('cDNA_ID', 'TimepointLabel', 'SubjectId')
  
  metadata <- merge(metadata, metadata2, by = 'SubjectId', all.y = T)
  metadata <- metadata[names(metadata) != 'SubjectId']
  
  if (errorIfUnknownIdsFound && !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID)) {
    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }
  
  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- 1:nrow(toAdd)
  toAdd <- merge(toAdd, metadata, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))

  seuratObj$TimepointLabel <- naturalsort::naturalfactor(seuratObj$TimepointLabel)
  seuratObj$TimepointLabel <- forcats::fct_relevel(seuratObj$TimepointLabel, 'Chall_05', after = Inf)

  seuratObj$GroupName <- naturalsort::naturalfactor(seuratObj$GroupName)

  return(seuratObj)
}

#' @title ApplyPC531Metadata
#' @description Applies standard metadata related to PC531
#'
#' @param seuratObj A Seurat object.
#' @param errorIfUnknownIdsFound If true, the function will fail if the seurat object contains unknown IDs
#' @param reApplyMetadata If true, QueryAndApplyCdnaMetadata will be re-run
#' @return A modified Seurat object.
#' @export
ApplyPC531Metadata <- function(seuratObj, errorIfUnknownIdsFound = TRUE, reApplyMetadata = TRUE) {
  if (reApplyMetadata) {
    seuratObj <- .ApplyMetadata(seuratObj)
  }
  
  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1297",
    schemaName="lists",
    queryName="PC531_Subjects",
    colNameOpt="rname",
    colSelect = 'subjectid,VaccinationDate,InfectionDate,Outcome,VaccineType',
  )
  names(metadata) <- c('SubjectId', 'VaccinationDate', 'InfectionDate', 'Outcome', 'VaccineType')
  
  metadata2 <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1297",
    schemaName="lists",
    queryName="PC531",
    colNameOpt="rname",
    colSelect = 'cDNA_ID,Label,DPV,PID,cDNA_ID/sortId/sampleId/subjectId,pvl'
  )
  names(metadata2) <- c('cDNA_ID', 'TimepointLabel', 'DPV', 'PID', 'SubjectId', 'PVL')
  
  metadata <- merge(metadata, metadata2, by = 'SubjectId', all.y = T)
  metadata <- metadata[names(metadata) != 'SubjectId']
  
  if (errorIfUnknownIdsFound && !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID)) {
    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }

  if (any(duplicated(metadata$cDNA_ID))) {
    dups <- metadata$cDNA_ID[duplicated(metadata$cDNA_ID)]
    stop(paste0('There were duplicated cDNA_IDs in the metadata: ', paste0(dups, collapse = ',')))
  }

  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- 1:nrow(toAdd)
  toAdd <- merge(toAdd, metadata, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]
  
  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))
  
  return(seuratObj)
}

#' @title ApplyAcuteNxMetadata
#' @description Applies standard metadata related to the PPG/AcuteNx project
#'
#' @param seuratObj A Seurat object.
#' @param errorIfUnknownIdsFound If true, the function will fail if the seurat object contains unknown IDs
#' @param reApplyMetadata If true, QueryAndApplyCdnaMetadata will be re-run
#' @return A modified Seurat object.
#' @export
ApplyAcuteNxMetadata <- function(seuratObj, errorIfUnknownIdsFound = TRUE, reApplyMetadata = TRUE) {
  if (reApplyMetadata) {
    seuratObj <- .ApplyMetadata(seuratObj)
  }

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1297",
    schemaName="lists",
    queryName="AcuteNx_Subjects",
    colNameOpt="rname",
    colSelect="SubjectId,NxTimepoint,NxPVL,NxDate,InfectionDate,Dose,Route,Cohort"
  )
  names(metadata) <- c('SubjectId', 'NxTimepoint', 'NxPVL', 'NxDate', 'InfectionDate', 'Dose', 'Route', 'PcCohort')

  #if (errorIfUnknownIdsFound && !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID)) {
  #  missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
  #  stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  #}

  #if (any(duplicated(metadata$cDNA_ID))) {
  #  dups <- metadata$cDNA_ID[duplicated(metadata$cDNA_ID)]
  #  stop(paste0('There were duplicated cDNA_IDs in the metadata: ', paste0(dups, collapse = ',')))
  #}

  toAdd <- data.frame(SubjectId = seuratObj$SubjectId, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- 1:nrow(toAdd)
  toAdd <- merge(toAdd, metadata, by.x = 'SubjectId', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))

  return(seuratObj)
}

.SetFieldsToUnknown <- function(seuratObj, fieldNames) {
  for (fieldName in fieldNames) {
    sel <- is.na(seuratObj@meta.data[[fieldName]])
    if (any(sel)) {
      if (is.character(seuratObj@meta.data[[fieldName]])) {
        seuratObj@meta.data[[fieldName]][sel] <- 'Unknown'
      }
    }
  }

  return(seuratObj)
}

.ApplyMetadata <- function(seuratObj) {
  if ('BarcodePrefix' %in% names(seuratObj@meta.data) && !any(is.na(as.integer(seuratObj$BarcodePrefix)))) {
    seuratObj <- QueryAndApplyCdnaMetadata(seuratObj)
  } else if ('cDNA_ID' %in% names(seuratObj@meta.data)) {
    seuratObj <- QueryAndApplyMetadataUsingCDNA(seuratObj)
  } else {
    stop('Unable to find either BarcodePrefix or cDNA_ID in meta.data')
  }

  return(seuratObj)
}