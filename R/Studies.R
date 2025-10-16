
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
  seuratObj <- .AppendDemographics(seuratObj)

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/587",
    schemaName="lists",
    queryName="PC475",
    colNameOpt="rname",
    colSelect = 'libraryid,label,timepoint,wpi,dpi,LibraryId/sortId/sampleId/subjectId,LibraryId/sortId/sampleId/sampledate,Nx_pVL,SIV_RNA,SIV_DNA,CAVL_SampleName,ART_Initiation,ART_Release,ChallengeDate',
  )
  names(metadata) <- c('cDNA_ID', 'TimepointLabel', 'Timepoint', 'WPI', 'DPI', 'SubjectId', 'SampleDate', 'Nx_pVL', 'SIV_RNA', 'SIV_DNA', 'CAVL_SampleName', 'ART_Initiation', 'ART_Release', 'ChallengeDate')

  metadata$TimepointLabel <- naturalsort::naturalfactor(metadata$TimepointLabel)
  metadata$Timepoint <- naturalsort::naturalfactor(metadata$Timepoint)

  metadata$DaysPostArtInitiation <- as.integer(as.Date(metadata$SampleDate) - as.Date(metadata$ART_Initiation))
  metadata$WeeksPostArtInitiation <- round(metadata$DaysPostArtInitiation / 7, 1)
  metadata$DaysPostArtRelease <- as.integer(as.Date(metadata$SampleDate) - as.Date(metadata$ART_Release))

  metadata2 <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber",
    schemaName="laboratory",
    queryName="project_usage",
    colNameOpt="rname",
    colSelect = 'subjectid,groupname,project',
    colFilter = makeFilter(c('project', 'IN', 'PC475 / HIV Cure;PC567'))
  )
  names(metadata2) <- c('SubjectId', 'PC475_Group', 'project')
  metadata2$PC475_Group[metadata2$project == 'PC567'] <- 'ART-Only'
  metadata2 <- metadata2[names(metadata2) != 'project']

  metadata <- merge(metadata, metadata2, by = 'SubjectId', all.x = T)
  metadata <- metadata[! names(metadata) %in% c('SubjectId', 'SampleDate')]

  metadata$ViremicCategory <- NA
  metadata$ViremicCategory[is.na(metadata$ChallengeDate) & !is.na(metadata$ART_Initiation)] <- 'ART-Only'
  metadata$ViremicCategory[!is.na(metadata$ChallengeDate) & metadata$Nx_pVL == 1] <- 'Aviremic'
  metadata$ViremicCategory[!is.na(metadata$ChallengeDate) & metadata$Nx_pVL > 1] <- 'Viremic'
  metadata$ViremicCategory[!is.na(metadata$ChallengeDate) & metadata$PC475_Group == 'ON-ART'] <- 'On-ART'

  if (errorIfUnknownIdsFound && (any(is.na(seuratObj$cDNA_ID)) || !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID))) {
    if (any(is.na(seuratObj$cDNA_ID))) {
      stop('There were missing cDNA_IDs in the seurat object')
    }

    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }
  
  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- seq_len(nrow(toAdd))
  toAdd <- merge(toAdd, metadata, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal on metadata')
  }

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
  seuratObj <- .AppendDemographics(seuratObj)

  metadata <- .GetTbMetadata()

  if (errorIfUnknownIdsFound && (any(is.na(seuratObj$cDNA_ID)) || !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID))) {
    if (any(is.na(seuratObj$cDNA_ID))) {
      stop('There were missing cDNA_IDs in the seurat object')
    }

    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }

  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- seq_len(nrow(toAdd))
  toAdd <- merge(toAdd, metadata, by = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal on metadata')
  }

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
    colSelect="cDNA_ID,SampleType,TimepointLabel,CFU_Homogenate,CFU_Tissue,Imputed_CFU_Tissue,cDNA_ID/sortid/sampleid/subjectid,cDNA_ID/sortid/sampleid/tissue",
    colNameOpt="rname"
  )
  names(cDNA) <- c('cDNA_ID', 'SampleType', 'TimepointLabel', 'CFU_Homogenate', 'CFU_Tissue', 'Imputed_CFU_Tissue', 'SubjectId', 'Tissue')

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/627",
    schemaName="lists",
    queryName="TB_Studies",
    colSelect="subjectid,study,vaccine,vaccinationdate,challenge,ChallengeType,NecropsyDate,PID,PathScore,LungPathScore,ChallengeDate,R_Caudal_Lung_1_Disease,LeftLungPathScore,ChallengeSide,RCaudalPathScoreLongTerm,TotalRightLungPathScoreLongTerm,TotalLeftLungPathScoreLongTerm",
    colNameOpt="rname"
  )
  names(metadata) <- c('SubjectId', 'TB_Study', 'Vaccine', 'VaccinationDate', 'Challenge', 'ChallengeType', 'NecropsyDate', 'PID', 'PathScore', 'LungPathScore', 'ChallengeDate', 'R_Caudal_Lung_1_Disease', 'LeftLungPathScore', 'ChallengeSide', 'RCaudalPathScoreLongTerm', 'TotalRightLungPathScoreLongTerm', 'TotalLeftLungPathScoreLongTerm')

  # CFU_Tissue_Rescaled
  cDNA$CFU_Homogenate[is.na(cDNA$CFU_Homogenate)] <- 0
  cDNA$CFU_Homogenate <- as.numeric(cDNA$CFU_Homogenate)
  cDNA$CFU_Homogenate_Rescaled <- scales::rescale(asinh(cDNA$CFU_Homogenate + 1), to = c(0,1))

  cDNA$CFU_Tissue[is.na(cDNA$CFU_Tissue)] <- 0
  cDNA$CFU_Tissue <- as.numeric(cDNA$CFU_Tissue)
  cDNA$CFU_Tissue_Rescaled <- scales::rescale(asinh(cDNA$CFU_Tissue + 1), to = c(0,1))

  metadata$IsMockChallenged <- !is.na(metadata$Challenge) & metadata$Challenge == 'Mock-challenged'

  #Round to week:
  metadata$Timepoint <- metadata$PID
  metadata$Timepoint[!is.na(metadata$Timepoint)] <- round(metadata$Timepoint[!is.na(metadata$Timepoint)]/7, 0)*7
  metadata$Timepoint[!is.na(metadata$Timepoint)] <- paste0('Day ', metadata$Timepoint[!is.na(metadata$Timepoint)])

  metadata$ChallengeSide[is.na(metadata$ChallengeSide)] <- 'Unknown'

  # NOTE: make a simpler timepoint/mock field that lumps Mock-challenged samples into a different timepoint level
  if (any(!is.na(metadata$IsMockChallenged) & metadata$IsMockChallenged)) {
    metadata$Timepoint[!is.na(metadata$IsMockChallenged) & metadata$IsMockChallenged] <- paste0('Mock / ', metadata$Vaccine[!is.na(metadata$IsMockChallenged) & metadata$IsMockChallenged])
  }

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

  metadata$VaccineAndChallenge <- as.character(metadata$Vaccine)
  if ('Mock-challenged' %in% metadata$Challenge) {
    metadata$VaccineAndChallenge[!is.na(metadata$Challenge) & metadata$Challenge == 'Mock-challenged'] <- paste0(metadata$Vaccine[!is.na(metadata$Challenge) & metadata$Challenge == 'Mock-challenged'], '-Mock')
  }
  metadata$VaccineAndChallenge <- naturalsort::naturalfactor(metadata$VaccineAndChallenge)

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
    
    'Delta pp71 RhCMV-TB/6Ag',
    'RhCMV-TB/6Ag',
    'RhCMV-TB/9Ag'
  )
  
  for (l in rev(expectedOrder)) {
    if (l %in% unique(metadata$VaccineAndChallenge)) {
      metadata$VaccineAndChallenge <- forcats::fct_relevel(metadata$VaccineAndChallenge, l, after = 0)
    }
  }

  cDNA <- merge(cDNA, metadata, by = 'SubjectId', all.x = T)

  cDNA$IsChallengeSide <- case_when(
    cDNA$ChallengeSide == 'Left' & grepl(cDNA$Tissue, pattern = '-L$') ~ TRUE,
    cDNA$ChallengeSide == 'Right' & grepl(cDNA$Tissue, pattern = '-R$') ~ TRUE,
    .default = FALSE
  )

  cDNA$SideType <- ifelse(cDNA$IsChallengeSide, yes = 'Challenged Lung', 'Contralateral Lung')

  cDNA <- cDNA %>% select(-Tissue)

  #Make a baseline timepoint exception for the timepoints
  timepointLevels <- levels(cDNA$Timepoint)
  cDNA$Timepoint <- as.character(cDNA$Timepoint)
  cDNA[!is.na(cDNA$SampleType) & cDNA$SampleType == 'Baseline', 'Timepoint'] <- "Baseline"
  cDNA$Timepoint <- naturalsort::naturalfactor(cDNA$Timepoint, levels = c('Baseline', timepointLevels))

  cDNA$TimepointWithMock <- cDNA$Timepoint
  cDNA$TimepointWithMock <- forcats::fct_expand(cDNA$TimepointWithMock, 'Mock', after = 1)
  cDNA$TimepointWithMock[cDNA$Challenge == 'Mock-challenged'] <- 'Mock'
  cDNA$TimepointWithMock <- forcats::fct_drop(cDNA$TimepointWithMock)

  cDNA$VaccineGroup <- cDNA$Vaccine
  cDNA$VaccineGroup <- forcats::fct_expand(cDNA$VaccineGroup, 'RhCMV/TB', after = 4)
  cDNA$VaccineGroup <- forcats::fct_expand(cDNA$VaccineGroup, 'ID-BCG', after = 1)
  cDNA$VaccineGroup[cDNA$Vaccine %in% c("RhCMV-TB/6Ag", "RhCMV-TB/9Ag")] <- "RhCMV/TB"
  cDNA$VaccineGroup[cDNA$Vaccine %in% c("BCG (adult)", "BCG (at birth)")] <- "ID-BCG"
  cDNA$VaccineGroup <- forcats::fct_drop(cDNA$VaccineGroup)

  cDNA$TB_Immunized <- ifelse(grepl(cDNA$VaccineGroup, pattern = 'BCG') | grepl(cDNA$VaccineGroup, pattern = '/TB'), yes = 'TB-Immunized', no = 'Not Immunized')

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
  seuratObj <- .AppendDemographics(seuratObj)
  
  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1172",
    schemaName="lists",
    queryName="MalariaAnimals",
    colNameOpt="rname",
    colSelect = 'subjectid,sex,groupname,mmr,d0,mmr_date,cvac1,cvac2,cvac3,challengedate,PreExposureDate,ChallengeType,Protection',
  )
  names(metadata) <- c('SubjectId', 'Sex', 'GroupName', 'MMR', 'D0', 'MMR_Date', 'CVac1', 'CVac2', 'CVac3', 'ChallengeDate', 'PreExposureDate', 'ChallengeType', 'Protection')
  
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

  if (errorIfUnknownIdsFound && (any(is.na(seuratObj$cDNA_ID)) || !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID))) {
    if (any(is.na(seuratObj$cDNA_ID))) {
      stop('There were missing cDNA_IDs in the seurat object')
    }

    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }
  
  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- seq_len(nrow(toAdd))
  toAdd <- merge(toAdd, metadata, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal on metadata')
  }

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))

  seuratObj$TimepointLabel <- as.character(seuratObj$TimepointLabel)
  seuratObj$TimepointLabel[seuratObj$TimepointLabel == 'Baseline-2' & seuratObj$GroupName == 'Malaria Preexposed'] <- 'Challenge-1'
  seuratObj$TimepointLabel <- naturalsort::naturalfactor(seuratObj$TimepointLabel)
  seuratObj$TimepointLabel <- forcats::fct_relevel(seuratObj$TimepointLabel, 'Chall_05', after = Inf)
  seuratObj$TimepointLabel <- forcats::fct_relevel(seuratObj$TimepointLabel, 'Chall_12', after = Inf)
  seuratObj$TimepointLabel <- forcats::fct_relevel(seuratObj$TimepointLabel, 'Challenge-1', after = 3)

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
  seuratObj <- .AppendDemographics(seuratObj)
  
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
  
  if (errorIfUnknownIdsFound && (any(is.na(seuratObj$cDNA_ID)) || !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID))) {
    if (any(is.na(seuratObj$cDNA_ID))) {
      stop('There were missing cDNA_IDs in the seurat object')
    }
    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }

  if (any(duplicated(metadata$cDNA_ID))) {
    dups <- metadata$cDNA_ID[duplicated(metadata$cDNA_ID)]
    stop(paste0('There were duplicated cDNA_IDs in the metadata: ', paste0(dups, collapse = ',')))
  }

  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- seq_len(nrow(toAdd))
  toAdd <- merge(toAdd, metadata, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal on metadata')
  }

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
  seuratObj <- .AppendDemographics(seuratObj)

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1297",
    schemaName="lists",
    queryName="AcuteNx_Subjects",
    colNameOpt="rname",
    colSelect="SubjectId,NxTimepoint,NxPVL,NxDate,InfectionDate,Dose,Route,Cohort"
  )
  names(metadata) <- c('SubjectId', 'NxTimepoint', 'NxPVL', 'NxDate', 'InfectionDate', 'Dose', 'Route', 'PcCohort')

  if (errorIfUnknownIdsFound && any(is.na(seuratObj$cDNA_ID))) {
    stop('There were blank cDNA_IDs in the seurat object')
  }

  if (any(duplicated(metadata$cDNA_ID))) {
   dups <- metadata$cDNA_ID[duplicated(metadata$cDNA_ID)]
   stop(paste0('There were duplicated cDNA_IDs in the metadata: ', paste0(dups, collapse = ',')))
  }

  toAdd <- data.frame(SubjectId = seuratObj$SubjectId, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- seq_len(nrow(toAdd))
  toAdd <- merge(toAdd, metadata, by.x = 'SubjectId', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]
  toAdd$NxDate <- as.Date(toAdd$NxDate)

  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal on metadata')
  }

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj$SampleDate <- as.Date(seuratObj$SampleDate)
  seuratObj$Timepoint <- dplyr::case_when(
    !is.na(seuratObj$NxDate) & seuratObj$SampleDate == seuratObj$NxDate ~ 'Necropsy',
    !is.na(seuratObj$NxDate) & seuratObj$SampleDate > seuratObj$NxDate ~ 'POST-NX!',
    !is.na(seuratObj$InfectionDate) & seuratObj$SampleDate > seuratObj$InfectionDate ~ 'Post-challenge',
    !is.na(seuratObj$InfectionDate) & seuratObj$SampleDate <= seuratObj$InfectionDate ~ 'Pre-infection',
    .default = 'UNKNOWN'
  )

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
  if ('BarcodePrefix' %in% names(seuratObj@meta.data) && !any(is.na(as.integer(as.character(seuratObj$BarcodePrefix))))) {
    seuratObj <- QueryAndApplyCdnaMetadata(seuratObj)
  } else if ('cDNA_ID' %in% names(seuratObj@meta.data)) {
    seuratObj <- QueryAndApplyMetadataUsingCDNA(seuratObj)
  } else {
    stop('Unable to find either BarcodePrefix or cDNA_ID in meta.data')
  }

  seuratObj$SampleDate <- as.Date(seuratObj$SampleDate)

  return(seuratObj)
}

.AppendDemographics <- function(seuratObj) {
  if (! 'SubjectId' %in% names(seuratObj@meta.data)) {
    stop('Missing SubjectId column!')
  }

  dat <- suppressWarnings(labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Internal/PMR",
    schemaName="study",
    queryName="Demographics",
    colSelect="Id,gender,species,birth,death",
    colFilter=makeFilter(c('Id', 'IN', paste0(unique(seuratObj$SubjectId), collapse = ';'))),
    containerFilter=NULL,
    colNameOpt="rname"
  ))
  names(dat) <- c('SubjectId', 'Sex', 'Species', 'Birth', 'Death')

  toAdd <- seuratObj@meta.data[,'SubjectId', drop = FALSE]
  toAdd$SortOrder <- seq_len(ncol(seuratObj))
  toAdd$CellBarcode <- rownames(seuratObj@meta.data)

  toAdd <- merge(toAdd, dat, by = 'SubjectId', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal')
  }

  if (any(!is.na(seuratObj$SubjectId) & !is.na(toAdd$SubjectId) & toAdd$SubjectId != seuratObj$SubjectId)) {
    stop('SubjectId not equal when adding metadata')
  }

  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'SubjectId')]
  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)

  return(seuratObj)
}

#' @title ApplyEC_Metadata
#' @description Applies standard metadata related to JS46 / EC Project
#'
#' @param seuratObj A Seurat object.
#' @param errorIfUnknownIdsFound If true, the function will fail if the seurat object contains unknown IDs
#' @param reApplyMetadata If true, QueryAndApplyCdnaMetadata will be re-run
#' @return A modified Seurat object.
#' @export
ApplyEC_Metadata <- function(seuratObj, errorIfUnknownIdsFound = TRUE, reApplyMetadata = TRUE) {
  if (reApplyMetadata) {
    seuratObj <- .ApplyMetadata(seuratObj)
  }
  seuratObj <- .AppendDemographics(seuratObj)

  cDNA <- seuratObj@meta.data %>%
    select(cDNA_ID, SubjectId, SampleDate) %>%
    unique()

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/Collaborations/EC_Project",
    schemaName="study",
    queryName="Demographics",
    colSelect="Id,mhcGenotypes/A01,mhcGenotypes/A02,mhcGenotypes/B08,mhcGenotypes/B17,outcomes/outcomes,sivART/infectionDate",
    colNameOpt="rname"
  ) %>%
    rename(
      SubjectId = id,
      Outcome = 'outcomes_outcomes',
      A01 = 'mhcgenotypes_a01',
      A02 = 'mhcgenotypes_a02',
      B08 = 'mhcgenotypes_b08',
      B17 = 'mhcgenotypes_b17',
      InfectionDate = 'sivart_infectiondate'
    )

  metadata$Genotype <- case_when(
    !is.na(metadata$B08) & !is.na(metadata$B17) ~ 'B08-B17',
    !is.na(metadata$B08) ~ 'B08',
    !is.na(metadata$B17) ~ 'B17',
    .default = 'OTHER'
  )

  if (errorIfUnknownIdsFound && (any(is.na(seuratObj$SubjectId)) || !all(seuratObj$SubjectId %in% metadata$SubjectId))) {
    if (any(is.na(seuratObj$SubjectId))) {
      stop('There were missing SubjectIds in the seurat object')
    }
    missing <- sort(unique(seuratObj$SubjectId[!seuratObj$SubjectId %in% metadata$SubjectId]))
    stop(paste0('There were SubjectId in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }

  cDNA <- cDNA %>%
    left_join(metadata, by = 'SubjectId')

  metadata2 <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/Collaborations/EC_Project",
    schemaName="study",
    queryName="additionalDatatypes",
    viewName="",
    colSelect="Id,date,timePostSivChallenge/daysPostInfection",
    colSort="Id",
    colFilter=makeFilter(c("category", "EQUAL", "scRNA-seq")),
    containerFilter=NULL,
    colNameOpt="rname"
  ) %>%
    rename(
      SubjectId = 'id',
      SampleDate = 'date',
      DPI = 'timepostsivchallenge_dayspostinfection'
    )

  metadata2$Timepoint <- case_when(
    metadata2$DPI < 0 ~ 'Baseline',
    metadata2$DPI == 0 ~ 'D0',
    metadata2$DPI < 35 ~ 'Acute',
    .default = 'OTHER'
  )

  cDNA <- cDNA %>%
    left_join(metadata2, by = c('SubjectId', 'SampleDate'))

  if (any(duplicated(cDNA$cDNA_ID))) {
    dups <- cDNA$cDNA_ID[duplicated(cDNA$cDNA_ID)]
    stop(paste0('There were duplicated cDNA_IDs in the metadata: ', paste0(dups, collapse = ',')))
  }

  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- seq_len(nrow(toAdd))
  toAdd <- merge(toAdd, cDNA, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId', 'SampleDate')]

  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal on metadata')
  }

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))

  return(seuratObj)
}

#' @title ApplyPPG_Stim_Metadata
#' @description Applies standard metadata related to the PPG T Cell Stim Data
#'
#' @param seuratObj A Seurat object.
#' @param errorIfUnknownIdsFound If true, the function will fail if the seurat object contains unknown IDs
#' @param reApplyMetadata If true, QueryAndApplyCdnaMetadata will be re-run
#' @return A modified Seurat object.
#' @export
ApplyPPG_Stim_Metadata <- function(seuratObj, errorIfUnknownIdsFound = TRUE, reApplyMetadata = TRUE) {
  if (reApplyMetadata) {
    seuratObj <- .ApplyMetadata(seuratObj)
  }
  seuratObj <- .AppendDemographics(seuratObj)

  metadata <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1297",
    schemaName="lists",
    queryName="PPG_TCR_Stim_Subjects",
    colNameOpt="rname",
    colSelect = 'subjectid,vaccinetype,challengedate,vaccinationdate,orfs,purpose',
  )
  names(metadata) <- c('SubjectId', 'VaccineType', 'ChallengeDate', 'VaccinationDate', 'ORFs', 'StimCohort')

  metadata2 <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber",
    schemaName="tcrdb",
    queryName="stims",
    colNameOpt="rname",
    colSelect = 'cdna_id,controlStimId,flowQuantification,cDNA_ID/sortId/sampleId/subjectId'
  )
  names(metadata2) <- c('cDNA_ID', 'NoStim_cDNA_ID', 'IFNG_TNF_ICS', 'SubjectId')

  metadata <- merge(metadata, metadata2, by = 'SubjectId', all.y = T)
  metadata <- metadata[names(metadata) != 'SubjectId']

  if (errorIfUnknownIdsFound && (any(is.na(seuratObj$cDNA_ID)) || !all(seuratObj$cDNA_ID %in% metadata$cDNA_ID))) {
    if (any(is.na(seuratObj$cDNA_ID))) {
      stop('There were missing cDNA_IDs in the seurat object')
    }
    missing <- sort(unique(seuratObj$cDNA_ID[!seuratObj$cDNA_ID %in% metadata$cDNA_ID]))
    stop(paste0('There were cDNA_IDs in the seurat object missing from the metadata, missing: ', paste0(missing, collapse = ',')))
  }

  if (any(duplicated(metadata$cDNA_ID))) {
    dups <- metadata$cDNA_ID[duplicated(metadata$cDNA_ID)]
    stop(paste0('There were duplicated cDNA_IDs in the metadata: ', paste0(dups, collapse = ',')))
  }

  toAdd <- data.frame(cDNA_ID = seuratObj$cDNA_ID, CellBarcode = colnames(seuratObj))
  toAdd$SortOrder <- seq_len(nrow(toAdd))
  toAdd <- merge(toAdd, metadata, by.x = 'cDNA_ID', all.x = TRUE)
  toAdd <- arrange(toAdd, SortOrder)
  rownames(toAdd) <- toAdd$CellBarcode
  toAdd <- toAdd[!names(toAdd) %in% c('CellBarcode', 'SortOrder', 'cDNA_ID', 'SubjectId')]

  if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
    stop('Row names not equal on metadata')
  }

  seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  seuratObj <- .SetFieldsToUnknown(seuratObj, names(toAdd))

  return(seuratObj)
}
