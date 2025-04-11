#' @include LabKeySettings.R
#' @include Utils.R
#' @import utils



.ScoreFns <- list(
  TandNK_Activation_UCell = function(df) {
    if (! 'TandNK_Activation_UCell' %in% names(df)) {
      stop('Missing field: TandNK_Activation_UCell')
    }

    return(df$TandNK_Activation_UCell > 0.5)
  }
)

#' @title PrepareTcrData
#' @description This will group TCR data by the supplied chain, and calculate fields used for filtering
#'
#' @param seuratObjOrDf Either a Seurat object or the meta.data dataframe from an object
#' @param subjectId The subjectId to process
#' @param enforceAllDataPresent If true, the function will error unless all stims defined in lists.tcr_stims are present in the input object
#' @param chain The chain to summarize
#' @param minEDS If provided, cells with EDS less than this value will be discarded
#' @param dropUnknownTNK_Type If true, cells with Ambiguous or Unknown TNK_Type will be dropped
#' @param lowFreqThreshold Clones that never have a per-sample fraction above this value will be labeled as Low Freq.
#' @export
#' @import Rlabkey
#' @import dplyr
PrepareTcrData <- function(seuratObjOrDf, subjectId, minEDS = 0, enforceAllDataPresent = TRUE, chain = 'TRB', dropUnknownTNK_Type = FALSE, lowFreqThreshold = 0.001) {
  groupingFields <- c('cDNA_ID', 'SubjectId')

  if (typeof(seuratObjOrDf) == 'S4') {
    dat <- seuratObj@meta.data
  } else if (typeof(seuratObjOrDf) == 'list') {
    dat <- seuratObjOrDf
  } else {
    stop('Unknown object provided for seuratObjOrDf')
  }

  dat <- dat %>%
    filter(SubjectId == subjectId)

  dat$SubjectId <- naturalsort::naturalfactor(dat$SubjectId)
  if (nrow(dat) == 0) {
    stop(paste0('No records found for the subject: ', subjectId))
  }

  if (! chain %in% names(dat)) {
    stop(paste0('Missing chain: ', chain))
  }

  if (!is.null(minEDS) && minEDS > 0) {
    if (!'Tcell_EffectorDifferentiation' %in% names(dat)) {
      stop('Missing field Tcell_EffectorDifferentiation from input table')
    }

    print(paste0('Filtering data to EDS>=', minEDS))
    origCells <- nrow(dat)
    dat <- dat %>% filter(Tcell_EffectorDifferentiation >= minEDS)
    print(paste0('cells after filter:', nrow(dat), ', original: ', origCells))
  }

  for (fieldName in c('cDNA_ID')) {
    if (!fieldName %in% names(dat)) {
      stop(paste0('Missing field: ', fieldName))
    }
  }

  allStims <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/1297",
    schemaName="lists",
    queryName="TCR_Stims",
    colSelect="cDNA_ID,NoStimId",
    colFilter=makeFilter(
      c("cDNA_ID/sortId/sampleId/subjectId", "EQUALS", subjectId),
      c("cDNA_ID/readsetId/totalFiles", "GT", 0)
    ),
    colNameOpt="rname"
  ) %>%
    rename(
      cDNA_ID = 'cdna_id',
      NoStimId = 'nostimid'
    )

  print(paste0('Found ', nrow(allStims), ' known stims'))

  if (enforceAllDataPresent && nrow(allStims) == 0) {
    stop('No matching stims found')
  }

  if (! all(allStims$cDNA_ID %in% dat$cDNA_ID)) {
    missing <- unique(allStims$cDNA_ID)
    missing <- missing[! missing %in% dat$cDNA_ID]

    if (enforceAllDataPresent) {
      stop(paste0('Missing cDNA_IDs: ', paste0(sort(missing), collapse = ', ')))
    } else {
      warning(paste0('Missing some cDNA_IDs: ', paste0(sort(missing), collapse = ', ')))
    }
  }

  allStims$IsNoStim <- allStims$cDNA_ID %in% allStims$NoStimId
  noStimData <- allStims %>%
    filter(IsNoStim)

  if (! all(noStimData$cDNA_ID %in% dat$cDNA_ID)) {
    missing <- unique(noStimData$cDNA_ID)
    missing <- missing[! missing %in% dat$cDNA_ID]

    if (enforceAllDataPresent) {
      stop(paste0('Missing cDNA_IDs for NoStims: ', paste0(sort(missing), collapse = ', ')))
    } else {
      warning(paste0('Missing cDNA_IDs for NoStims: ', paste0(sort(missing), collapse = ', ')))
    }
  }

  if (any(duplicated(allStims$cDNA_ID))) {
    stop('Duplicates found in stim data')
  }

  dat$Label <- dat[[chain]]
  dat <- dat %>% filter(!is.na(Label))

  if (dropUnknownTNK_Type) {
    if (! 'TNK_Type'  %in% names(dat)) {
      stop('Missing TNK_Type field')
    }

    dat <- dat %>%
      filter(TNK_Type != 'Ambiguous') %>%
      filter(TNK_Type != 'Unknown')
  }

  dat <- dat %>%
    group_by(across(all_of(c(groupingFields)))) %>%
    mutate(TotalCellsForSample = n()) %>%
    group_by(across(all_of(c(groupingFields, 'IsActive')))) %>%
    mutate(TotalCellsForSampleByActivation = n()) %>%
    ungroup() %>%
    group_by(across(all_of(c(groupingFields, 'Label')))) %>%
    mutate(TotalCellsForClone = n()) %>%
    ungroup() %>%
    group_by(across(all_of(c(groupingFields, 'Label', 'TotalCellsForSample', 'TotalCellsForSampleByActivation', 'TotalCellsForClone', 'IsActive')))) %>%
    summarize(TotalCellsForCloneByActive = n()) %>%
    as.data.frame() %>%
    mutate(
      Fraction = TotalCellsForCloneByActive / TotalCellsForSample,
      FractionOfCloneActive = TotalCellsForCloneByActive / TotalCellsForClone,
      FractionOfSampleActive = TotalCellsForSampleByActivation / TotalCellsForSample
    ) %>%
    group_by(across(all_of(c('SubjectId', 'Label')))) %>%
    mutate(MaxFractionInSubject = max(Fraction))

  dat$Label <- as.character(dat$Label)
  dat$Label[dat$MaxFractionInSubject < lowFreqThreshold] <- 'Low Freq'
  dat$Label <- as.factor(dat$Label)
  dat$Label <- forcats::fct_reorder(dat$Label, dat$MaxFractionInSubject, .desc = TRUE)

  dat$IsActiveLabel <- ifelse(dat$IsActive, yes = 'Activated', no = 'Not Activated')

  dat <- dat %>%
    group_by(across(all_of(c(groupingFields, 'Label')))) %>%
    mutate(IsShared = n_distinct(IsActive) > 1)
  dat$IsShared[dat$Label == 'Low Freq'] <- FALSE
  dat$IsShared <- ifelse(dat$IsShared, yes = 'Yes', no = 'No')

  dat <- dat %>%
    left_join(allStims, by = 'cDNA_ID')
  
  # Now merge with the active NoStim values:
  noStimSummary <- dat %>%
    as.data.frame() %>%
    filter(cDNA_ID %in% noStimData$cDNA_ID) %>%
    filter(Label != 'Low Freq') %>%
    filter(IsActive) %>%
    select(cDNA_ID, IsActive, Label, TotalCellsForClone, Fraction) %>%
    rename(
      NoStimId = cDNA_ID,
      NoStimFractionActive = Fraction,
      NoStimCellsActive = TotalCellsForClone
    )

  dat <- dat %>%
    left_join(noStimSummary, by = c('NoStimId', 'Label', 'IsActive'))

  # Now merge with the NoStim values:
  noStimSummary <- dat %>%
    as.data.frame() %>%
    filter(cDNA_ID %in% noStimData$cDNA_ID) %>%
    filter(Label != 'Low Freq') %>%
    group_by(cDNA_ID, Label) %>%
    mutate(TotalCellsForClone = sum(TotalCellsForClone), Fraction = sum(Fraction)) %>%
    as.data.frame() %>%
    select(cDNA_ID, Label, Fraction, TotalCellsForClone) %>%
    unique() %>%
    rename(
      NoStimId = cDNA_ID,
      NoStimFraction = Fraction,
      NoStimCells = TotalCellsForClone
    )

  dat <- dat %>%
    left_join(noStimSummary, by = c('NoStimId', 'Label'))

  # If this clone is present but not active, treat the active frequency as zero, not NA
  dat$NoStimFractionActive[is.na(dat$NoStimFractionActive) & !is.na(dat$NoStimFraction)] <- 0
  dat$NoStimCellsActive[is.na(dat$NoStimCellsActive) & !is.na(dat$NoStimCells)] <- 0

  meta <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath="/Labs/Bimber/",
    schemaName="singlecell",
    queryName="cdna_libraries",
    colSelect="rowid,sortId/sampleId/subjectId,sortId/sampleId/sampledate,sortId/sampleId/stim,sortId/sampleId/assayType,sortId/sampleId/tissue,sortId/population",
    colFilter=makeFilter(
      c("rowid", "IN", paste0(unique(dat$cDNA_ID), collapse = ';'))
    ),
    colNameOpt="rname"
  ) %>%
    rename(
      cDNA_ID = 'rowid',
      SubjectId = 'sortid_sampleid_subjectid',
      SampleDate = 'sortid_sampleid_sampledate',
      Stim = 'sortid_sampleid_stim',
      AssayType = 'sortid_sampleid_assaytype',
      Tissue = 'sortid_sampleid_tissue',
      Population = 'sortid_population'
    )

  dat <- dat %>%
    left_join(meta, by = c('cDNA_ID', 'SubjectId'))

  return(dat)
}

#' @title GenerateTcrPlot
#' @description This will plot TCR data created by PrepareTcrData
#'
#' @param dat A dataframe, generated by PrepareTcrData
#' @param plotTitle An optional string used as the plot title
#' @param xFacetField An optional field used for faceting
#' @param yFacetField An optional field used for faceting
#' @param patternField An optional field used for the pattern aesthetic
#' @param dropInactive If true, cells where IsActive=FALSE will be dropped
#' @export
#' @import dplyr
GenerateTcrPlot <- function(dat, xFacetField = NA, plotTitle = NULL, yFacetField = 'IsActiveLabel', patternField = 'IsShared', dropInactive = FALSE) {
  if (dropInactive) {
    dat <- dat %>%
      filter(IsActive)
  }

  dat$Label <- forcats::fct_drop(dat$Label)
  colorSteps <- max(min(length(unique(dat$Label[dat$Label != 'Low Freq'])), 9), 3)
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(colorSteps, 'Set1'))

  if ('Low Freq' %in% dat$Label) {
    dat$Label <- forcats::fct_relevel(dat$Label, 'Low Freq', after = 0)
  }

  cols <- getPalette(length(unique(dat$Label[dat$Label != 'Low Freq'])))
  cols <- sample(cols, size = length(cols))
  if ('Low Freq' %in% dat$Label) {
    cols <- c('#ECECEC', cols)
  }
  names(cols) <- levels(dat$Label)

  if (!is.na(patternField)) {
    dat$PatternField <- dat[[patternField]]
    if (is.logical(dat$PatternField)) {
      dat$PatternField <- ifelse(dat$PatternField, yes = 'Yes', no = 'No')
    }
  } else {
    dat$PatternField <- NA
  }
  patternValues <- c('Yes' = 'stripe', 'No' = 'none')

  if (is.null(yFacetField) || is.na(yFacetField)) {
    yFacetField <- '.'
  }

  if (is.null(xFacetField) || is.na(xFacetField)) {
    xFacetField <- '.'
  }

  groupFields <- c('Stim')
  if (yFacetField != '.') {
    groupFields <- c(groupFields, yFacetField)
  }

  if (xFacetField != '.') {
    groupFields <- c(groupFields, xFacetField)
  }

  groupFields <- unique(groupFields)

  labelData <- dat %>%
    filter(IsActive) %>%
    group_by(across(all_of(c(groupFields)))) %>%
    summarize(Fraction = sum(Fraction), TotalCellsForCloneByActive = sum(TotalCellsForCloneByActive), TotalCellsForSample = max(TotalCellsForSample)) %>%
    as.data.frame() %>%
    mutate(LabelText = paste0(TotalCellsForCloneByActive, ' / ', TotalCellsForSample))

  PT <- ggplot(dat, aes(x = Stim, y = Fraction)) +
    ggpattern::geom_col_pattern(aes(pattern = PatternField, fill = Label), pattern_fill = "black",
                                color = 'black',
                                pattern_density = 0.2,
                                pattern_spacing = 0.05,
                                pattern_key_scale_factor = 0.6
    ) +
    ggpattern::scale_pattern_manual(values = patternValues) +
    scale_fill_manual(values = cols) +
    labs(y = 'Pct of Cells', x = '', fill = 'Clone') +
    geom_text(data = labelData, aes(label = LabelText), position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
    scale_y_continuous(label = scales::percent, expand = expansion(add = c(0, min(0.01, max(dat$FractionOfSampleActive[dat$IsActive])*0.2)))) +
    egg::theme_article(base_size = 14) +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if (xFacetField != '.' || yFacetField != '.') {
    wrap_by <- function(xFacetField, yFacetField) {
      facet_grid(vars(!!sym(yFacetField)), vars(!!sym(xFacetField)), scales = 'free', space = 'free_x')
    }

    PT <- PT + wrap_by(xFacetField, yFacetField)
  }

  if (!is.null(plotTitle)) {
    PT <- PT + ggtitle(label = plotTitle)
  }

  return(PT)
}


#' @title GenerateTcrPlot
#' @description This will plot TCR data created by PrepareTcrData
#'
#' @param dat A dataframe, generated by PrepareTcrData
#' @param minCellsPerClone Any clone with fewer than this many cells in the IsActive category will be filtered
#' @param minFoldChangeAboveNoStim The fold change will be calculated as the fraction of cells/clone active relative to the no-stim background. Any clone below this value will be filtered
#' @param minFractionOfCloneActive If any clone has fewer than this fraction of cells active, it will be filtered
#' @export
ApplyCloneFilters <- function(dat, minCellsPerClone = 2, minFoldChangeAboveNoStim = 2, minFractionOfCloneActive = 0.01) {
  dat$Filter <- NA

  if (any(dat$Label == 'Low Freq')) {
    dat$Filter[dat$Label == 'Low Freq'] <- 'Low Freq'
  }

  if (!is.na(minCellsPerClone)) {
    if ( ! 'TotalCellsForCloneByActive' %in% names(dat)) {
      stop('Missing field: TotalCellsForCloneByActive')
    }

    dat$Filter[dat$TotalCellsForCloneByActive < minCellsPerClone] <- 'Insufficient Cells'
  }

  if (!is.na(minFoldChangeAboveNoStim) && minFoldChangeAboveNoStim > 0) {
    if ( ! 'NoStimFractionActive' %in% names(dat)) {
      stop('Missing field: NoStimFractionActive')
    }

    foldChangeData <- dat$Fraction / dat$NoStimFractionActive
    toFilter <- !is.na(foldChangeData) & foldChangeData < minFoldChangeAboveNoStim
    dat$Filter[toFilter] <- 'Below NoStim Background'
  }

  if (!is.na(minFractionOfCloneActive) && minFractionOfCloneActive > 0) {
    if ( ! 'FractionOfCloneActive' %in% names(dat)) {
      stop('Missing field: FractionOfCloneActive')
    }

    dat$Filter[dat$FractionOfCloneActive < minFractionOfCloneActive] <- 'Below Min Fraction Active'
  }

  # Only apply to IsActive:
  dat$Filter[!dat$IsActive] <- NA

  print(sort(table(dat$Filter[dat$IsActive], useNA = 'ifany'), decreasing = TRUE))

  dat$IsFiltered <- !is.na(dat$Filter)

  return(dat)
}