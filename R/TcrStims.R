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
#' @param retainRowsWithoutCDR3 If true, rather than dropping rows without CDR3 data, these will be assigned as 'No TCR' as retained
#' @export
#' @import Rlabkey
#' @import dplyr
PrepareTcrData <- function(seuratObjOrDf, subjectId, minEDS = 0, enforceAllDataPresent = TRUE, chain = 'TRB', dropUnknownTNK_Type = FALSE, lowFreqThreshold = 0.001, retainRowsWithoutCDR3 = FALSE) {
  groupingFields <- c('cDNA_ID', 'SubjectId')

  if (typeof(seuratObjOrDf) == 'S4') {
    dat <- seuratObj@meta.data
  } else if (typeof(seuratObjOrDf) == 'list') {
    dat <- seuratObjOrDf
  } else {
    stop('Unknown object provided for seuratObjOrDf')
  }

  for (fieldName in c(groupingFields, 'IsActive')) {
    if (!fieldName %in% names(dat)) {
      stop(paste0('Missing field: ', fieldName))
    }
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
    print(paste0('cells after filter: ', nrow(dat), ', original: ', origCells))
  }

  # TODO: generalize!
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

  allStims$IsControlSample <- allStims$cDNA_ID %in% allStims$NoStimId
  noStimData <- allStims %>%
    filter(IsControlSample)

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

  dat$Clonotype <- dat[[chain]]
  dat$V_Gene <- dat[[paste0(chain, '_V')]]

  if (retainRowsWithoutCDR3) {
    dat$Clonotype <- as.character(dat$Clonotype)
    dat$Clonotype[is.na(dat$Clonotype)] <- 'No TCR'
  } else {
    origRows <- nrow(dat)
    dat <- dat %>% filter(!is.na(Clonotype))
    print(paste0('Dropping rows without CDR3. Before/after filter: ', origRows, ' / ', nrow(dat)))
  }

  if (dropUnknownTNK_Type) {
    if (! 'TNK_Type'  %in% names(dat)) {
      stop('Missing TNK_Type field')
    }

    origRows <- nrow(dat)
    dat <- dat %>%
      filter(TNK_Type != 'Ambiguous') %>%
      filter(TNK_Type != 'Unknown')

    print(paste0('Dropping Ambiguous/Unknown TNK_Type cells. Before/after filter: ', origRows, ' / ', nrow(dat)))
  }

  dat <- dat %>%
    group_by(across(all_of(groupingFields))) %>%
    mutate(TotalCellsForSample = n()) %>%
    ungroup() %>%
    group_by(across(all_of(c(groupingFields, 'IsActive')))) %>%
    mutate(TotalCellsForSampleAndState = n()) %>%
    ungroup() %>%
    group_by(across(all_of(c(groupingFields, 'Clonotype')))) %>%
    mutate(TotalCellsForClone = n()) %>%
    ungroup() %>%
    group_by(across(all_of(c(groupingFields, 'Clonotype', 'TotalCellsForSample', 'TotalCellsForSampleAndState', 'TotalCellsForClone', 'IsActive')))) %>%
    summarize(
      TotalCellsForCloneAndState = n(),
      V_Gene = paste0(sort(unique(V_Gene)), collapse = ',')
    ) %>%
    as.data.frame() %>%
    mutate(
      FractionOfCloneWithStateInSample = TotalCellsForCloneAndState / TotalCellsForSample,
      FractionOfCloneWithState = TotalCellsForCloneAndState / TotalCellsForClone,
      FractionOfSampleWithState = TotalCellsForSampleAndState / TotalCellsForSample
    ) %>%
    group_by(across(all_of(c('SubjectId', 'Clonotype')))) %>%
    mutate(MaxFractionInSubject = max(FractionOfCloneWithStateInSample))

  dat$OrigClonotype <- dat$Clonotype
  dat$Clonotype <- as.character(dat$Clonotype)
  dat$Clonotype[dat$MaxFractionInSubject < lowFreqThreshold] <- 'Low Freq'
  dat$Clonotype <- as.factor(dat$Clonotype)
  dat$Clonotype <- forcats::fct_reorder(dat$Clonotype, dat$MaxFractionInSubject, .desc = TRUE)

  dat$IsActiveLabel <- ifelse(dat$IsActive, yes = 'Activated', no = 'Not Activated')

  dat <- dat %>%
    group_by(across(all_of(c(groupingFields, 'Clonotype')))) %>%
    mutate(IsSharedAcrossStates = n_distinct(IsActive) > 1)
  dat$IsSharedAcrossStates[dat$Clonotype == 'Low Freq'] <- FALSE
  dat$IsSharedAcrossStates <- ifelse(dat$IsSharedAcrossStates, yes = 'Yes', no = 'No')

  dat <- dat %>%
    left_join(allStims, by = 'cDNA_ID')
  
  # Now merge with the active NoStim values:
  noStimSummary <- dat %>%
    as.data.frame() %>%
    filter(cDNA_ID %in% noStimData$cDNA_ID) %>%
    filter(Clonotype != 'Low Freq') %>%
    filter(IsActive) %>%
    select(cDNA_ID, IsActive, Clonotype, TotalCellsForClone, FractionOfCloneWithStateInSample) %>%
    rename(
      NoStimId = cDNA_ID,
      NoStimTotalCellsActive = TotalCellsForClone,
      NoStimFractionActive = FractionOfCloneWithStateInSample
    )

  dat <- dat %>%
    left_join(noStimSummary, by = c('NoStimId', 'Clonotype', 'IsActive'))

  # Now merge with the total NoStim values:
  noStimSummary <- dat %>%
    as.data.frame() %>%
    filter(cDNA_ID %in% noStimData$cDNA_ID) %>%
    filter(Clonotype != 'Low Freq') %>%
    group_by(cDNA_ID, Clonotype, TotalCellsForSample) %>%
    summarize(TotalCellsForClone = sum(TotalCellsForClone)) %>%
    as.data.frame() %>%
    mutate(FractionOfCloneInSample = TotalCellsForClone / TotalCellsForSample) %>%
    select(cDNA_ID, Clonotype, FractionOfCloneInSample, TotalCellsForClone) %>%
    unique() %>%
    rename(
      NoStimId = cDNA_ID,
      NoStimFractionOfCloneInSample = FractionOfCloneInSample,
      NoStimTotalCells = TotalCellsForClone
    )

  dat <- dat %>%
    left_join(noStimSummary, by = c('NoStimId', 'Clonotype'))

  # If this clone is present but not active, treat the active frequency as zero, not NA
  dat$NoStimTotalCells[is.na(dat$NoStimTotalCells)] <- 0
  dat$NoStimFractionOfCloneInSample[dat$NoStimTotalCells == 0] <- 0
  dat$NoStimFractionActive[is.na(dat$NoStimFractionActive) & dat$NoStimTotalCells == 0] <- 0
  dat$NoStimFractionActive[is.na(dat$NoStimFractionActive) & dat$NoStimTotalCells ] <- 0
  
  if (any(is.na(dat$NoStimTotalCells))) {
    stop('Found NAs in NoStimTotalCells')
  }
  
  if (any(is.na(dat$NoStimFractionOfCloneInSample))) {
    stop('Found NAs in NoStimFractionOfCloneInSample')
  }
  
  dat$NoStimTotalCellsActive[is.na(dat$NoStimTotalCellsActive) & !is.na(dat$NoStimTotalCells)] <- 0
  if (any(is.na(dat$NoStimTotalCellsActive))) {
    stop('Found NAs in NoStimTotalCellsActive')
  }

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
    left_join(meta, by = c('cDNA_ID', 'SubjectId')) %>%
    as.data.frame()

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
#' @param labelUsingCounts If true, the cell total will be listed above the bars
#' @export
#' @import dplyr
GenerateTcrPlot <- function(dat, xFacetField = NA, plotTitle = NULL, yFacetField = 'IsActiveLabel', patternField = 'IsShared', dropInactive = FALSE, labelUsingCounts = TRUE) {
  if (dropInactive) {
    dat <- dat %>%
      filter(IsActive)
  }

  dat$Clonotype <- naturalsort::naturalfactor(dat$Clonotype)
  dat$Clonotype <- forcats::fct_drop(dat$Clonotype)
  colorSteps <- max(min(length(unique(dat$Clonotype[dat$Clonotype != 'Low Freq'])), 9), 3)
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(colorSteps, 'Set1'))

  if ('No TCR' %in% dat$Clonotype) {
    dat$Clonotype <- forcats::fct_relevel(dat$Clonotype, 'No TCR', after = 0)
  }

  if ('Low Freq' %in% dat$Clonotype) {
    dat$Clonotype <- forcats::fct_relevel(dat$Clonotype, 'Low Freq', after = 0)
  }

  cols <- getPalette(length(unique(dat$Clonotype[! dat$Clonotype %in% c('Low Freq', 'No TCR')])))
  cols <- sample(cols, size = length(cols))
  if ('No TCR' %in% dat$Clonotype) {
    cols <- c('#FFFFFF', cols)
  }

  if ('Low Freq' %in% dat$Clonotype) {
    cols <- c('#ECECEC', cols)
  }
  names(cols) <- levels(dat$Clonotype)

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

  groupFields <- 'Stim'
  if (yFacetField != '.') {
    groupFields <- c(groupFields, yFacetField)
  }

  if (xFacetField != '.') {
    groupFields <- c(groupFields, xFacetField)
  }

  groupFields <- unique(groupFields)

  labelData <- dat %>%
    filter(IsActive) %>%
    group_by(across(all_of(groupFields)))

  if (nrow(labelData) == 0) {
    labelData$FractionOfCloneWithStateInSample <- numeric()
    labelData$TotalCellsForCloneAndState <- integer()
    labelData$TotalCellsForSample <- integer()
    labelData$LabelText <- character()
  } else {
    labelData <- labelData %>%
      summarize(FractionOfCloneWithStateInSample = sum(FractionOfCloneWithStateInSample), TotalCellsForCloneAndState = sum(TotalCellsForCloneAndState), TotalCellsForSample = max(TotalCellsForSample)) %>%
      as.data.frame() %>%
      mutate(LabelText = paste0(TotalCellsForCloneAndState, ' /\n', TotalCellsForSample))
  }

  if (nrow(dat) == 0) {
    print('No passing data, skipping plot')
    return()
  }

  PT <- ggplot(dat, aes(x = Stim, y = FractionOfCloneWithStateInSample)) +
    ggpattern::geom_col_pattern(aes(pattern = PatternField, fill = Clonotype), pattern_fill = "black",
                                color = 'black',
                                pattern_density = 0.2,
                                pattern_spacing = 0.05,
                                pattern_key_scale_factor = 0.6
    ) +
    ggpattern::scale_pattern_manual(values = patternValues) +
    scale_fill_manual(values = cols) +
    labs(y = 'Pct of Cells', x = '', fill = 'Clone') +
    scale_y_continuous(label = scales::percent, expand = expansion(add = c(0, min(0.02, max(dat$FractionOfCloneWithStateInSample[dat$IsActive])*0.2)))) +
    egg::theme_article(base_size = 14) +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if (labelUsingCounts) {
    PT <- PT + geom_text(data = labelData, aes(label = LabelText), position=position_dodge(width=0.9), vjust=-0.25, size = 3)
  }

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


#' @title GroupOverlappingClones
#' @description This will plot TCR data created by PrepareTcrData
#' @export
#'
#' @param dat A dataframe, generated by PrepareTcrData
#' @param groupingFields A list of fields, typically sample-variables, on which to group
#' @param dataMask An optional logical vector. If provided, the 'dat' dataframe is filtered to include only rows where dataMask is TRUE.
#'
GroupOverlappingClones <- function(dat, groupingFields, dataMask = NULL) {
  if (is.null(dataMask)) {
    dupes <- unique(grep(dat$Clonotype, pattern = ',', value = TRUE))
  } else {
    if (length(dataMask) != nrow(dat)) {
      stop('The dataMask should be a vector of the same length as the dat dataframe')
    }

    if (!is.logical(dataMask)) {
      stop('The dataMask should be a logical vector')
    }

    dupes <- unique(grep(dat$Clonotype[dataMask], pattern = ',', value = TRUE))
  }

  if (length(dupes) > 0) {
    print(paste0('Total multi-chain clonotypes: ', length(dupes)))
    dat$OriginalClone <- dat$Clonotype
    dat$Clonotype <- as.character(dat$Clonotype)
    for (cn in dupes) {
      tokens <- unlist(strsplit(cn, split = ','))
      for (c in tokens) {
        if (any(dat$Clonotype == c)) {
          print(paste0('Combining: ', c, ' into ', cn))
          dat$Clonotype[dat$Clonotype == c] <- cn
        }
      }
    }

    # Re-group:
    dat <- dat %>%
      group_by(across(all_of(unique(c(groupingFields, 'TotalCellsForSample', 'TotalCellsForSampleAndState', 'IsActive', 'Clonotype'))))) %>%
      summarise(
        TotalCellsForClone = sum(TotalCellsForClone),
        TotalCellsForCloneAndState = sum(TotalCellsForCloneAndState),
        NoStimTotalCellsActive = sum(NoStimTotalCellsActive),
        NoStimFractionOfCloneInSample = sum(NoStimFractionOfCloneInSample),
        NoStimTotalCells = sum(NoStimTotalCells),
        OrigCloneNames = paste0(sort(unique(OriginalClone)), collapse = '|'),
        V_Gene = paste0(sort(unique(V_Gene)), collapse = ',')
      ) %>%
      as.data.frame() %>%
      mutate(
        FractionOfCloneWithStateInSample = TotalCellsForCloneAndState / TotalCellsForSample,
        FractionOfCloneWithState = TotalCellsForCloneAndState / TotalCellsForClone,
        FractionOfSampleWithState = TotalCellsForSampleAndState / TotalCellsForSample,
        NoStimFractionActive = NoStimTotalCellsActive / NoStimTotalCells
      ) %>%
      group_by(across(all_of(c('SubjectId', 'Clonotype')))) %>%
      mutate(MaxFractionInSubject = max(FractionOfCloneWithStateInSample))

    dat <- dat %>%
      select(all_of(c(groupingFields, 'TotalCellsForSample', 'TotalCellsForSampleAndState', 'IsActive', 'Clonotype', 'TotalCellsForClone', 'TotalCellsForCloneAndState',
                      'NoStimTotalCellsActive', 'NoStimTotalCells', 'NoStimFractionOfCloneInSample', 'OrigCloneNames', 'V_Gene',
                      'FractionOfCloneWithStateInSample', 'FractionOfCloneWithState', 'FractionOfSampleWithState', 'NoStimFractionActive', 'MaxFractionInSubject'
      )))
  } else {
    print('No multi-CDR3 clonotypes present')
  }

  return(dat)
}


#' @title ApplyCloneFilters
#' @description This will filter TCR data created by PrepareTcrData
#'
#' @param dat A dataframe, generated by PrepareTcrData
#' @param minCellsPerClone Any clone with fewer than this many cells in the IsActive category will be filtered
#' @param minFoldChangeAboveNoStim The fold change will be calculated as the fraction of cells/clone active relative to the no-stim background. Any clone below this value will be filtered
#' @param minFractionOfCloneActive If any clone has fewer than this fraction of cells active, it will be filtered
#' @export
ApplyCloneFilters <- function(dat, minCellsPerClone = 2, minFoldChangeAboveNoStim = 2, minFractionOfCloneActive = 0.005) {
  dat$Filter <- NA

  if (any(dat$Clonotype == 'Low Freq')) {
    dat$Filter[dat$Clonotype == 'Low Freq'] <- 'Low Freq'
  }

  if (!is.na(minCellsPerClone)) {
    if ( ! 'TotalCellsForCloneAndState' %in% names(dat)) {
      stop('Missing field: TotalCellsForCloneAndState')
    }

    dat$Filter[dat$TotalCellsForCloneAndState < minCellsPerClone] <- 'Insufficient Cells'
  }

  if (!is.na(minFoldChangeAboveNoStim) && minFoldChangeAboveNoStim > 0) {
    if ( ! 'NoStimFractionActive' %in% names(dat)) {
      stop('Missing field: NoStimFractionActive')
    }

    foldChangeData <- dat$FractionOfCloneWithStateInSample / dat$NoStimFractionActive
    toFilter <- !is.na(foldChangeData) & foldChangeData < minFoldChangeAboveNoStim
    dat$Filter[toFilter] <- 'Below NoStim Background'
  }

  if (!is.na(minFractionOfCloneActive) && minFractionOfCloneActive > 0) {
    if ( ! 'FractionOfCloneWithState' %in% names(dat)) {
      stop('Missing field: FractionOfCloneWithState')
    }

    dat$Filter[dat$IsActive & dat$FractionOfCloneWithState < minFractionOfCloneActive] <- 'Below Min Fraction Active'
  }

  # Only apply to IsActive:
  dat$Filter[!dat$IsActive] <- NA

  print(sort(table(dat$Filter[dat$IsActive], useNA = 'ifany'), decreasing = TRUE))

  dat$IsFiltered <- !is.na(dat$Filter)
  dat <- as.data.frame(dat)

  return(dat)
}

#' @title CalculateClonotypeEnrichment
#' @description This will filter calculate enrichment of CDR3s in one population relative to a control
#'
#' @param dataToTest A dataframe, generated by PrepareTcrData. This should lack the control data.
#' @param controlData A dataframe, generated by PrepareTcrData, containing the control CDR3 data.
#' @param groupingField A field, present in dataToTest+controlData, used to group samples.
#' @param cloneFieldName A field, present in dataToTest+controlData, holding the clonotype value
#' @param contrastField A field, present in dataToTest+controlData, used to contrast data. Typically this is isActive.
#' @param cloneSizeField A field, present in dataToTest+controlData, denoting the number of cells from this clonotype
#' @param showProgress If TRUE, a progress bar will be shown
#' @export
CalculateClonotypeEnrichment <- function(dataToTest, controlData, groupingField = 'cDNA_ID', cloneFieldName = 'Clonotype', contrastField = 'IsActive', cloneSizeField = 'TotalCellsForCloneAndState', showProgress = FALSE) {
  for (fn in c(cloneFieldName, groupingField, contrastField, cloneSizeField)) {
    if (!fn %in% names(dataToTest)) {
      stop(paste0('Missing field in dataToTest: ', fn))
    }

    if (!fn %in% names(controlData)) {
      stop(paste0('Missing field in controlData: ', fn))
    }

    if (any(is.na(dataToTest[[fn]]))) {
      stop(paste0('dataToTest has NA values for: ', fn))
    }

    if (any(is.na(controlData[[fn]]))) {
      stop(paste0('controlData has NA values for: ', fn))
    }
  }

  if (length(unique(dataToTest$SubjectId)) != 1) {
    stop('Expected a single subjectId in dataToTest')
  }
  subjectId <- unique(dataToTest$SubjectId)

  if (length(unique(controlData$SubjectId)) != 1) {
    stop('Expected a single subjectId in controlData')
  }

  if (subjectId != unique(controlData$SubjectId)) {
    stop('SubjectId do not match between dataToTest and controlData')
  }

  if (any(is.na(dataToTest[[cloneFieldName]]))) {
    stop('Found NA values for clonotype in dataToTest')
  }

  if (any(is.na(controlData[[cloneFieldName]]))) {
    stop('Found NA values for clonotype in controlData')
  }

  uniqueClones <- sort(unique(c(dataToTest[[cloneFieldName]], controlData[[cloneFieldName]])))
  print(paste0('Total clonotypes: ', length(uniqueClones)))

  df <- rbind(
    dataToTest %>% select(all_of(c(cloneFieldName, groupingField, contrastField, cloneSizeField))) %>% mutate(Category = 'Experiment'),
    controlData %>% select(all_of(c(cloneFieldName, groupingField, contrastField, cloneSizeField))) %>% mutate(Category = 'Control')
  )

  df <- as.data.frame(lapply(df, rep, df[[cloneSizeField]]))

  df <- df %>% rename(
    'Clonotype' = cloneFieldName,
    'GroupName' = groupingField
  )

  ctlGroup <- unique(controlData[[groupingField]])
  if (length(ctlGroup) > 1) {
    stop('controlData can only have one value for groupingField')
  }

  # Set control group to intercept:
  df$GroupName <- as.factor(df$GroupName)
  df$GroupName <- forcats::fct_expand(df$GroupName, as.character(ctlGroup), after = 0)
  if (ctlGroup != levels(df$GroupName)[1]) {
    stop('intercept not properly set!')
  }

  doWork <- function(df, uniqueClones, pb = NULL) {
    results <- future.apply::future_lapply(seq_along(uniqueClones), future.seed = CellMembrane::GetSeed(), FUN = function(i) {
      if (!is.null(pb)) {
        #update the progress bar
        pb()
      }

      cdr3 <- uniqueClones[i]
      dat <- df %>%
        filter(Clonotype == cdr3)

      if (nrow(dat) <= 1) {
        print(paste0('Skipping CDR3: ', cdr3,' with 1 or fewer rows/subject'))

        return(data.frame(Clonotype = cdr3,
                          GroupName = dat$GroupName,
                          coefficients = NA,
                          p_values = NA,
                          error = "Too few cells"

        ))
      }

      if (length(unique(dat$GroupName)) == 1) {
        return(data.frame(Clonotype = cdr3,
                          GroupName = dat$GroupName,
                          coefficients = NA,
                          p_values = NA,
                          error = "Single Group"

        ))
      }

      model <- logistf::logistf(IsActive ~ GroupName, family = "binomial", data = dat)
      groupNames <- gsub(names(coef(model)), pattern = 'GroupName', replacement = '')
      groupNames[1] <- 'CONTROL'

      model_results <- data.frame(Clonotype = cdr3,
                                  GroupName = groupNames,
                                  coefficients = coef(model),
                                  p_values = model$prob,
                                  error = "No"
      )

      return(model_results)
    })
  }

  if (showProgress) {
    progressr::with_progress({
      pb <- progressr::progressor(steps = length(uniqueClones))

      results <- doWork(df = df, uniqueClones = uniqueClones, pb = pb)
    })
  } else {
    results <- doWork(df = df, uniqueClones = uniqueClones, pb = NULL)
  }

  results <- do.call(rbind, results)
  rownames(results) <- NULL

  results$GroupName[results$GroupName == 'CONTROL'] <- ctlGroup

  results$FDR <- NA
  for (groupName in unique(results$GroupName)) {
      results$FDR[results$GroupName == groupName] <- p.adjust(results$p_values[results$GroupName == groupName], method = "fdr")
  }

  results <- results %>%
    select(Clonotype, GroupName, coefficients, p_values, FDR) %>%
    unique()

  # Restore original datatype
  if (is.integer(dataToTest[[groupingField]])) {
    results$GroupName <- as.integer(results$GroupName)
  }

  names(results)[names(results) == 'GroupName'] <- groupingField
  names(results)[names(results) == 'Clonotype'] <- cloneFieldName

  return(results)
}

#' @title AppendClonotypeEnrichmentPVals
#' @description This will iterate the results from PrepareTcrData, calculating clones enriched in the activated cells relative to negative controls
#'
#' @param dat A dataframe, generated by PrepareTcrData.
#' @param showProgress If TRUE, a progress bar will be shown
#' @export
AppendClonotypeEnrichmentPVals <- function(dat, showProgress = FALSE) {
  dataWithPVal <- NULL
  for (ctlId in unique(dat$NoStimId)) {
    dataToTest <- dat %>%
      filter(NoStimId == ctlId & cDNA_ID != ctlId) %>%
      select(SubjectId, cDNA_ID, Clonotype, IsActive, TotalCellsForCloneAndState)

    controlData <- dat %>%
      filter(cDNA_ID == ctlId) %>%
      select(SubjectId, cDNA_ID, Clonotype, IsActive, TotalCellsForCloneAndState)

    if (nrow(dataToTest) == 0) {
      print(paste0('No rows, skipping: '))
      next
    }

    if (length(unique(dataToTest$SubjectId)) > 1){
      stop(paste0('More than one subject for: ', ctlId))
    }

    results <- Rdiscvr::CalculateClonotypeEnrichment(dataToTest, controlData, showProgress = showProgress)

    dataWithPVal <- rbind(dataWithPVal, dat %>% filter(IsActive) %>% left_join(results, by = c('cDNA_ID', 'Clonotype')))
  }

  # Set a floor:
  dataWithPVal$FDR[dataWithPVal$FDR == 0] <- 0.0001
  P1 <- ggplot(dataWithPVal %>% filter(!is.na(FDR) & !IsControlSample), aes(x = coefficients, y = -log10(FDR), color = Stim, label = Clonotype)) +
    geom_point() +
    ggrepel::geom_label_repel(show.legend = FALSE, size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    egg::theme_article() +
    xlab("Log Odds Ratio") +
    ylab("-log10(FDR)") +
    ggtitle('Clonotype Enrichment')

  print(P1)

  return(dataWithPVal)
}