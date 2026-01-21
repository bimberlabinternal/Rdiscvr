#' @include LabKeySettings.R
#' @include Utils.R
#' @import utils


utils::globalVariables(
  names = c('Clonotype','FDR','FractionOfCloneInSample','FractionOfCloneWithStateInSample','GroupName','IsControlSample','LabelText',
            'NoStimFractionOfCloneInSample','NoStimId','NoStimTotalCells','NoStimTotalCellsActive','OriginalClone','PatternField','Stim','TNK_Type',
            'Tcell_EffectorDifferentiation','TotalCellsForClone','TotalCellsForCloneAndState','TotalCellsForSample','TotalCellsForSampleAndState','V_Gene', 'J_Gene', 'cdr3WithSegments',
            'cDNA_ID','coefficients', 'p_val', 'error', 'FractionOfCloneWithState', 'Antigens', 'Chain', 'ChainsForAntigenMatch', 'HasIE', 'HasNoStim', 'IsIE', 'IsNoStim',
            'SampleDate', 'Tissue', 'fractioncloneactivated', 'maxFractionCloneActivated', 'maxTotalCloneSize', 'meanCloneSize', 'meanFractionCloneActivated', 'totalclonesize',
            'TandNK_ActivationCore_UCell', 'TandNK_Activation_UCell', 'TandNK_Activation3_UCell', 'IsFiltered', 'FailedEnrichment', 'FractionOfSampleWithState', 'container'),
  package = 'Rdiscvr',
  add = TRUE
)


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
#' @param enforceAllDataPresent If true, the function will error unless all stims listed in tcrdb.stims for the subjectId are present in the input object
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
    dat <- seuratObjOrDf@meta.data
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

  for (suffix in c('_V', '_J', '_Segments')) {
    fn <- paste0(chain, suffix)
    if (! fn %in% names(dat)) {
      stop(paste0('Missing field: ', fn))
    }
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

  allStims <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="tcrdb",
    queryName="stims",
    colSelect="cdna_id,controlStimId,status",
    colFilter=makeFilter(
      c("cDNA_ID/sortId/sampleId/subjectId", "EQUALS", subjectId),
      c("cDNA_ID/readsetId/totalFiles", "GT", 0)
    ),
    colNameOpt="rname"
  ) %>%
    rename(
      cDNA_ID = 'cdna_id',
      NoStimId = 'controlstimid',
      StimStatus = 'status'
    )

  print(paste0('Found ', nrow(allStims), ' known stims'))

  if (enforceAllDataPresent && nrow(allStims) == 0) {
    stop('No matching stims found')
  }

  if (! all(allStims$cDNA_ID %in% dat$cDNA_ID)) {
    missing <- unique(allStims$cDNA_ID)
    missing <- missing[! missing %in% dat$cDNA_ID]

    # Allow failed stims to be missing:
    failed <- allStims$cDNA_ID[grepl(allStims$StimStatus, pattern = 'Fail')]
    missing <- missing[! missing %in% failed]

    if (length(missing) > 0) {
      if (enforceAllDataPresent) {
        stop(paste0('Missing cDNA_IDs: ', paste0(sort(missing), collapse = ', ')))
      } else {
        warning(paste0('Missing some cDNA_IDs: ', paste0(sort(missing), collapse = ', ')))
      }
    }
  }

  allStims$IsControlSample <- allStims$cDNA_ID %in% allStims$NoStimId
  noStimData <- allStims %>%
    filter(IsControlSample) %>%
    filter(!grepl(StimStatus, pattern = 'Fail'))

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
  dat$J_Gene <- dat[[paste0(chain, '_J')]]
  dat$cdr3WithSegments <- dat[[paste0(chain, '_Segments')]]

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
      V_Gene = paste0(sort(unique(V_Gene)), collapse = ','),
      J_Gene = paste0(sort(unique(J_Gene)), collapse = ','),
      cdr3WithSegments = paste0(sort(unique(cdr3WithSegments)), collapse = ','),
    ) %>%
    as.data.frame() %>%
    mutate(
      FractionOfCloneWithStateInSample = TotalCellsForCloneAndState / TotalCellsForSample,
      FractionOfCloneWithState = TotalCellsForCloneAndState / TotalCellsForClone,
      FractionOfSampleWithState = TotalCellsForSampleAndState / TotalCellsForSample
    ) %>%
    group_by(across(all_of(c('SubjectId', 'Clonotype')))) %>%
    mutate(MaxFractionInSubject = max(FractionOfCloneWithStateInSample)) %>%
    as.data.frame()

  dat$V_Gene <- sapply(dat$V_Gene, function(x){
    return(paste0(sort(unique(unlist(strsplit(x, split = ',')))), collapse = ','))
  })

  dat$J_Gene <- sapply(dat$J_Gene, function(x){
    return(paste0(sort(unique(unlist(strsplit(x, split = ',')))), collapse = ','))
  })

  dat$cdr3WithSegments <- sapply(dat$cdr3WithSegments, function(x){
    return(paste0(sort(unique(unlist(strsplit(x, split = ',')))), collapse = ','))
  })

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

  dat <- dat %>%
    filter(!grepl(StimStatus, pattern = 'Fail'))

  if (nrow(dat) == 0) {
    stop('No rows remain after dropping failed stims')
  }

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
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
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
#' @param groupLowFreq If true, all rows where clonotype is 'Low Freq' will be grouped to simplify plotting
#' @export
#' @import dplyr
GenerateTcrPlot <- function(dat, xFacetField = NA, plotTitle = NULL, yFacetField = 'IsActiveLabel', patternField = 'IsShared', dropInactive = FALSE, labelUsingCounts = TRUE, groupLowFreq = FALSE) {
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

  if (groupLowFreq && 'Low Freq' %in% dat$Clonotype) {
    fields <- unique(c(groupFields, 'PatternField', 'Clonotype'))
    lf <- dat %>%
      filter(Clonotype == 'Low Freq') %>%
      group_by(across(all_of(fields))) %>%
      summarise(FractionOfCloneWithStateInSample = sum(FractionOfCloneWithStateInSample))

    dat <- dat %>%
      filter(Clonotype != 'Low Freq')

    fields <- c(fields, 'FractionOfCloneWithStateInSample')
    lf <- lf %>%
      select(all_of(fields))

    dat <- dat %>%
      select(all_of(fields))

    dat <- rbind(dat, lf)
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
    scale_y_continuous(labels = scales::percent, expand = expansion(add = c(0, min(0.02, max(dat$FractionOfCloneWithStateInSample[dat$IsActive])*0.2)))) +
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
#' @param maxRatioToCombine If provided, the function will calculate the total number of cells in the parent group divided by the total number of cells in the group to combine. It will only merge groups is that ratio is above this value
#' @param dataMask An optional logical vector. If provided, the 'dat' dataframe is filtered to include only rows where dataMask is TRUE.
#'
GroupOverlappingClones <- function(dat, groupingFields, maxRatioToCombine = 0.5, dataMask = NULL) {
  origClonotypes <- NULL
  if (is.null(dataMask)) {
    dupes <- unique(grep(dat$Clonotype, pattern = ',', value = TRUE))
    origClonotypes <- unique(dat$Clonotype)
  } else {
    if (length(dataMask) != nrow(dat)) {
      stop('The dataMask should be a vector of the same length as the dat dataframe')
    }

    if (!is.logical(dataMask)) {
      stop('The dataMask should be a logical vector')
    }

    dupes <- unique(grep(dat$Clonotype[dataMask], pattern = ',', value = TRUE))
    origClonotypes <- unique(dat$Clonotype[dataMask])
  }

  joinedClones <- FALSE
  if (length(dupes) > 0) {
    print(paste0('Total multi-chain clonotypes: ', length(dupes)))
    dat$OriginalClone <- dat$Clonotype
    dat$Clonotype <- as.character(dat$Clonotype)
    for (cn in dupes) {
      tokens <- unlist(strsplit(cn, split = ','))
      for (c in tokens) {
        if (any(dat$Clonotype == c)) {
          sum1 <- sum(dat$TotalCellsForCloneAndState[dat$Clonotype == c])
          sum2 <- sum(dat$TotalCellsForCloneAndState[dat$Clonotype == cn])
          rat <- sum1 / sum2

          doCombine <- is.na(maxRatioToCombine) || is.null(maxRatioToCombine) || rat < maxRatioToCombine
          prefix <- ifelse(doCombine, yes = 'Combining', no = 'Skipping')
          print(paste0(prefix, ': ', c, ' (', sum1, ') into ', cn, ' (', sum2, '). Ratio: ', rat))
          if (doCombine) {
            joinedClones <- TRUE
            dat$Clonotype[dat$Clonotype == c] <- cn
          }
        }
      }
    }

    # Re-group. Separate Low Freq so we dont artificially lump these:
    lowFreq <- dat %>% filter(Clonotype == 'Low Freq')
    dat <- dat %>%
      filter(Clonotype != 'Low Freq') %>%
      group_by(across(all_of(unique(c(groupingFields, 'TotalCellsForSample', 'TotalCellsForSampleAndState', 'IsActive', 'Clonotype'))))) %>%
      summarise(
        TotalCellsForClone = sum(TotalCellsForClone),
        TotalCellsForCloneAndState = sum(TotalCellsForCloneAndState),
        NoStimTotalCellsActive = sum(NoStimTotalCellsActive),
        NoStimFractionOfCloneInSample = sum(NoStimFractionOfCloneInSample),
        NoStimTotalCells = sum(NoStimTotalCells),
        OrigClonotype = paste0(sort(unique(OriginalClone)), collapse = '|'),
        V_Gene = paste0(sort(unique(V_Gene)), collapse = ','),
        J_Gene = paste0(sort(unique(J_Gene)), collapse = ','),
        cdr3WithSegments = paste0(sort(unique(cdr3WithSegments)), collapse = ',')
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

    dat$V_Gene <- sapply(dat$V_Gene, function(x){
      return(paste0(sort(unique(unlist(strsplit(x, split = ',')))), collapse = ','))
    })

    dat$J_Gene <- sapply(dat$J_Gene, function(x){
      return(paste0(sort(unique(unlist(strsplit(x, split = ',')))), collapse = ','))
    })

    dat$cdr3WithSegments <- sapply(dat$cdr3WithSegments, function(x){
      return(paste0(sort(unique(unlist(strsplit(x, split = ',')))), collapse = ','))
    })

    dat <- dat %>%
      select(all_of(c(groupingFields, 'TotalCellsForSample', 'TotalCellsForSampleAndState', 'IsActive', 'Clonotype', 'TotalCellsForClone', 'TotalCellsForCloneAndState',
                      'NoStimTotalCellsActive', 'NoStimTotalCells', 'NoStimFractionOfCloneInSample', 'OrigClonotype', 'V_Gene', 'J_Gene', 'cdr3WithSegments',
                      'FractionOfCloneWithStateInSample', 'FractionOfCloneWithState', 'FractionOfSampleWithState', 'NoStimFractionActive', 'MaxFractionInSubject'
      )))

    lowFreq <- lowFreq %>% select(all_of(names(dat)))
    dat <- rbind(dat, lowFreq)
  } else {
    print('No multi-CDR3 clonotypes present')
  }

  if (joinedClones) {
    print('Some clones were merged, so repeating GroupOverlappingClones()')
    if (!all(is.null(dataMask))) {
      # This needs to be reset since the row length is changed:
      dataMask <- dat$Clonotype %in% origClonotypes
    }
    return(GroupOverlappingClones(dat = dat, groupingFields = groupingFields, maxRatioToCombine = maxRatioToCombine, dataMask = dataMask))
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
  dat$IsFiltered <- !is.na(dat$Filter)
  dat <- as.data.frame(dat)

  return(dat)
}

#' @title CalculateClonotypeEnrichment
#' @description This will filter calculate enrichment of CDR3s in one population relative to a control
#'
#' @param dataToTest A dataframe, such as generated by PrepareTcrData. This should lack the control data, with one row per clonotype
#' @param controlData A dataframe, such as generated by PrepareTcrData, containing the control CDR3 data, with one row per clonotype
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
      stop(paste0('dataToTest has ', sum(is.na(dataToTest[[fn]])), ' NA values for: ', fn))
    }

    if (any(is.na(controlData[[fn]]))) {
      stop(paste0('controlData has ', sum(is.na(controlData[[fn]])),' NA values for: ', fn))
    }
  }

  if (any(is.na(dataToTest[[cloneFieldName]]))) {
    stop('Found NA values for clonotype in dataToTest')
  }

  if (nrow(controlData) == 0) {
    stop('There are no rows in controlData')
  }

  if (any(is.na(controlData[[cloneFieldName]]))) {
    stop('Found NA values for clonotype in controlData')
  }

  uniqueClones <- sort(unique(c(dataToTest[[cloneFieldName]], controlData[[cloneFieldName]])))
  print(paste0('Total clonotypes: ', length(uniqueClones)))

  df <- rbind(
    dataToTest %>% select(all_of(c(cloneFieldName, groupingField, contrastField, cloneSizeField))),
    controlData %>% select(all_of(c(cloneFieldName, groupingField, contrastField, cloneSizeField)))
  )

  # This converts the dataframe with one row/clonotype into one row/cell
  df <- as.data.frame(lapply(df, rep, df[[cloneSizeField]]))

  df <- df %>% rename(
    'Clonotype' = cloneFieldName,
    'GroupName' = groupingField
  )

  ctlGroup <- unique(controlData[[groupingField]])
  if (all(is.na(ctlGroup)) || all(is.null(ctlGroup))) {
    stop(paste0('controlData did not have any values for groupingField: ', groupingField))
  }

  if (length(ctlGroup) > 1) {
    stop('controlData can only have one value for groupingField')
  }

  # Set control group to intercept:
  df$GroupName <- as.factor(df$GroupName)
  if (all(is.na(levels(df$GroupName))) || length(levels(df$GroupName)) == 0) {
    stop(paste0('The grouping variable did not have any values: ', groupingField))
  }

  if (!ctlGroup %in% levels(df$GroupName)) {
    df$GroupName <- forcats::fct_expand(df$GroupName, as.character(ctlGroup), after = 0)
  } else {
    df$GroupName <- forcats::fct_relevel(df$GroupName, as.character(ctlGroup), after = 0)
  }

  if (as.character(ctlGroup) != levels(df$GroupName)[1]) {
    stop(paste0('intercept not properly set! levels: ', paste0(levels(df$GroupName), collapse = ','), ', expected: ', ctlGroup))
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
        return(data.frame(Clonotype = cdr3,
                          GroupName = dat$GroupName,
                          coefficients = NA,
                          p_val = NA,
                          error = paste0("Too few cells: ", nrow(dat))
        ))
      }

      if (length(unique(dat$GroupName)) == 1) {
        return(data.frame(Clonotype = cdr3,
                          GroupName = dat$GroupName,
                          coefficients = NA,
                          p_val = NA,
                          error = "Single Group"
        ))
      }

      model <- logistf::logistf(IsActive ~ GroupName, family = "binomial", data = dat)
      groupNames <- gsub(names(stats::coef(model)), pattern = 'GroupName', replacement = '')
      groupNames[1] <- 'CONTROL'

      return(data.frame(Clonotype = cdr3,
                                  GroupName = groupNames,
                                  coefficients = stats::coef(model),
                                  p_val = model$prob,
                                  error = NA
      ))
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
      results$FDR[results$GroupName == groupName] <- stats::p.adjust(results$p_val[results$GroupName == groupName], method = "fdr")
  }

  results <- results %>%
    select(Clonotype, GroupName, coefficients, p_val, FDR, error) %>%
    unique()

  # Restore original datatype
  if (is.integer(dataToTest[[groupingField]])) {
    results$GroupName <- as.integer(as.character(results$GroupName))
    if (any(is.na(results$GroupName))) {
      stop('NAs introduced into GroupName after conversion')
    }
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
  if (nrow(dat) == 0) {
    return(NULL)
  }

  dataWithPVal <- NULL
  for (ctlId in unique(dat$NoStimId)) {
    inputData <- dat %>%
      filter(NoStimId == ctlId & cDNA_ID != ctlId) %>%
      filter(Clonotype != 'Low Freq')

    dataToTest <- inputData %>%
      select(SubjectId, cDNA_ID, Clonotype, IsActive, TotalCellsForCloneAndState)

    controlData <- dat %>%
      filter(cDNA_ID == ctlId) %>%
      select(SubjectId, cDNA_ID, Clonotype, IsActive, TotalCellsForCloneAndState)

    if (nrow(dataToTest) == 0) {
      print(paste0('No experimental rows, skipping control ID: ', ctlId))
      next
    }

    if (nrow(controlData) == 0) {
      print(paste0('There are no rows in controlData, skipping: ' , ctlId))
      next
    }

    if (length(unique(dataToTest$SubjectId)) > 1){
      stop(paste0('More than one subject for: ', ctlId))
    }

    results <- CalculateClonotypeEnrichment(dataToTest, controlData, showProgress = showProgress)

    dataWithPVal <- rbind(dataWithPVal, inputData %>% filter(IsActive) %>% left_join(results, by = c('cDNA_ID', 'Clonotype')))
  }

  if (all(is.null(dataWithPVal))) {
    return(NULL)
  }

  # Set a floor:
  maxNonInfiniteTransformedFDR <- max(-log10(dataWithPVal$FDR[!is.infinite(-log10(dataWithPVal$FDR))]), na.rm = TRUE)
  dataWithPVal$FDR[dataWithPVal$FDR == 0] <- maxNonInfiniteTransformedFDR + 0.0001
  P1 <- ggplot(dataWithPVal %>% filter(!is.na(FDR) & !IsControlSample), aes(x = coefficients, y = -log10(FDR), color = Stim, label = Clonotype)) +
    geom_point() +
    ggrepel::geom_label_repel(show.legend = FALSE, size = 3) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    egg::theme_article() +
    xlab("Log Odds Ratio") +
    ylab("-log10(FDR)") +
    ggtitle('Clonotype Enrichment')

  print(P1)

  return(dataWithPVal)
}

.GenerateTcrQcPlots <- function(dat) {
  subjectId <- paste0(dat %>% select(SubjectId) %>% unique() %>% as.character(), collapse = ',')
  toPlot <- dat %>%
    filter(IsActive) %>%
    mutate(Filter = {ifelse(is.na(Filter), yes = 'Pass', no = Filter)}) %>%
    mutate(Fraction = TotalCellsForCloneAndState / TotalCellsForSample)

  if (nrow(toPlot) == 0) {
    filterPlot1 <- patchwork::plot_spacer()
  } else {
    filterPlot1 <- ggplot(toPlot, aes(x = Stim, y = Fraction, fill = Filter)) +
      geom_col(color = 'black') +
      facet_grid(. ~ SampleDate, scales = 'free', space = 'free_x') +
      scale_y_continuous(labels = scales::percent) +
      egg::theme_article() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        y = '% Activated',
        x = ''
      ) +
      ggtitle('Filter Summary')
  }

  toPlot <- dat %>%
    filter(IsActive & !IsControlSample) %>%
    mutate(Filter = {ifelse(is.na(Filter), yes = 'Pass', no = Filter)})

  if (nrow(toPlot) == 0) {
    filterPlot2 <- patchwork::plot_spacer()
  } else {
    filterPlot2 <- ggplot(toPlot, aes(x = Stim, y = FractionOfCloneWithState)) +
      geom_boxplot(color = 'black', outlier.shape = NA) +
      geom_jitter(aes(color = Filter, size = TotalCellsForCloneAndState)) +
      facet_grid(IsActiveLabel ~ SampleDate, scales = 'free', space = 'free') +
      scale_y_continuous(labels = scales::percent) +
      egg::theme_article() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        y = '% Active (of clonotype)',
        x = '',
        size = '# Cells',
        color = 'Filter',
        title = '% Activated By Clone and Stim'
      )
  }

  toPlot <- dat %>%
    filter(IsActive & !IsControlSample) %>%
    mutate(Filter = {ifelse(is.na(Filter), yes = 'Pass', no = Filter)})

  if (nrow(toPlot) == 0) {
    filterPlot3 <- patchwork::plot_spacer()
  } else {
    filterPlot3 <- ggplot(toPlot, aes(x = FractionOfCloneWithStateInSample, y = FractionOfCloneWithState)) +
      geom_jitter(aes(color = Filter, size = TotalCellsForCloneAndState)) +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      egg::theme_article() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        y = '% Active (of clonotype)',
        x = '% Activated (of total cells)',
        size = '# Cells',
        color = 'Filter',
        title = '% Active/Clone vs. % Activation'
      )
  }

  toPlot <- dat %>%
    filter(IsActive & !IsControlSample) %>%
    mutate(Filter = {ifelse(is.na(Filter), yes = 'Pass', no = Filter)})

  if (nrow(toPlot) == 0) {
    filterPlot4 <- patchwork::plot_spacer()
  } else {
    filterPlot4 <- ggplot(toPlot, aes(x = TotalCellsForClone/TotalCellsForSample, y = FractionOfCloneWithState)) +
      geom_jitter(aes(color = Filter, size = TotalCellsForCloneAndState)) +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      egg::theme_article() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        y = '% Active (of clonotype)',
        x = '% of Total Cells',
        size = '# Cells',
        color = 'Filter',
        title = '% Active/Clone vs. Total Clone Size'
      )
  }

  P <- ((filterPlot1 + filterPlot2) / (filterPlot3 + filterPlot4)) +
    patchwork::plot_layout(guides = 'collect') +
    patchwork::plot_annotation(title = paste0(subjectId, ': Filter QC'))

  return(P)
}

.getNumAntigensFieldName <- function(fieldPrefix){
  return(ifelse(is.null(fieldPrefix), yes = 'NumAntigens', no = paste0(fieldPrefix, 'NumAntigens')))
}

.getAntigensFieldName <- function(fieldPrefix){
  return(ifelse(is.null(fieldPrefix), yes = 'Antigens', no = paste0(fieldPrefix, 'Antigens')))
}

#' @title ApplyKnownClonotypicData
#' @description This will query the clone_responses table and append a column tagging each cell for matching antigens (based on clonotype)
#'
#' @param seuratObj The seurat object
#' @param antigenInclusionList If provided, only antigens from this list will be considered
#' @param antigenExclusionList If provided, antigens on this list will be omitted
#' @param minActivationFrequency If provided, only responses with activationFrequency (of the parent population) will be included
#' @param minFractionCloneActivated If provided, only responses where fractionCloneActivated is above this value will be will be included
#' @param fieldPrefix If provided, this will be appended to the beginning of the output field names
#' @param doNotPruneUsingCognateChain By default, matching is performed using the primary chain (as defined in the DB). If this is provided, this additional step is skipped
#' @export
ApplyKnownClonotypicData <- function(seuratObj, antigenInclusionList = NULL, antigenExclusionList = NULL, minActivationFrequency = 0, minFractionCloneActivated = 0, fieldPrefix = NULL, doNotPruneUsingCognateChain = FALSE) {
  numAntigensFieldName <- .getNumAntigensFieldName(fieldPrefix)
  antigensFieldName <- .getAntigensFieldName(fieldPrefix)

  subjectIds <- sort(unique(seuratObj$SubjectId))
  responseData <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="tcrdb",
    queryName="clone_responses",
    colSelect="cDNA_ID/sortId/sampleId/subjectId,cDNA_ID/sortId/sampleId/stim,chain,clonotype,totalclonesize,fractioncloneactivated,activationfrequency,clonename,cognatecdr3s",
    colFilter=makeFilter(
      c("cDNA_ID/sortId/sampleId/subjectId", "IN", paste0(subjectIds, collapse = ';')),
      c('clonotype', "NEQ", "No TCR"),
      c('cDNA_ID/status', 'DOES_NOT_CONTAIN', 'Failed'),
      c('status', "NOT_EQUAL_OR_MISSING", "Below Threshold")
    ),
    colNameOpt="rname"
  )

  names(responseData) <- c('SubjectId', 'Stim', 'Chain', 'Clonotype', 'totalclonesize', 'fractioncloneactivated', 'activationfrequency', 'CloneName', 'cognatecdr3s')

  if (nrow(responseData) == 0) {
    print('No matching clones found in DB, skipping')
    seuratObj[[numAntigensFieldName]] <- 0
    seuratObj[[antigensFieldName]] <- NA
    return(seuratObj)
  }

  return(.ApplyKnownClonotypicData(seuratObj,
                                   responseData = responseData,
                                   antigenInclusionList = antigenInclusionList,
                                   antigenExclusionList = antigenExclusionList,
                                   minActivationFrequency = minActivationFrequency,
                                   minFractionCloneActivated = minFractionCloneActivated,
                                   fieldPrefix = fieldPrefix,
                                   doNotPruneUsingCognateChain = doNotPruneUsingCognateChain
  ))
}

# responseData should be a data.frame with the columns: c('SubjectId', 'Stim', 'Chain', 'Clonotype', 'totalclonesize', 'fractioncloneactivated', 'activationfrequency')
.ApplyKnownClonotypicData <- function(seuratObj, responseData, antigenInclusionList = NULL, antigenExclusionList = NULL, minActivationFrequency = 0, minFractionCloneActivated = 0, fieldPrefix = NULL, doNotPruneUsingCognateChain = FALSE) {
  numAntigensFieldName <- .getNumAntigensFieldName(fieldPrefix)
  antigensFieldName <- .getAntigensFieldName(fieldPrefix)

  if (!all(is.null(antigenInclusionList))) {
    responseData <- responseData |>
      filter( Stim %in% antigenInclusionList)

    if (nrow(responseData) == 0) {
      print('No matching clones after applying inclusionList, skipping')
      seuratObj[[numAntigensFieldName]] <- 0
      seuratObj[[antigensFieldName]] <- NA
      return(seuratObj)
    }
  }

  if (!all(is.null(antigenExclusionList))) {
    responseData <- responseData |>
      filter(! Stim %in% antigenExclusionList)

    if (nrow(responseData) == 0) {
      print('No matching clones after applying exclusionList, skipping')
      seuratObj[[numAntigensFieldName]] <- 0
      seuratObj[[antigensFieldName]] <- NA
      return(seuratObj)
    }
  }

  if (!is.na(minActivationFrequency) && minActivationFrequency > 0) {
    responseData <- responseData |>
      filter(activationfrequency > minActivationFrequency)

    if (nrow(responseData) == 0) {
      print('No matching clones after applying minActivationFrequency, skipping')
      seuratObj[[numAntigensFieldName]] <- 0
      seuratObj[[antigensFieldName]] <- NA
      return(seuratObj)
    }
  }

  if (!is.na(minFractionCloneActivated) && minFractionCloneActivated > 0) {
    responseData <- responseData |>
      filter(fractioncloneactivated > minFractionCloneActivated)

    if (nrow(responseData) == 0) {
      print('No matching clones after applying minFractionCloneActivated, skipping')
      seuratObj[[numAntigensFieldName]] <- 0
      seuratObj[[antigensFieldName]] <- NA
      return(seuratObj)
    }
  }

  responseData <- responseData %>%
    mutate(
      IsNoStim = grepl(Stim, pattern = '^NoStim')
    ) %>%
    group_by(SubjectId, Clonotype) %>%
    summarize(
      Chain = paste0(sort(unique(Chain)), collapse = ','),
      CloneName = paste0(sort(unique(CloneName)), collapse = ','),
      Antigens = paste0(sort(unique(Stim)), collapse = ','),
      HasNoStim = sum(IsNoStim)>0,
      MaxTotalCloneSize = max(totalclonesize),
      MaxFractionCloneActivated = max(fractioncloneactivated),
      MaxActivationFrequency = max(activationfrequency),
      cognatecdr3s = paste0(sort(unique(cognatecdr3s)), collapse = ',')
    ) %>%
    as.data.frame() %>%
    mutate(
      DetectedAsSingleCDR3 = !grepl(Clonotype, pattern = ','),
      ClonotypesUsedForClonotypeMatch = paste0(Chain, ':', Clonotype)
    ) %>%
    tidyr::separate_longer_delim(Clonotype, delim = ',')

  if (nrow(responseData) == 0) {
    print('No matching clones, skipping')
    seuratObj[[numAntigensFieldName]] <- 0
    seuratObj[[antigensFieldName]] <- NA
    return(seuratObj)
  }

  toAppend <- NULL
  subjectIds <- sort(unique(seuratObj$SubjectId))
  for (subjectId in subjectIds) {
    responseDataForSubject <- responseData %>%
      filter(!is.na(SubjectId) & SubjectId == subjectId)
    if (nrow(responseDataForSubject) == 0) {
      next
    }

    for (idx in seq_len(nrow(responseDataForSubject))) {
      clonotype <- responseDataForSubject$Clonotype[idx]
      chain  <- NA
      if (any(grepl(x = clonotype, pattern = ':'))) {
        tokens <- unlist(strsplit(x = clonotype, split = ':'))
        chain <- tokens[1]
        clonotype <- tokens[2]
      } else {
        chain <- responseDataForSubject$Chain[idx]
        if (is.na(chain) || is.null(chain)) {
          stop(paste0('Clone lacks chain info: ', clonotype))
        }
      }

      cognateChain <- .GetCognateChain(chain)
      if (! cognateChain %in% colnames(seuratObj@meta.data)) {
        stop(paste0('Missing column: ', cognateChain))
      }

      cognateCDR3s <- responseDataForSubject$cognatecdr3s[idx]
      if (!is.na(cognateCDR3s) && !is.null(cognateCDR3s) && cognateCDR3s != '') {
        cognateCDR3s <- unlist(strsplit(cognateCDR3s, split = ','))
        cognateCDR3s <- unique(unlist(sapply(cognateCDR3s, function(x){
          # Strip frequency information:
          x <- unlist(strsplit(x, split = ':'))[1]
          return(x)
        })))
      }

      # Ensure consistent value when empty:
      if (all(is.na(cognateCDR3s)) || all(is.null(cognateCDR3s)) || all(cognateCDR3s == '')) {
        cognateCDR3s <- NULL
      }

      sel <- !is.na(seuratObj$SubjectId) & seuratObj$SubjectId == subjectId & grepl(pattern = paste0("(?:^|,)", clonotype, "(?:$|,)"), x = seuratObj[[chain]][[1]])
      if (any(sel)) {
        toAdd <- responseDataForSubject[rep(idx, sum(sel)),]
        toAdd$ChainsUsedForClonotypeMatch <- chain
        toAdd$ClonotypesUsedForClonotypeMatch <- responseDataForSubject$ClonotypesUsedForClonotypeMatch[idx]
        toAdd$CellBarcode <- rownames(seuratObj@meta.data)[sel]

        if (!doNotPruneUsingCognateChain && !all(is.null(cognateCDR3s))) {
          toRetain <- unlist(sapply(seuratObj@meta.data[[cognateChain]][sel], function(x){
            if (is.na(x) || x == '') {
              return(TRUE)
            }

            x <- unlist(strsplit(x, split = ','))

            return(length(intersect(x, cognateCDR3s)) > 0)
          }))

          if (sum(toRetain) != nrow(toAdd)) {
            toDrop <- nrow(toAdd) - sum(toRetain)
            print(paste0('Dropping ', toDrop, ' hits out of ', nrow(toAdd), ' because of mismatched second chain'))
            toAdd <- toAdd[toRetain,]
          }
        }

        toAppend <- rbind(toAppend, toAdd)
      }
    }
  }

  if (any(duplicated(toAppend$CellBarcode))) {
    toAppend <- toAppend %>%
      group_by(CellBarcode) %>%
      summarize(
        Antigens = paste0(sort(unique(Antigens)), collapse = ','),
        CloneName = paste0(sort(unique(CloneName)), collapse = ','),
        ChainsUsedForClonotypeMatch = paste0(sort(unique(ChainsUsedForClonotypeMatch)), collapse = ','),
        ClonotypesUsedForClonotypeMatch = paste0(unique(sort(ClonotypesUsedForClonotypeMatch)), collapse = ';'),
        HasNoStim = max(HasNoStim),
        MaxTotalCloneSize = max(MaxTotalCloneSize),
        MaxFractionCloneActivated = max(MaxFractionCloneActivated),
        MaxActivationFrequency = max(MaxActivationFrequency),
        DetectedAsSingleCDR3 = max(DetectedAsSingleCDR3)
      ) %>%
      as.data.frame()

    toAppend$Antigens <- unlist(sapply(toAppend$Antigens, function(x){
      if (is.na(x)) {
        return(NA)
      }

      x <- sort(unique(unlist(strsplit(x, split = ','))))

      return(paste0(x, collapse = ','))
    }))

    toAppend$CloneName <- unlist(sapply(toAppend$CloneName, function(x){
      if (is.na(x)) {
        return(NA)
      }

      x <- sort(unique(unlist(strsplit(x, split = ','))))

      return(paste0(x, collapse = ','))
    }))

    toAppend$ChainsUsedForClonotypeMatch <- unlist(sapply(toAppend$ChainsUsedForClonotypeMatch, function(x){
      if (is.na(x)) {
        return(NA)
      }

      x <- sort(unique(unlist(strsplit(x, split = ','))))

      return(paste0(x, collapse = ','))
    }))

    toAppend$ClonotypesUsedForClonotypeMatch <- unlist(sapply(toAppend$ClonotypesUsedForClonotypeMatch, function(x){
      if (is.na(x)) {
        return(NA)
      }

      x <- sort(unique(unlist(strsplit(x, split = ';'))))

      return(paste0(x, collapse = ';'))
    }))
  }

  rownames(toAppend) <- toAppend$CellBarcode
  toAppend$CellBarcode <- NULL
  toAppend$Clonotype <- NULL
  toAppend$Chain <- NULL

  toAppend$NumAntigens <- sapply(toAppend$Antigens, function(x){
    if (is.na(x) || x == '') {
      return(0)
    }

    x <- unique(unlist(strsplit(x, split = ',')))

    return(length(x))
  })
  toAppend$NumAntigens[is.na(toAppend$NumAntigens)] <- 0

  if (!is.null(fieldPrefix)) {
    names(toAppend) <- paste0(fieldPrefix, names(toAppend))
  }

  # Clear existing values:
  for (n in names(toAppend)) {
    if (n %in% names(seuratObj@meta.data)) {
      seuratObj[[n]] <- NULL
    }
  }
  seuratObj <- Seurat::AddMetaData(seuratObj, toAppend)

  # Always provide a value
  seuratObj[[numAntigensFieldName]][is.na(seuratObj[[numAntigensFieldName]])] <- 0

  if (length(names(seuratObj@reductions)) > 0) {
    print(DimPlot(seuratObj, group.by = antigensFieldName))
    print(DimPlot(seuratObj, group.by = numAntigensFieldName))
  }

  return(seuratObj)
}


#' @title IdentifyAndStoreActiveClonotypes
#' @description This is a very specialized method designed to score and summarize activated TCR clonotypes
#'
#' @param seuratObj The seurat object. This expects
#' @param chain The chain to summarize
#' @param method Either 'Cluster-Based' or 'sPLS'
#' @param storeStimLevelData If true, activation levels will be stord in tcrtb.stims
#' @param maxRatioToCombine Passed to GroupOverlappingClones
#' @param minEDS If provided, only cells with EDS>minEDS will be included
#' @export
#' @import Rlabkey
#' @import dplyr
IdentifyAndStoreActiveClonotypes <- function(seuratObj, chain = 'TRB', method = 'sPLS', storeStimLevelData = TRUE, maxRatioToCombine = 1.0, minEDS = 2) {
  allDataWithPVal <- .IdentifyActiveClonotypes(seuratObj, chain = chain, method = method, maxRatioToCombine = maxRatioToCombine, minEDS = minEDS)
  allDataWithPVal$chain <- chain

  # Calculate/store frequencies for clones that responded in at least one sample:
  allDataWithPVal$Status <- NA
  allClones <- unique(allDataWithPVal$Clonotype)
  allClones <- allClones[allClones != 'No TCR']
  for (subjectId in sort(unique(allDataWithPVal$SubjectId))) {
    datForSubject <- allDataWithPVal %>% filter(SubjectId == subjectId)
    for (cdnaId in sort(unique(datForSubject$cDNA_ID))) {
      datForStim <- datForSubject %>% filter(cDNA_ID == cdnaId)
      missingClonotypes <- allClones[! allClones %in% datForStim$Clonotype ]
      if (length(missingClonotypes) == 0) {
        next
      }

      chainWithSegmentsField <- paste0(chain, '_Segments')
      vField <- paste0(chain, '_V')
      jField <- paste0(chain, '_J')
      toAppend <- seuratObj@meta.data %>%
        filter(cDNA_ID == cdnaId) %>%
        filter(dplyr::if_all(dplyr::all_of(c(chain)), ~ !is.na(.x))) %>%
        group_by(cDNA_ID) %>%
        mutate(TotalCellsForSample = n()) %>%
        ungroup() %>%
        select(all_of(c('cDNA_ID', chain, vField, jField, chainWithSegmentsField))) %>%
        rename(c(
          Clonotype = !!chain,
          cdr3WithSegments = !!chainWithSegmentsField,
          V_Gene = !!vField,
          J_Gene = !!jField
        )) %>%
        as.data.frame() %>%
        filter(Clonotype %in% missingClonotypes) %>%
        group_by(cDNA_ID, Clonotype, V_Gene, J_Gene, cdr3WithSegments) %>%
        summarize(TotalCellsForClone = n()) %>%
        as.data.frame() %>%
        mutate(
          method = method,
          TotalCellsForCloneAndState = 0,
          FractionOfCloneWithState = 0,
          chain = chain,
          IsActive = TRUE
        )

      stillNeeded <- missingClonotypes[! missingClonotypes %in% toAppend$Clonotype]
      if (length(stillNeeded) > 0) {
        toAppend2 <- seuratObj@meta.data %>%
          select(all_of(c(chain, vField, jField, chainWithSegmentsField))) %>%
          rename(c(
            Clonotype = !!chain,
            cdr3WithSegments = !!chainWithSegmentsField,
            V_Gene = !!vField,
            J_Gene = !!jField
          )) %>%
          as.data.frame() %>%
          filter(Clonotype %in% stillNeeded) %>%
          unique() %>%
          mutate(
            method = method,
            TotalCellsForCloneAndState = 0,
            FractionOfCloneWithState = 0,
            chain = chain,
            cDNA_ID = cdnaId,
            IsActive = TRUE
          )

        if (nrow(toAppend2) > 0) {
          toAppend <- plyr::rbind.fill(toAppend, toAppend2)
        }
      }

      if (nrow(toAppend) > 0) {
        print(paste0('cDNA ', cdnaId, ': adding ', nrow(toAppend), ' placeholder records for clonotypes without activation'))
        toAppend$Status <- 'Below Threshold'

        allDataWithPVal <- plyr::rbind.fill(allDataWithPVal, toAppend)
      }
    }
  }

  # Add cognate chains:
  allDataWithPVal <- .AddCognateChains(chain, seuratObj, allDataWithPVal)

  .UpdateTcrStimDb(allDataWithPVal, chain = chain, methodName = method, storeStimLevelData = storeStimLevelData, allCDNA_IDs = unique(seuratObj$cDNA_ID))
}

.IdentifyActivatedCluster <- function(dat, resolutions = c('0.2', '0.4', '0.6', '0.8', '1.2')) {
  foundHit <- FALSE
  for (threshold in c(0.3, 0.25)) {
    for (res in resolutions) {
      fieldName <- paste0('ClusterNames_', res)
      if (!fieldName %in% names(dat)) {
        stop(paste0('Missing: ', fieldName))
      }

      if (!'TandNK_ActivationCore_UCell' %in% names(dat)) {
        stop('Missing TandNK_ActivationCore_UCell')
      }

      if (!'TandNK_Activation_UCell' %in% names(dat)) {
        stop('Missing TandNK_Activation_UCell')
      }

      if (!'TandNK_Activation3_UCell' %in% names(dat)) {
        stop('Missing TandNK_Activation3_UCell')
      }

      df <- dat %>%
        group_by(across(all_of(fieldName))) %>%
        summarize(
          TandNK_ActivationCore_UCell = mean(TandNK_ActivationCore_UCell),
          TandNK_Activation_UCell = mean(TandNK_Activation_UCell),
          TandNK_Activation3_UCell = mean(TandNK_Activation3_UCell)
        )

      df$CombinedScore <- pmax(df$TandNK_Activation_UCell, df$TandNK_Activation3_UCell, df$TandNK_ActivationCore_UCell)

      df[[fieldName]] <- naturalsort::naturalfactor(df[[fieldName]])
      P1 <- ggplot(df, aes_string(x = fieldName, y = 'TandNK_Activation3_UCell'), fill = fieldName) +
        geom_col(color = 'black') +
        egg::theme_article() +
        geom_hline(yintercept = threshold, color = 'red') +
        ggtitle('TandNK_Activation3_UCell') +
        NoLegend()

      P2 <- ggplot(df, aes_string(x = fieldName, y = 'TandNK_Activation_UCell'), fill = fieldName) +
        geom_col(color = 'black') +
        egg::theme_article() +
        geom_hline(yintercept = threshold, color = 'red') +
        ggtitle('TandNK_Activation_UCell') +
        NoLegend()

      P4 <- ggplot(df, aes_string(x = fieldName, y = 'TandNK_ActivationCore_UCell'), fill = fieldName) +
        geom_col(color = 'black') +
        egg::theme_article() +
        geom_hline(yintercept = threshold, color = 'red') +
        ggtitle('TandNK_ActivationCore_UCell') +
        NoLegend()

      P3 <- ggplot(df, aes_string(x = fieldName, y = 'CombinedScore'), fill = fieldName) +
        geom_col(color = 'black') +
        egg::theme_article() +
        geom_hline(yintercept = threshold, color = 'red') +
        ggtitle('CombinedScore') +
        NoLegend()

      P <- P1 + P2 + P4 + P3 + patchwork::plot_annotation(title = res)
      print(P)

      if (max(df$CombinedScore) < threshold) {
        next
      }

      return(data.frame(resolution = fieldName, cluster = as.character(df[[fieldName]][which(df$CombinedScore == max(df$CombinedScore))]), value = max(df$CombinedScore), threshold = threshold))
    }

    if (foundHit) {
      return(NULL)
    }
  }
}

.GetCognateChain <- function(chain) {
  return(switch(chain,
         "TRA" = "TRB",
         "TRB" = "TRA",
         "TRG" = "TRD",
         "TRD" = "TRG",
         stop(paste0('Uknown chain: ', chain))
  ))
}

.IdentifyActiveClonotypes <- function(seuratObj, chain = 'TRB', method = 'sPLS', maxRatioToCombine = 1.0, minOddsRatio = 0.5, minEDS = 2) {
  if (method == 'Cluster-Based') {
    activatedCluster <- .IdentifyActivatedCluster(dat)
    if (all(is.null(activatedCluster))) {
      print('Unable to find activated cluster')
      return(NULL)
    }

    print(paste0('Activated cluster: ', activatedCluster$resolution, ', ', activatedCluster$cluster))
    seuratObj$IsActive <- seuratObj@meta.data[[activatedCluster$resolution]] == activatedCluster$cluster
    print(paste0('Active cluster: ', activatedCluster$cluster, '. Total cells: ', sum(seuratObj$IsActive)))
  } else if (method == 'sPLS') {
    # This is create by RIRA::PredictTcellActivation()
    if (! 'GeneralizedTCR_sPLSDA_ConsensusClass' %in% names(seuratObj@meta.data)) {
        stop('Missing field: GeneralizedTCR_sPLSDA_ConsensusClass')
    }

    if (! 'Is_TCR_Stimulated' %in% names(seuratObj@meta.data)) {
      stop('Missing field: Is_TCR_Stimulated. This should be created by RIRA::PredictTcellActivation()')
    }

    if (any(!is.na(seuratObj@meta.data[[chain]]) & is.na(seuratObj$Is_TCR_Stimulated))) {
      stop('There were NA values for Is_TCR_Stimulated that contained TCR data')
    }

    seuratObj$IsActive <- seuratObj$Is_TCR_Stimulated
    print(paste0('Total active cells: ', sum(seuratObj$IsActive)))
  } else {
    stop(paste0('Unknown method: ', method))
  }

  allDataWithPVal <- NULL
  for (subjectId in sort(unique(seuratObj$SubjectId))) {
    dat <- seuratObj@meta.data %>%
      filter(SubjectId == subjectId)

    if (!is.integer(dat$cDNA_ID)) {
      print('Converting cDNA_ID to an integer')
      dat$cDNA_ID <- as.character(dat$cDNA_ID)
      converted <- as.integer(dat$cDNA_ID)
      if (any(is.na(converted))) {
        stop('Non-numeric cDNA_IDs found: ', paste0(unique(dat$cDNA_ID[is.na(converted)]), collapse = ','))
      }

      dat$cDNA_ID <- converted
    }

    dat <- PrepareTcrData(dat, subjectId = subjectId, minEDS = minEDS, retainRowsWithoutCDR3 = TRUE, chain = chain, enforceAllDataPresent = FALSE)
    dat$MethodString <- method

    dat <- GroupOverlappingClones(dat, maxRatioToCombine = maxRatioToCombine, dataMask = dat$IsActive, groupingFields = c('cDNA_ID', 'SubjectId', 'SampleDate', 'Stim', 'NoStimId', 'IsControlSample', 'IsActiveLabel', 'MethodString', 'AssayType'))
    dat <- ApplyCloneFilters(dat, minCellsPerClone = 2, minFractionOfCloneActive = 0.025, minFoldChangeAboveNoStim = NA)
    dat$IsFiltered <- ifelse(is.na(dat$Filter), yes = 'No', no = 'Yes')

    filterPlots <- .GenerateTcrQcPlots(dat)
    print(filterPlots)

    subjectId <- paste0(sort(unique(dat$SubjectId)), collapse = ',')
    P1 <- GenerateTcrPlot(dat, xFacetField = 'SampleDate', dropInactive = FALSE, patternField = 'IsFiltered', plotTitle = paste0(subjectId, ": EDS > ", minEDS, ", Unfiltered"), groupLowFreq = TRUE)

    P1 <- P1 + patchwork::plot_annotation(title = paste0(subjectId, ': Unfiltered Clonotypes and ICS'))
    print(P1)

    passingRows <- dat %>%
      filter(IsFiltered == 'No')
    if (nrow(passingRows) == 0) {
      print('No rows passed filters, skipping')
      next
    }

    dataWithPVal <- AppendClonotypeEnrichmentPVals(passingRows)
    if (all(is.null(dataWithPVal))){
      next
    }

    dataWithPVal$FailedEnrichment <- !is.na(dataWithPVal$coefficients) & dataWithPVal$coefficients < minOddsRatio

    tryCatch({
      passingClones <- GenerateTcrPlot(dataWithPVal, xFacetField = 'SampleDate', dropInactive = TRUE, patternField = 'FailedEnrichment', plotTitle = paste0(subjectId, ": EDS > ", minEDS, ", Passing Enrichment"), groupLowFreq = FALSE)
      if (!all(is.null(passingClones))){
        passingClones <- passingClones +
          geom_hline(yintercept = 0.005, linetype = 'dotted', colour = 'red', linewidth = 1) +
          theme(
            legend.position = ifelse(n_distinct(dataWithPVal$Clonotype) > 15, yes = 'none', no = 'right'),
            legend.text = element_text(size = rel(0.5)),
            legend.title = element_text(size = rel(0.75))
          ) +
          labs(pattern = 'Failed Enrichment?')

        # Reset the fill colors:
        if (n_distinct(dataWithPVal$Clonotype) < 10) {
          passingClones <- passingClones + scale_fill_discrete()
        }

        print(passingClones)

      }
    }, error = function(e){
      warning('Error generating passingClones plot')
      print(conditionMessage(e))
      traceback()
    })

    allDataWithPVal <- rbind(allDataWithPVal, dataWithPVal)
  }

  return(allDataWithPVal)
}

.AddCognateChains <- function(chain, seuratObj, allDataWithPVal, minFraction = 0.05) {
  cognateChain <- .GetCognateChain(chain)

  chainData <- data.frame(SourceChain = allDataWithPVal$Clonotype, ToJoin = allDataWithPVal$Clonotype) %>%
    group_by(SourceChain, ToJoin) %>%
    summarize(CellForSource = n()) %>%
    as.data.frame() %>%
    tidyr::separate_longer_delim(cols = ToJoin, delim = ',') %>%
    filter(!is.na(SourceChain)) %>%
    unique()

  cognateData <- seuratObj@meta.data %>%
    select(all_of(c(chain, cognateChain))) %>%
    rename(
      'TargetChain' = !!chain,
      'CognateChain' = !!cognateChain
    ) %>%
    filter(!is.na(TargetChain)) %>%
    group_by(TargetChain, CognateChain) %>%
    summarize(CellForCognate = n()) %>%
    as.data.frame() %>%
    mutate(ToJoin = TargetChain) %>%
    tidyr::separate_longer_delim(cols = ToJoin, delim = ',')

  joinedData <- chainData %>%
    left_join(cognateData, by = 'ToJoin', relationship = 'many-to-many') %>%
    group_by(SourceChain) %>%
    mutate(TotalForSource = sum(CellForCognate[!is.na(CognateChain)])) %>%
    group_by(SourceChain, CognateChain, TotalForSource) %>%
    summarize(Total = sum(CellForCognate[!is.na(CognateChain)])) %>%
    as.data.frame() %>%
    mutate(Fraction = Total / TotalForSource) %>%
    mutate(Weight = 1/(stringr::str_count(CognateChain, ",")+1)) %>%
    mutate(WeightedTotal = Total * Weight) %>%
    tidyr::separate_longer_delim(cols = CognateChain, delim = ',') %>%
    group_by(SourceChain, CognateChain, TotalForSource) %>%
    summarise(WeightedTotal = sum(WeightedTotal)) %>%
    as.data.frame() %>%
    mutate(Fraction = WeightedTotal / TotalForSource) %>%
    filter(Fraction > minFraction) %>%
    group_by(SourceChain) %>%
    mutate(numCognates = n()) %>%
    mutate(ChainAndWeight = paste0(CognateChain, ':', round(Fraction, 3))) %>%
    group_by(SourceChain) %>%
    summarize(cognateCdr3s = paste0(ChainAndWeight, collapse = ',')) %>%
    as.data.frame()

  allDataWithPVal <- allDataWithPVal %>%
    left_join(joinedData, by = c('Clonotype' = 'SourceChain'))

  return(allDataWithPVal)
}

.UpdateTcrStimDb <- function(allDataWithPVal, chain, methodName = NULL, storeStimLevelData = TRUE, allCDNA_IDs = NULL) {
  if (all(is.null(allDataWithPVal))) {
    print('No data, skipping import')

    if (!all(is.null(allCDNA_IDs))) {
      toDelete <- suppressWarnings(labkey.selectRows(
        baseUrl=.getBaseUrl(),
        folderPath=.getLabKeyDefaultFolder(),
        schemaName="tcrdb",
        queryName="clone_responses",
        colSelect="rowid,container",
        colFilter=makeFilter(
          c("cdna_id", "IN", paste0(unique(allCDNA_IDs), collapse = ';')),
          c('chain', 'EQUALS', chain)
        ),
        containerFilter=NULL,
        colNameOpt="rname"
      ))

      if (nrow(toDelete) > 0) {
        print(paste0('Deleting ', nrow(toDelete), ' rows in tcrdb.clone_responses'))
        deleted <- labkey.deleteRows(
          baseUrl=.getBaseUrl(),
          folderPath=.getLabKeyDefaultFolder(),
          schemaName="tcrdb",
          queryName="clone_responses",
          toDelete = toDelete
        )
      }
    }

    return()
  }

  cloneNames <- allDataWithPVal %>%
    group_by(SubjectId, Clonotype) %>%
    summarize() %>%
    as.data.frame() %>%
    group_by(SubjectId) %>%
    mutate(CloneId = row_number()) %>%
    mutate(CloneName = paste0(SubjectId, '-', chain, '-', CloneId)) %>%
    select(-CloneId)

  allDataWithPVal <- allDataWithPVal %>%
    left_join(cloneNames, by = c('SubjectId', 'Clonotype'))
  rm(cloneNames)

  if (!is.integer(allDataWithPVal$cDNA_ID)) {
    print('Converting cDNA_ID to an integer')
    allDataWithPVal$cDNA_ID <- as.character(allDataWithPVal$cDNA_ID)
    converted <- as.integer(allDataWithPVal$cDNA_ID)
    if (any(is.na(converted))) {
      stop('Non-numeric cDNA_IDs found: ', paste0(unique(allDataWithPVal$cDNA_ID[is.na(converted)]), collapse = ','))
    }
    
    allDataWithPVal$cDNA_ID <- converted
  }
  
  toUpdate <- allDataWithPVal %>%
    group_by(cDNA_ID) %>%
    filter(IsActive) %>%
    filter(is.na(FailedEnrichment) | !FailedEnrichment) %>%
    summarise(quantification = unique(FractionOfSampleWithState)*100, nClones = n()) %>%
    mutate(QuantificationMethod = methodName)

  allCDNA <- toUpdate$cDNA_ID
  if (!all(is.null(allCDNA_IDs))) {
    allCDNA <- sort(unique(c(allCDNA, allCDNA_IDs)))
  }

  containerInfo <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="tcrdb",
    queryName="stims",
    colSelect="rowid,cdna_id,container",
    colFilter=makeFilter(
      c("cdna_id/rowid", "IN", paste0(allCDNA, collapse = ';'))
    ),
    colNameOpt="rname"
  )

  toUpdate <- toUpdate %>%
    left_join(containerInfo, by = c('cDNA_ID' = 'cdna_id'))

  if (storeStimLevelData) {
    print(paste0('Updating ', nrow(toUpdate), ' rows in tcrdb.stims'))

    if (nrow(toUpdate) > 0) {
      updated <- labkey.updateRows(
        baseUrl=.getBaseUrl(),
        folderPath=.getLabKeyDefaultFolder(),
        schemaName="tcrdb",
        queryName="stims",
        toUpdate = toUpdate
      )
    } else {
      print('No updates needed for tcrdb.stims')
    }
  }

  # Now clone data:
  toInsertOrUpdate <- allDataWithPVal %>%
    filter(IsActive) %>%
    filter(is.na(FailedEnrichment) | !FailedEnrichment) %>%
    rename(
      activationFrequency = 'FractionOfCloneWithStateInSample',
      totalCells = 'TotalCellsForCloneAndState',
      vGene = 'V_Gene',
      jGene = 'J_Gene',
      cdr3WithSegments = 'cdr3WithSegments',
      totalCloneSize = 'TotalCellsForClone',
      fractionCloneActivated = 'FractionOfCloneWithState',
      totalCellsForSample = 'TotalCellsForSample',
      oddsRatio = 'coefficients',
      enrichmentFDR = 'FDR',
      comments = 'MethodString'
    ) %>%
    mutate(chain = chain)

  existingRows <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
    schemaName="tcrdb",
    queryName="clone_responses",
    colSelect="rowid,clonotype,cdna_id,container",
    colFilter=makeFilter(
      c("cdna_id", "IN", paste0(unique(allCDNA), collapse = ';')),
      c('chain', 'EQUALS', chain)
    ),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  # This seems to occur if there are zero rows:
  if (is.numeric(existingRows$clonotype)) {
    existingRows$clonotype <- as.character(existingRows$clonotype)
  }

  # And Clones/Responses:
  toDelete <- existingRows %>%
    left_join(toInsertOrUpdate, by = c('cdna_id' = 'cDNA_ID', 'clonotype' = 'Clonotype')) %>%
    filter(is.na(SubjectId))

  toInsertOrUpdate <- toInsertOrUpdate %>%
    left_join(existingRows, by = c('cDNA_ID' = 'cdna_id', 'Clonotype' = 'clonotype'))

  toInsert <- toInsertOrUpdate %>%
    filter(is.na(rowid)) %>%
    select(-rowid) %>%
    select(-container)

  if (nrow(toInsert) > 0) {
    existingLibraries <- labkey.selectRows(
      baseUrl=.getBaseUrl(),
      folderPath=.getLabKeyDefaultFolder(),
      schemaName="singlecell",
      queryName="cdna_libraries",
      colSelect="rowid,container",
      colFilter=makeFilter(
        c("rowid", "IN", paste0(unique(toUpdate$cDNA_ID), collapse = ';'))
      ),
      containerFilter=NULL,
      colNameOpt="rname"
    )
    toInsert <- toInsert %>% left_join(existingLibraries, by = c('cDNA_ID' = 'rowid'))
    for (fieldName in names(toInsert)) {
      if (all(is.na(toInsert[[fieldName]]))) {
        print(paste0('dropping all-NA field: ', fieldName))
        toInsert[[fieldName]] <- NULL
      }
    }

    for (fn in c('enrichmentFDR', 'oddsRatio')) {
      if (any(is.na(toInsert[[fn]]))) {
        toInsert[[fn]][is.na(toInsert[[fn]])] <- ''
      }
    }

    print(paste0('Inserting ', nrow(toInsert), ' rows in tcrdb.clone_responses'))
    added <- labkey.insertRows(
      baseUrl=.getBaseUrl(),
      folderPath=.getLabKeyDefaultFolder(),
      schemaName="tcrdb",
      queryName="clone_responses",
      toInsert = toInsert
    )
  }

  toUpdate <- toInsertOrUpdate %>% filter(!is.na(container))
  if (nrow(toUpdate) > 0) {
    for (fn in c('enrichmentFDR', 'oddsRatio')) {
      if (any(is.na(toUpdate[[fn]]))) {
        toUpdate[[fn]][is.na(toUpdate[[fn]])] <- ''
      }
    }

    print(paste0('Updating ', nrow(toUpdate), ' rows in tcrdb.clone_responses'))
    added <- labkey.updateRows(
      baseUrl=.getBaseUrl(),
      folderPath=.getLabKeyDefaultFolder(),
      schemaName="tcrdb",
      queryName="clone_responses",
      toUpdate = toUpdate
    )
  }

  if (nrow(toDelete) > 0) {
    toDelete <- toDelete %>%
      select(rowid, container)

    print(paste0('Deleting ', nrow(toDelete), ' rows in tcrdb.clone_responses'))
    deleted <- labkey.deleteRows(
      baseUrl=.getBaseUrl(),
      folderPath=.getLabKeyDefaultFolder(),
      schemaName="tcrdb",
      queryName="clone_responses",
      toDelete = toDelete
    )
  }
}
