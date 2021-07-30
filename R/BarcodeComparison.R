#' @include LabKeySettings.R
#' @import Rlabkey
#' @import utils

utils::globalVariables(
	names = c('dataset1', 'type1', 'type2', 'max_intersect', 'max_intersect_by_type', 'HTO', 'readsetid', 'sortid_hto', 'TRAJ', 'TRBJ', 'TRDJ', 'TRGJ', 'TRAC', 'TRBC', 'TRDC', 'TRGC', 'TRBD', 'TRDD'),
	package = 'Rdiscvr',
	add = TRUE
)

utils::globalVariables(
names = c('chain', 'cdr3', 'LabelCol', 'barcode', 'ChainCDR3s', 'TRA', 'TRB', 'TRD', 'TRG', 'TRAV', 'TRBV', 'TRDV', 'TRGV', 'CloneName'),
package = 'Rdiscvr',
add = TRUE
)

#' @title CompareCellBarcodeSets
#'
#' @description This iterates the provided workbooks, identifying every 10x cDNA library and the GEX/TCR/HTO libraries.  It generates summaries of the cell barcode intersect between them, which can help debug sample swaps.
#' @param workbooks A vector of workbook IDs
#' @param savePath The folder path to use for saving output files
#' @export
#' @importFrom dplyr %>% mutate group_by
CompareCellBarcodeSets <- function(workbooks, savePath = './') {
  summary <- .GenerateDataToCompareBarcodeSets(workbooks, savePath)

  htoBC <- list()
  gexBC <- list()
  tcrBC <- list()
  citeBC <- list()

  for (i in 1:nrow(summary)) {
    row <- summary[i,]
    name <- as.character(row$Name)

    if (!is.na(row$HTO_Top_BarcodesFile) & !is.na(row$GEX_CallsFile) & !is.na(row$TCR_CallsFile)){
      bc <- read.table(row$HTO_Top_BarcodesFile, header = T, sep = '\t')
      htoBC[[name]] <- bc$cellbarcode

      bc <- read.table(row$GEX_CallsFile, header = T, sep = '\t')$cellbarcode
      gexBC[[name]] <- bc

      bc <- read.table(row$TCR_CallsFile, header = T, sep = '\t')$cellbarcode
      tcrBC[[name]] <- bc

    } else {
      print(paste0('missing one or more files: ', name, ':'))
      if (is.na(row$HTO_Top_BarcodesFile)){
        print('HTO missing')
      }

      if (is.na(row$GEX_CallsFile)){
        print('GEX missing')
      }

      if (is.na(row$TCR_CallsFile)){
        print('TCR missing')
      }
    }
  }

  df <- data.frame(dataset1 = character(), type1 = character(), dataset2 = character(), type2 = character(), intersect = integer(), length1 = integer(), length2 = integer())

  df <- .ProcessSet(df, htoBC[1:length(gexBC)], gexBC, 'HTO', 'GEX')
  df <- .ProcessSet(df, htoBC[1:length(tcrBC)], tcrBC, 'HTO', 'TCR')

  df <- .ProcessSet(df, gexBC, htoBC[1:length(gexBC)], 'GEX', 'HTO')
  df <- .ProcessSet(df, gexBC, tcrBC, 'GEX', 'TCR')

  df <- .ProcessSet(df, tcrBC, htoBC[1:length(tcrBC)], 'TCR', 'HTO')
  df <- .ProcessSet(df, tcrBC, gexBC, 'TCR', 'GEX')

  df$fraction <- df$intersect / df$length1
  df <- df[order(df$dataset1, df$type1, df$type2, -df$intersect),]
  df <- df %>% group_by(dataset1, type1, type2) %>% mutate(max_intersect_by_type = max(intersect))
  df <- df %>% group_by(dataset1, type1) %>% mutate(max_intersect = max(intersect))

  write.table(df, file = file.path(savePath, 'cell_barcode_comparisons.txt'), sep = '\t', row.names = F, quote = F)

  self <- df[df$dataset1 == df$dataset2, c('dataset1', 'type1', 'type2', 'intersect', 'fraction')]
  names(self) <- c('dataset1', 'type1', 'type2', 'self_intersect', 'self_intersect_fraction')

  df2 <- df[df$intersect == df$max_intersect_by_type,]
  write.table(df2, file = 'top_intersect_by_type.txt', sep = '\t', row.names = F, quote = F)

  df3 <- df2[df2$dataset1 != df2$dataset2,]
  df3 <- merge(df3, self, by = c('dataset1', 'type1', 'type2'), all.x = T, all.y = F)
  #filter ties:
  df3 <- df3[df3$intersect != df3$self_intersect,]
  write.table(df3, file = file.path(savePath, 'conflicting_intersect.txt'), sep = '\t', row.names = F, quote = F)

  #Now look for instances where the raw HTO calls find more HTOs than either TCR or GEX
  dfSummary <- NA
  for (i in 1:nrow(summary)) {
    row <- summary[i,]
    name <- as.character(row$Name)

    if (is.na(row$HTO_Top_BarcodesFile)){
      print(paste0('No HTO call file for: ', name))
      next
    }

    df <- data.frame(HTO = character(), Type = character(), Count.HTO = integer(), Fraction.HTO = numeric(), Count.Compare = integer(), Fraction.Compare = numeric(), Difference = numeric())
    #GEX:
    dfG <- .CompareHtosByCall(name, row$HTO_Top_BarcodesFile, row$GEX_CallsFile, 'GEX')
    if (!all(is.na(dfG))){
      df <- rbind(df, dfG)
    }

    #TCR
    dfT <- .CompareHtosByCall(name, row$HTO_Top_BarcodesFile, row$TCR_CallsFile, 'TCR')
    if (!all(is.na(dfT))){
      df <- rbind(df, dfT)
    }

    if (all(is.na(dfSummary))) {
      dfSummary <- df
    } else {
      dfSummary <- rbind(dfSummary, df)
    }
  }

  dfSummary <- merge(dfSummary, summary[c('Name', 'ExpectedHTOs')], all.x = T, by.x = c('Dataset'), by.y = c('Name'))
  dfSummary$Unexpected <- apply(dfSummary, 1, function(r){
    htos <- unlist(strsplit(r['ExpectedHTOs'], ','))

    return(!(r['HTO'] %in% htos))
  })


  write.table(dfSummary, file = file.path(savePath, paste0('htoCompare.txt')), sep = '\t', row.names = F, quote = F)
  write.table(dfSummary[dfSummary$Unexpected & dfSummary$Fraction.HTO > 0.025,], file = file.path(savePath, paste0('UnexpectedHTOs.txt')), sep = '\t', row.names = F, quote = F)

  return(dfSummary)
}

.GenerateDataToCompareBarcodeSets <- function(workbooks, savePath = './') {
  metrics <- labkey.selectRows(
    baseUrl=Rdiscvr:::.getBaseUrl(),
    folderPath=Rdiscvr:::.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="quality_metrics",
    viewName="",
    colSort="-rowid",
    colFilter=makeFilter(
      c("category", "IN", "Cell Hashing;Cell Hashing Concordance"),
      c("workbook/workbookId", "IN", paste0(workbooks, collapse = ';'))
    ),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  summary <- data.frame(
    Name = character(),
    GEX_ReadsetId = integer(),
    HTO_Reads = integer(),
    GEX_Reads = integer(),
    CITE_Reads = integer(),
    GEX_FractionOfInputCalled = numeric(),
    GEX_InputBarcodes = numeric(),
    GEX_TotalCalledNotInInput = numeric(),
    TCR_FractionOfInputCalled = numeric(),
    TCR_InputBarcodes = numeric(),
    TCR_TotalCalledNotInInput = numeric(),
    HTO_Top_BarcodesFile = character(),
    CITE_Top_BarcodesFile = character(),
    GEX_CallsFile = character(),
    TCR_CallsFile = character(),
    ExpectedHTOs = character()
  )

  for (wb in workbooks){
    print(paste0('processing: ', wb))
    localPath <- file.path(savePath, wb)
    if (!dir.exists(localPath)){
      dir.create(localPath)
    }

    cDNAs <- labkey.selectRows(
      baseUrl=Rdiscvr:::.getBaseUrl(),
      folderPath=paste0(Rdiscvr:::.getLabKeyDefaultFolder(), wb),
      schemaName="singlecell",
      queryName="cdna_libraries",
      viewName="",
      colSort="-rowid",
      colSelect = 'rowid,readsetid,readsetid/name,hashingreadsetid,tcrreadsetid,citeseqreadsetid,hashingreadsetid/totalforwardReads,readsetid/totalforwardReads,sortId/hto,citeseqreadsetid/totalforwardReads',
      containerFilter=NULL,
      colNameOpt="rname"
    )
    print(paste0('total cDNA records: ', nrow(cDNAs)))

    callFiles <- labkey.selectRows(
      baseUrl=Rdiscvr:::.getBaseUrl(),
      folderPath=paste0(Rdiscvr:::.getLabKeyDefaultFolder(), '/', wb),
      schemaName="sequenceanalysis",
      queryName="outputfiles",
      viewName="",
      colSort="-rowid",
      colSelect="rowid,name,description,readset,readset/name,category,dataid/RowId,workbook,dataid/WebDavUrlRelative,dataid/WebDavUrlRelative,created",
      colFilter=makeFilter(c("category", "IN", "Cell Hashing Calls (VDJ);Seurat Cell Hashing Calls;Cell Hashing Counts;CITE-seq Counts")),
      containerFilter=NULL,
      colNameOpt="rname"
    )

    htoSummary <- cDNAs %>% group_by(readsetid) %>% summarise(ExpectedHTOs = paste0(sort(unique(sortid_hto)), collapse = ","))
    uniqueRs <- c()
    for (i in 1:nrow(cDNAs)) {
      row <- cDNAs[i,]
      if (row$readsetid %in% uniqueRs) {
        next
      }

      uniqueRs <- c(uniqueRs, row$readsetid)

      n <- gsub(x = row[['readsetid_name']], pattern = '-GEX', replacement = '')
      toAdd <- data.frame(Name = n, GEX_ReadsetId = row[['readsetid']], HTO_Reads = row[['hashingreadsetid_totalforwardreads']], GEX_Reads = row[['readsetid_totalforwardreads']])

      #metrics:
      htoMetrics <- metrics[metrics$readset == row[['hashingreadsetid']],]
      if (nrow(htoMetrics) > 0) {
        m <- htoMetrics[htoMetrics$metricname == 'Singlet',]
        if (nrow(m) > 0) {
          m <- m[m$dataid == max(m$dataid),]

          toAdd$HTO_Singlet <- m$metricvalue
        }
      } else {
        toAdd$HTO_Singlet <- NA
      }
      toAdd <- Rdiscvr:::.AppendMetrics(row, toAdd, metrics, 'GEX', 'readsetid')
      toAdd <- Rdiscvr:::.AppendMetrics(row, toAdd, metrics, 'TCR', 'tcrreadsetid')

      htos <- unique(htoSummary$ExpectedHTOs[htoSummary$readsetid == row[['readsetid']]])
      toAdd <- .DownloadRawCountFile(wb, callFiles, row[['hashingreadsetid']], toAdd, 'HTO_Top_BarcodesFile', savePath, category = 'Cell Hashing Counts', barcodeWhitelist = unlist(strsplit(htos, split = ',')))
      toAdd <- .DownloadRawCountFile(wb, callFiles, row[['citeseqreadsetid']], toAdd, 'CITESEQ_Top_BarcodesFile', savePath, category = 'CITE-seq Counts')

      toAdd <- Rdiscvr:::.DownloadCallFile(wb, callFiles, row[['readsetid']], toAdd, 'GEX_CallsFile', savePath, category = 'Cell Hashing Calls')
      toAdd <- Rdiscvr:::.DownloadCallFile(wb, callFiles, row[['tcrreadsetid']], toAdd, 'TCR_CallsFile', savePath, category = 'Cell Hashing Calls (VDJ)')

      #now merge HTOs:
      toAdd <- merge(toAdd, htoSummary, by.x = c('GEX_ReadsetId'), by.y = c('readsetid'), all.x = T)

      summary <- rbind(summary, toAdd)
    }
  }

  return(summary)
}

.AppendMetrics <- function(row, toAdd, metrics, type, field) {
  for (name in c('FractionOfInputCalled', 'InputBarcodes', 'TotalCalledNotInInput')) {
    toAdd[[paste0(type, '_', name)]] <- NA
  }

  xMetrics <- metrics[metrics$readset == row[[field]],]
  if (nrow(xMetrics) > 0) {
    latest <- xMetrics[xMetrics$dataid == max(xMetrics$dataid),]

    for (name in c('FractionOfInputCalled', 'InputBarcodes', 'TotalCalledNotInInput')) {
      toAdd[[paste0(type, '_', name)]] <- latest[latest$metricname == name,]$metricvalue
    }
  }

  return(toAdd)
}

.ProcessSet <- function(df, set1, set2, type1, type2) {
  for (name1 in names(set1)) {
    h <- set1[[name1]]

    for (name2 in names(set2)) {
      g <- set2[[name2]]
      i <- length(intersect(h, g))

      df <- rbind(df, data.frame(dataset1 = name1, type1 = type1, dataset2 = name2, type2 = type2, intersect = i, length1 = length(h), length2 = length(g)))
    }
  }

  return(df)
}

.CompareHtosByCall <- function(name, htoCallFile, compareFile, type) {
  if (is.na(compareFile)){
    print(paste0('Missing ', type, ' call file for: ', name))
    return(NA)
  }

  t1 <- read.table(htoCallFile, header = T, sep = '\t')
  t1 <- t1 %>% group_by(HTO) %>% summarise(Count = dplyr::n())
  t1 <- t1[!(t1$HTO %in% c('Doublet', 'Negative', 'Discordant')),]
  t1$Fraction <- t1$Count / sum(t1$Count)

  t2 <- read.table(compareFile, header = T, sep = '\t')
  t2 <- t2 %>% group_by(HTO) %>% summarise(Count = dplyr::n())
  t2 <- t2[!(t2$HTO %in% c('Doublet', 'Negative', 'Discordant')),]
  t2$Fraction <- t2$Count / sum(t2$Count)

  df <- merge(x = t1, y = t2, by = 'HTO', all = T, suffixes = c('.HTO', '.Compare'))
  df$Count.HTO[is.na(df$Count.HTO)] <- 0
  df$Fraction.HTO[is.na(df$Fraction.HTO)] <- 0

  df$Count.Compare[is.na(df$Count.Compare)] <- 0
  df$Fraction.Compare[is.na(df$Fraction.Compare)] <- 0

  df$Type <- type
  df$Difference <- abs(df$Fraction.HTO - df$Fraction.Compare)
  df$Dataset <- name

  df <- df[c('Dataset', 'HTO', 'Type', 'Count.HTO', 'Fraction.HTO', 'Count.Compare', 'Fraction.Compare', 'Difference')]

  return(df)
}

.DownloadRawCountFile <- function(wb, callFiles, readsetId, toAdd, fieldName, localPath, category, barcodeWhitelist = NULL) {
  toAdd[fieldName] <- NA
  
  cf <- callFiles[callFiles$readset == readsetId & callFiles$category == category,]
  if (nrow(cf) > 1) {
    print(paste0('Multiple call files, using latest: ', toAdd$Name, ' ', fieldName))
  }
  
  #use the most recent (highest rowId)
  if (nrow(cf) > 0) {
    row <- cf[cf$rowid == max(cf$rowid),]
    
    suffix <- '.rawCountDir'

    fn <- paste0(row['readset_name'], '.', row['readset'], '.', row['rowid'], '.', row$category, suffix)
    fn <- gsub('\\(', '_', fn)
    fn <- gsub('\\)', '_', fn)
    fn <- gsub(' ', '_', fn)
    fn <- gsub('__', '_', fn)
    
    lp <- paste0(localPath, '/', wb, '/', fn)
    if (!dir.exists(lp)) {
      dir.create(lp)
    }
    
    expectedDir <- paste0(lp, '/raw_feature_bc_matrix')
    outFile <- paste0(lp, '/countsPerCell.txt')
    if (file.exists(outFile)) {
      print(paste0('file exists, reusing: ', outFile))
    } else {
      Rdiscvr::DownloadOutputDirectoryFromOutputFile(outputFileId = row[['rowid']], outFile = lp, overwrite = T, pathTranslator = function(x){
        return(dirname(x))
      })
    }
    
    # Read dir, find top barcodes, save to file:
    mat <- as.matrix(cellhashR::ProcessCountMatrix(rawCountData = expectedDir, barcodeWhitelist = barcodeWhitelist))

    sortedMat <- colSums(mat)
    names(sortedMat) <- sapply(colnames(mat), function(x){
      x <- unlist(strsplit(x, split = '-'))[1]
      return(x)
    })
    sortedMat <- sort(sortedMat, decreasing = T)
    
    write.table(data.frame(cellbarcode = names(sortedMat), count = unname(sortedMat)), file = outFile, sep = '\t', row.names = F)
    
    toAdd[fieldName] <- outFile
  } else {
    print('No count dir found')
    toAdd[fieldName] <- NA
  }
  
  return(toAdd)
}

.DownloadCallFile <- function(wb, callFiles, readsetId, toAdd, fieldName, localPath, category) {
  toAdd[fieldName] <- NA

  cf <- callFiles[callFiles$readset == readsetId & callFiles$category == category,]
  if (nrow(cf) > 1) {
    print(paste0('Multiple call files, using latest: ', toAdd$Name, ' ', fieldName))
  }

  #use the most recent (highest rowId)
  if (nrow(cf) > 0) {
    row <- cf[cf$rowid == max(cf$rowid),]

    suffix <- '.calls'

    fn <- paste0(row['readset_name'], '.', row['readset'], '.', row['rowid'], '.', row$category, suffix, '.txt')
    fn <- gsub('\\(', '_', fn)
    fn <- gsub('\\)', '_', fn)
    fn <- gsub(' ', '_', fn)
    fn <- gsub('__', '_', fn)

    remotePath <- row[['dataid_webdavurlrelative']]
    lp <- paste0(localPath, '/', wb, '/', fn)
    if (file.exists(lp)) {
      print(paste0('file exists, reusing: ', lp))
    } else {
      success <- labkey.webdav.get(
        baseUrl=Rdiscvr:::.getBaseUrl(),
        folderPath=paste0(Rdiscvr:::.getLabKeyDefaultFolder(), wb),
        remoteFilePath = remotePath,
        overwrite = T,
        localFilePath = lp
      )

      if (!success) {
        warning(paste0('Unable to download whitelist for readset: ', row['readset_name'], ' from: ', remotePath))
      }
    }

    toAdd[fieldName] <- lp
  }

  return(toAdd)
}
