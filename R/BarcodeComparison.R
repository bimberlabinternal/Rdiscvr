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
#' @param filePrefix A prefix appended to the name of all summary output files
#' @export
#' @importFrom dplyr %>% mutate group_by
CompareCellBarcodeSets <- function(workbooks, savePath = tempdir(), filePrefix = '') {
  summary <- .GenerateDataToCompareBarcodeSets(workbooks, savePath)
  write.table(summary, file = file.path(savePath, paste0(filePrefix, 'rawComparisonData.txt')), sep = '\t', row.names = F, quote = F)

  htoBC <- list()
  htoBC_All <- list()
  gexBC <- list()
  tcrBC <- list()
  citeBC <- list()

  for (i in 1:nrow(summary)) {
    row <- summary[i,]
    name <- as.character(row$Name)
    if (is.na(name)) {
      print('skipping NA row')
      next
    }

    if (!is.na(row$HTO_Top_BarcodesFile) && !is.na(row$GEX_CallsFile)){
      bc <- read.table(row$HTO_Top_BarcodesFile, header = T, sep = '\t')
      htoBC[[name]] <- bc$cellbarcode

      bc <- read.table(row$HTO_Top_BarcodesFile_All, header = T, sep = '\t')
      htoBC_All[[name]] <- bc$cellbarcode

      bc <- read.table(row$GEX_CallsFile, header = T, sep = '\t')$cellbarcode
      gexBC[[name]] <- bc
    } else {
      print(paste0('Missing one or more required files, skipping: ', name, ':'))
      if (is.na(row$HTO_Top_BarcodesFile)){
        print(paste0('HTO missing: ', name))
      }

      if (is.na(row$GEX_CallsFile)){
        print(paste0('GEX missing: ', name))
      }

      next
    }

    if (!is.na(row$TCR_CallsFile)) {
      print(paste0('processing TCR:', name))
      bc <- read.table(row$TCR_CallsFile, header = T, sep = '\t')$cellbarcode
      tcrBC[[name]] <- bc
    }

    if (!is.na(row$CITE_Top_BarcodesFile)) {
      print(paste0('processing CITE-seq:', name))
      bc <- read.table(row$CITE_Top_BarcodesFile, header = T, sep = '\t')
      citeBC[[name]] <- bc$cellbarcode
    }
  }

  print('starting summary')
  print(paste0('Total GEX: ', length(gexBC)))
  print(paste0('Total VDJ: ', length(tcrBC)))
  print(paste0('Total HTO: ', length(htoBC)))
  print(paste0('Total CITE: ', length(citeBC)))

  df <- data.frame(dataset1 = character(), type1 = character(), dataset2 = character(), type2 = character(), intersect = integer(), length1 = integer(), length2 = integer())

  df <- .ProcessSet(df, gexBC, htoBC, 'GEX', 'HTO')
  df <- .ProcessSet(df, gexBC, htoBC, 'GEX', 'HTO-All')
  if (length(tcrBC) > 0) {
    df <- .ProcessSet(df, gexBC, tcrBC, 'GEX', 'TCR')
  }
  if (length(citeBC) > 0) {
    df <- .ProcessSet(df, gexBC, citeBC, 'GEX', 'CITE')
  }

  if (length(tcrBC) > 0) {
    df <- .ProcessSet(df, tcrBC, htoBC, 'TCR', 'HTO')
    df <- .ProcessSet(df, tcrBC, htoBC, 'TCR', 'HTO-All')
    df <- .ProcessSet(df, tcrBC, gexBC, 'TCR', 'GEX')
  }

  if (length(citeBC) > 0 && length(tcrBC) > 0) {
    df <- .ProcessSet(df, tcrBC, citeBC, 'TCR', 'CITE')
  }

  df$fraction <- df$intersect / df$length1
  df <- df[order(df$dataset1, df$type1, df$type2, -df$intersect),]
  df <- df %>% group_by(dataset1, type1, type2) %>% mutate(max_intersect_by_type = max(intersect))
  df <- df %>% group_by(dataset1, type1) %>% mutate(max_intersect = max(intersect))

  write.table(df, file = file.path(savePath, paste0(filePrefix, 'cell_barcode_comparisons.txt')), sep = '\t', row.names = F, quote = F)

  self <- df[df$dataset1 == df$dataset2, c('dataset1', 'type1', 'type2', 'intersect', 'fraction')]
  names(self) <- c('dataset1', 'type1', 'type2', 'self_intersect', 'self_intersect_fraction')

  df2 <- df[df$intersect == df$max_intersect_by_type,]
  write.table(df2, file = file.path(savePath, paste0(filePrefix, 'top_intersect_by_type.txt')), sep = '\t', row.names = F, quote = F)

  df3 <- df2[df2$dataset1 != df2$dataset2,]
  df3 <- merge(df3, self, by = c('dataset1', 'type1', 'type2'), all.x = T, all.y = F)
  #filter ties:
  df3 <- df3[df3$intersect != df3$self_intersect,]
  write.table(df3, file = file.path(savePath, paste0(filePrefix, 'conflicting_intersect.txt')), sep = '\t', row.names = F, quote = F)
}

.GenerateDataToCompareBarcodeSets <- function(workbooks, savePath = tempdir()) {
  metrics <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=.getLabKeyDefaultFolder(),
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
    HTO_Top_BarcodesFile_All = character(),
    CITE_Top_BarcodesFile = character(),
    GEX_CallsFile = character(),
    TCR_CallsFile = character(),
    ExpectedHTOs = character()
  )

  allWorkbooks <- labkey.selectRows(
    baseUrl=.getBaseUrl(),
    folderPath=paste0(.getLabKeyDefaultFolder()),
    schemaName="core",
    queryName="containers",
    colFilter = makeFilter(c('type', 'EQUALS', 'workbook')),
    colSelect = 'name',
    containerFilter=NULL,
    colNameOpt="rname"
  )$name

  for (wb in workbooks){
    print(paste0('processing: ', wb))
    if (!wb %in% allWorkbooks) {
      print('Doesnt exist, skipping')
      next
    }

    localPath <- file.path(savePath, wb)
    if (!dir.exists(localPath)){
      dir.create(localPath)
    }

    cDNAs <- suppressWarnings(labkey.selectRows(
      baseUrl=.getBaseUrl(),
      folderPath=paste0(.getLabKeyDefaultFolder(), wb),
      schemaName="singlecell",
      queryName="cdna_libraries",
      viewName="",
      colSort="-rowid",
      colSelect = 'plateid,readsetid,readsetid/name,hashingreadsetid,tcrreadsetid,citeseqreadsetid,hashingreadsetid/totalforwardReads,readsetid/totalforwardReads,sortId/hto,citeseqreadsetid/totalforwardReads',
      colFilter = makeFilter(c('readsetid', 'NOT_MISSING', '')),
      colNameOpt="rname"
    ))
    print(paste0('total cDNA records: ', nrow(cDNAs)))

    callFiles <- suppressWarnings(labkey.selectRows(
      baseUrl=.getBaseUrl(),
      folderPath=paste0(.getLabKeyDefaultFolder(), '/', wb),
      schemaName="sequenceanalysis",
      queryName="outputfiles",
      viewName="",
      colSort="-rowid",
      colSelect="rowid,name,description,readset,readset/name,category,dataid/RowId,workbook,dataid/WebDavUrlRelative,dataid/WebDavUrlRelative,created",
      colFilter=makeFilter(
        c("category", "IN", "Cell Hashing Calls (VDJ);Seurat Cell Hashing Calls;Cell Hashing Counts;CITE-seq Counts;Seurat Object Prototype;10x Loupe File")
      ),
      containerFilter=NULL,
      colNameOpt="rname"
    ))

    htoSummary <- cDNAs %>% group_by(readsetid) %>% summarise(ExpectedHTOs = paste0(sort(unique(sortid_hto)), collapse = ","))
    uniqueRs <- c()

    cDNAs <- cDNAs[names(cDNAs) != 'sortid_hto']
    cDNAs <- unique(cDNAs)
    print(paste0('unique cDNAs after collapse: ', nrow(cDNAs)))
    if (nrow(cDNAs) == 0) {
      print('no cDNAs')
      next
    }

    for (i in 1:nrow(cDNAs)) {
      row <- cDNAs[i,]
      print(paste0('processing: ', row$plateid))
      if (row$readsetid %in% uniqueRs) {
        next
      }

      uniqueRs <- c(uniqueRs, row$readsetid)

      n <- gsub(x = row[['readsetid_name']], pattern = '-GEX', replacement = '')
      toAdd <- data.frame(Name = n, GEX_ReadsetId = row[['readsetid']], HTO_Reads = row[['hashingreadsetid_totalforwardreads']], GEX_Reads = row[['readsetid_totalforwardreads']])

      #metrics:
      if (!is.null(row[['hashingreadsetid']]) && !is.na(row[['hashingreadsetid']])) {
        htoMetrics <- metrics[metrics$readset == row[['hashingreadsetid']],]
        if (nrow(htoMetrics) > 0) {
          print('adding hashing metrics')
          m <- htoMetrics[htoMetrics$metricname == 'Singlet',]
          if (nrow(m) > 0) {
            m <- m[m$dataid == max(m$dataid),]

            toAdd$HTO_Singlet <- m$metricvalue
          }
        } else {
          toAdd$HTO_Singlet <- NA
        }
      } else {
        toAdd$HTO_Singlet <- NA
      }

      toAdd <- .AppendMetrics(row, toAdd, metrics, 'GEX', 'readsetid')
      toAdd <- .AppendMetrics(row, toAdd, metrics, 'TCR', 'tcrreadsetid')

      htos <- unique(htoSummary$ExpectedHTOs[htoSummary$readsetid == row[['readsetid']]])
      toAdd <- .DownloadRawCountFile(wb, callFiles, row[['hashingreadsetid']], toAdd, 'HTO_Top_BarcodesFile', savePath, category = 'Cell Hashing Counts', barcodeWhitelist = unlist(strsplit(htos, split = ',')))
      toAdd <- .DownloadRawCountFile(wb, callFiles, row[['hashingreadsetid']], toAdd, 'HTO_Top_BarcodesFile_All', savePath, category = 'Cell Hashing Counts')

      toAdd <- .DownloadRawCountFile(wb, callFiles, row[['citeseqreadsetid']], toAdd, 'CITE_Top_BarcodesFile', savePath, category = 'CITE-seq Counts')

      toAdd <- .DownloadCallFile(wb, callFiles, row[['readsetid']], toAdd, 'GEX_CallsFile', savePath, category = 'Seurat Cell Hashing Calls')
      if (is.na(toAdd$GEX_CallsFile)) {
        print('Downloading GEX barcodes from seurat prototype')
        .DownloadBarcodesForSeurat(wb, callFiles, row[['readsetid']], toAdd, 'GEX_CallsFile', savePath, category = 'Seurat Object Prototype')
      }

      if (is.na(toAdd$GEX_CallsFile)) {
        print('Downloading GEX barcodes from loupe file')
        toAdd <- .DownloadBarcodesForLoupe(wb, callFiles, row[['readsetid']], toAdd, 'GEX_CallsFile', savePath, category = '10x Loupe File')
      }

      if (is.na(toAdd$GEX_CallsFile)) {
        print('Missing GEX Calls')
      }

      toAdd <- .DownloadCallFile(wb, callFiles, row[['tcrreadsetid']], toAdd, 'TCR_CallsFile', savePath, category = 'Cell Hashing Calls (VDJ)')
      if (is.na(row[['tcrreadsetid']]) && is.na(toAdd$TCR_CallsFile)) {
        print('Missing TCR Calls')
      }

      #now merge HTOs:
      if (!is.null(row[['hashingreadsetid']]) && !is.na(row[['hashingreadsetid']])) {
        toAdd <- merge(toAdd, htoSummary, by.x = c('GEX_ReadsetId'), by.y = c('readsetid'), all.x = T)
      } else {
        toAdd$ExpectedHTOs <- NA
      }

      summary <- rbind(summary, toAdd)
    }
  }

  return(summary)
}

.AppendMetrics <- function(row, toAdd, metrics, type, field) {
  for (name in c('FractionOfInputCalled', 'InputBarcodes', 'TotalCalledNotInInput')) {
    toAdd[[paste0(type, '_', name)]] <- NA
  }

  # readset type not present
  if (is.na(row[[field]]) || is.null(row[[field]])) {
    return(toAdd)
  }

  xMetrics <- metrics[metrics$readset == row[[field]],]
  if (nrow(xMetrics) > 0) {
    latest <- xMetrics[xMetrics$dataid == max(xMetrics$dataid),]
    if (nrow(latest) == 0) {
      return(toAdd)
    }

    for (name in c('FractionOfInputCalled', 'InputBarcodes', 'TotalCalledNotInInput')) {
      toAdd[[paste0(type, '_', name)]] <- latest[latest$metricname == name,]$metricvalue
    }
  }

  return(toAdd)
}

.ProcessSet <- function(df, set1, set2, type1, type2) {
  print(paste0('Processing: ', type1, ' / ', type2))

  for (name1 in names(set1)) {
    h <- set1[[name1]]
    if (all(is.na(h)) || all(is.null(h))) {
      print(paste0('skipping: ', name1))
      next
    }

    for (name2 in names(set2)) {
      g <- set2[[name2]]
      if (all(is.na(g)) || all(is.null(g))) {
        print(paste0('skipping: ', name2))
        next
      }

      if ('CITE' == type1 || 'HTO' == type1 || 'HTO-All' == type1) {
        h <- h[1:length(g)]
      } else if ('CITE' == type2 || 'HTO' == type2 || 'HTO-All' == type2) {
        g <- g[1:length(h)]
      }

      i <- length(intersect(h, g))

      toAdd <- data.frame(dataset1 = name1, type1 = type1, dataset2 = name2, type2 = type2, intersect = i, length1 = length(h), length2 = length(g))
      if (any(is.na(toAdd$dataset1))) {
        print('NA Value')
        print(name1)
        print(name2)
        print(toAdd)
      }

      df <- rbind(df, toAdd)
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
  t1$HTO <- t1$consensuscall

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
  if (is.null(readsetId) || is.na(readsetId)) {
    return(toAdd)
  }
  
  cf <- callFiles[!is.na(callFiles$readset) & callFiles$readset == readsetId & callFiles$category == category,]
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

    if (all(is.null(barcodeWhitelist))) {
      outFile <- paste0(lp, '/countsPerCell.all.txt')
    } else {
      outFile <- paste0(lp, '/countsPerCell.subset.txt')
    }

    if (file.exists(outFile)) {
      print(paste0('file exists, reusing: ', outFile))
    } else {
      Rdiscvr::DownloadOutputDirectoryFromOutputFile(outputFileId = row[['rowid']], outFile = lp, overwrite = T, pathTranslator = function(x){
        return(dirname(x))
      })

      # Read dir, find top barcodes, save to file:
      print('reading matrix')
      origMat <- Seurat::Read10X(expectedDir, gene.column=1, strip.suffix = TRUE)
      origMat <- as.matrix(origMat)

      if (!is.null(barcodeWhitelist)) {
        mat <- origMat[barcodeWhitelist,,drop=FALSE]
      } else {
        mat <- origMat
      }

      sortedMat <- colSums(mat)
      names(sortedMat) <- colnames(mat)
      sortedMat <- sort(sortedMat, decreasing = T)
    
      write.table(data.frame(cellbarcode = names(sortedMat), count = unname(sortedMat)), file = outFile, sep = '\t', row.names = F, quote = F)
    }

    toAdd[fieldName] <- outFile
  } else {
    print('No count dir found')
    toAdd[fieldName] <- NA
  }
  
  return(toAdd)
}

.DownloadCallFile <- function(wb, callFiles, readsetId, toAdd, fieldName, localPath, category) {
  toAdd[fieldName] <- NA

  cf <- callFiles[!is.na(callFiles$readset) & callFiles$readset == readsetId & callFiles$category == category,]
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
        baseUrl=.getBaseUrl(),
        folderPath=paste0(.getLabKeyDefaultFolder(), wb),
        remoteFilePath = remotePath,
        overwrite = T,
        localFilePath = lp
      )

      if (!success) {
        warning(paste0('Unable to download whitelist for readset: ', row['readset_name'], ' from: ', remotePath))
      }
    }

    toAdd[fieldName] <- lp
  } else {
    print(paste0('No file of type: ', category))
  }

  return(toAdd)
}

.DownloadBarcodesForSeurat <- function(wb, callFiles, readsetId, toAdd, fieldName, localPath, category) {
  toAdd[fieldName] <- NA

  cf <- callFiles[!is.na(callFiles$readset) & callFiles$readset == readsetId & callFiles$category == category,]
  if (nrow(cf) > 1) {
    print(paste0('Multiple call files, using latest: ', toAdd$Name, ' ', fieldName))
  }

  if (nrow(cf) > 0) {
    row <- cf[cf$rowid == max(cf$rowid),]

    suffix <- '.seuratBarcodes'
    fn <- paste0(row['readset_name'], '.', row['readset'], '.', row['rowid'], '.', row$category, suffix, '.txt')
    fn <- gsub('\\(', '_', fn)
    fn <- gsub('\\)', '_', fn)
    fn <- gsub(' ', '_', fn)
    fn <- gsub('__', '_', fn)
    lp <- paste0(localPath, '/', wb, '/', fn)
    if (file.exists(lp)) {
      print(paste0('file exists, reusing: ', lp))
    } else {

      barcodeOut <- paste0(row['readset_name'], '.', row['readset'], '.', row['rowid'], '.', row$category, '.cellBarcodes.csv')
      DownloadOutputFile(row[['rowid']], overwrite = T, outFile = barcodeOut, pathTranslator = function(x){
        x <- gsub(x, pattern = 'seurat.rds', replacement = 'cellBarcodes.csv')
        return(x)
      })

      data <- read.table(barcodeOut, sep = ',', header = F)
      names(data) <- c('cellbarcode')
      data$cellbarcode <- sapply(data$cellbarcode, function(x){
        x <- unlist(strsplit(x, split = '_'))[2]
        return(x)
      })

      write.table(data, file = lp, sep = '\t', row.names = F, quote = F)
    }

    toAdd[fieldName] <- lp
  } else {
    print(paste0('No file of type: ', category))
  }

  return(toAdd)
}

.DownloadBarcodesForLoupe <- function(wb, callFiles, readsetId, toAdd, fieldName, localPath, category) {
  toAdd[fieldName] <- NA

  cf <- callFiles[!is.na(callFiles$readset) & callFiles$readset == readsetId & callFiles$category == category,]
  if (nrow(cf) > 1) {
    print(paste0('Multiple call files, using latest: ', toAdd$Name, ' ', fieldName))
  }

  #use the most recent (highest rowId)
  if (nrow(cf) > 0) {
    row <- cf[cf$rowid == max(cf$rowid),]

    suffix <- '.loupe'

    fn <- paste0(row['readset_name'], '.', row['readset'], '.', row['rowid'], '.', row$category, suffix, '.txt')
    fn <- gsub('\\(', '_', fn)
    fn <- gsub('\\)', '_', fn)
    fn <- gsub(' ', '_', fn)
    fn <- gsub('__', '_', fn)

    lp <- paste0(localPath, '/', wb, '/', fn)
    if (!dir.exists(lp)) {
      dir.create(lp)
    }

    expectedDir <- paste0(lp, '/filtered_feature_bc_matrix')
    outFile <- paste0(lp, '/countsPerCell.loupe.txt')
    if (file.exists(outFile)) {
      print(paste0('file exists, reusing: ', outFile))
    } else {
      Download10xRawDataForLoupeFile(outputFileId = row[['rowid']], overwrite = T, outFile = lp, countType = 'filtered_feature_bc_matrix')

      # Read dir, find top barcodes, save to file:
      print('reading matrix from loupe')
      mat <- Seurat::Read10X(expectedDir, gene.column=1, strip.suffix = TRUE)
      mat <- as.matrix(mat)

      sortedMat <- colSums(mat)
      names(sortedMat) <- colnames(mat)
      sortedMat <- sort(sortedMat, decreasing = T)

      write.table(data.frame(cellbarcode = names(sortedMat), count = unname(sortedMat)), file = outFile, sep = '\t', row.names = F, quote = F)
    }

    toAdd[fieldName] <- outFile
  } else {
    print(paste0('No file of type: ', category))
  }

  return(toAdd)
}