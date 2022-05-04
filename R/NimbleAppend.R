#' @import Seurat tidyr Matrix
<<<<<<< HEAD
#' 
#' @title AppendNimbleCounts
#' @description Reads a given seurat object and a nimble file, and appends the nimble data to the object
#'
#' @param seuratObject A Seurat object.
#' @param file A nimble file
#' @param appendToCounts If true, append the nimble data to the count matrix. If false, make a new Nimble assay
#' @return A modified Seurat object.
#' @export
=======

>>>>>>> f80769446904478a9703977ec0f82736fdc66508
AppendNimbleCounts <- function(seuratObject, file, appendToCounts=FALSE) {
  if (!file.exists(file)) {
    stop("Nimble file not found.")
  }
  
  # Read file and construct df
  df <- read.table(file, sep="\t", header=FALSE)
  
<<<<<<< HEAD
  if(length(df[!(df$V1==""), ]) != 0) {
    # Remove blank feature names, just in case
    df <- df[!(df$V1==""), ]
    
    warning("Nimble data contains blank feature names")
  }
  
  #Remove ambiguous features
  df <- df[!grepl(",", df$V1), ]

=======
  # Remove blank feature names, just in case
  df <- df[!(df$V1==""), ]
  
  #Remove ambiguous features
  df <- df[!grepl(",", df$V1), ]
>>>>>>> f80769446904478a9703977ec0f82736fdc66508
  
  df <- tidyr::pivot_wider(df, names_from=V3, values_from=V2, values_fill=0)
  
  # Trim seurat object and get barcodes
  if(appendToCounts == TRUE) {
    seuratObject <- DietSeurat(seuratObject)
  }
  
  seuratBarcodes <- seuratObject@assays$RNA@counts@Dimnames[2][[1]]
  
  # Remove barcodes from nimble that aren't in seurat
  barcodeDiff <- colnames(df) %in% seuratBarcodes
  barcodeDiff[1] <- TRUE
  df <- df[barcodeDiff]
  
  # Fill zeroed barcodes that are in seurat but not in nimble
  zeroedBarcodes <- setdiff(seuratBarcodes, colnames(df)[-1])
  
  for(barcode in zeroedBarcodes) {
    df[barcode] <- 0
  }
  
  # Cast nimble df to matrix
  featureNames <- df$V1
  barcodes <- colnames(df)[-1]
  df <- subset(df, select=-(V1))
  m <- Reduce(cbind2, lapply(df, Matrix, sparse = TRUE))
  dimnames(m) <- list(featureNames, barcodes)
  
  
  if(appendToCounts == TRUE) {
    # Append nimble matrix to seurat count matrix
    seuratObject@assays$RNA@counts <- rbind(seuratObject@assays$RNA@counts, m)
  } else {
    # Add nimble as separate assay
    seuratObject[["Nimble"]] <- CreateAssayObject(counts = m)
  }
  
  return(seuratObject)
}
