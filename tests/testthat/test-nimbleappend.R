context("Nimble")

test_that("Nimble Append detects missing input file", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  expect_error(AppendNimbleCounts(seuratObj, "../testdata/nonexistent.tsv", targetAssayName = 'RNA'), "Nimble file does not exist: ../testdata/nonexistent.tsv", fixed=TRUE)
})

test_that("Nimble Append deletes blank feature names when appending", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'Nimble', maxLibrarySizeRatio = 1)
  expect_false('' %in% rownames(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')))
  expect_true('FeatureSource' %in% names(seuratObj@assays$Nimble@meta.features))
  expect_equal(unique(seuratObj@assays$Nimble@meta.features$FeatureSource), c('Nimble'))
})

test_that("Nimble Append works with empty input", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  fn <- tempfile(fileext = '.txt')
  file.create(fn)
  seuratObj <- AppendNimbleCounts(seuratObj, fn, targetAssayName = 'Nimble')
  unlink(fn)
  expect_false('Nimble' %in% names(seuratObj@assays))
})

test_that("Nimble Append deletes all barcodes not in Seurat when appending", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  nimbleExclusiveBarcodes <- c("12345_CCAGCGAAGTCCGTAT", "12345_CCAGCGAAGTCCGTAC")
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'RNA', maxLibrarySizeRatio = 1, replaceExistingAssayData = FALSE)
  expect_equal(nimbleExclusiveBarcodes %in% colnames(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')), c(FALSE, FALSE))
  
  expect_true('FeatureSource' %in% names(seuratObj@assays$RNA@meta.features))
  expect_equal(seuratObj@assays$RNA@meta.features$FeatureSource, c("CellRanger", "CellRanger", "CellRanger", "CellRanger", "Nimble", "Nimble", "Nimble", "Nimble"))
})

test_that("Nimble Append fills all barcodes in Seurat but not in Nimble when appending", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  seuratExclusiveBarcode <- "12345_TAAAAAAAAAAAAAAA"
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'RNA', maxLibrarySizeRatio = 1, replaceExistingAssayData = FALSE)
  nimbleFeatureCounts <- GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, seuratExclusiveBarcode][c('E', 'F', 'G', 'H')]
  expect_false(FALSE %in% (nimbleFeatureCounts == c(0, 0, 0, 0)))
})

test_that("Nimble Append output Seurat object is valid when appending", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  expectedBarcodes <- c("12345_AAAAAAAAAAAAAAAA", "12345_CAAAAAAAAAAAAAA", "12345_GAAAAAAAAAAAAAAA", "12345_TAAAAAAAAAAAAAAA")
  expectedFeatures <- c("A", "B", "C", "D", "E", "F", "G", "H")
  expectedValues <- list(c(1, 1, 1, 1, 1, 1, 1, 1), c(1, 1, 1, 1, 2, 2, 2, 2), c(1, 1, 1, 1, 4, 4, 4, 4), c(1, 1, 1, 1, 0, 0, 0, 0))
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'RNA', maxLibrarySizeRatio = 1, replaceExistingAssayData = FALSE)
  expect_equal(colnames(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')), expectedBarcodes)
  expect_equal(rownames(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')), expectedFeatures)
  
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[1]]), expectedValues[[1]])
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[2]]), expectedValues[[2]])
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[3]]), expectedValues[[3]])
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[4]]), expectedValues[[4]])
})

test_that("Nimble Append deletes blank feature names when creating new assay", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'Nimble', maxLibrarySizeRatio = 1)
  expect_false('' %in% rownames(GetAssayData(seuratObj, assay = 'Nimble', slot = 'counts')))
})

test_that("Nimble Append deletes all barcodes not in Seurat when creating new assay", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  nimbleExclusiveBarcodes <- c("12345_CCAGCGAAGTCCGTAT", "12345_CCAGCGAAGTCCGTAC")
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'Nimble', maxLibrarySizeRatio = 1)
  expect_equal(nimbleExclusiveBarcodes %in% colnames(GetAssayData(seuratObj, assay = 'Nimble', slot = 'counts')), c(FALSE, FALSE))
})

test_that("Nimble Append fills all barcodes in Seurat but not in Nimble when creating new assay", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  seuratExclusiveBarcode <- "12345_TAAAAAAAAAAAAAAA"
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'Nimble', maxLibrarySizeRatio = 1)
  nimbleFeatureCounts <- GetAssayData(seuratObj, assay = 'Nimble', slot = 'counts')[, seuratExclusiveBarcode]
  expect_equal(FALSE %in% (nimbleFeatureCounts == c(0, 0, 0, 0)), FALSE)
})

test_that("Nimble Append output Seurat object is valid when creating new assay", {
  seuratObj <- Seurat::UpdateSeuratObject(readRDS("../testdata/nimbleTest.rds"))
  expectedBarcodes <- c("12345_AAAAAAAAAAAAAAAA", "12345_CAAAAAAAAAAAAAA", "12345_GAAAAAAAAAAAAAAA", "12345_TAAAAAAAAAAAAAAA")
  expectedFeatures <- c("E", "F", "G", "H")
  expectedValues <- list(c(1, 1, 1, 1), c(2, 2, 2, 2), c(4, 4, 4, 4), c(0, 0, 0, 0))
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", targetAssayName = 'Nimble', maxLibrarySizeRatio = 1)
  expect_equal(colnames(GetAssayData(seuratObj, assay = 'Nimble', slot = 'counts')), expectedBarcodes)
  expect_equal(rownames(GetAssayData(seuratObj, assay = 'Nimble', slot = 'counts')), expectedFeatures)
  
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[1]]), expectedValues[[1]])
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[2]]), expectedValues[[2]])
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[3]]), expectedValues[[3]])
  expect_equal(unname(GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[, expectedBarcodes[4]]), expectedValues[[4]])
})
