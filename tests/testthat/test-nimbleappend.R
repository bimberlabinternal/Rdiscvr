context("Nimble")


test_that("Nimble Append detects missing input file", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  expect_error(AppendNimbleCounts(seuratObj, "../testdata/nonexistent.tsv"), "Nimble file not found.", fixed=TRUE)
})

test_that("Nimble Append deletes blank feature names when appending", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=FALSE)
  expect_equal('' %in% rownames(seuratObj@assays$RNA@counts), FALSE)
})

test_that("Nimble Append deletes all barcodes not in Seurat when appending", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  nimbleExclusiveBarcodes = c("12345_CCAGCGAAGTCCGTAT", "12345_CCAGCGAAGTCCGTAC")
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=TRUE)
  expect_equal(nimbleExclusiveBarcodes %in% colnames(seuratObj@assays$RNA@counts), c(FALSE, FALSE))
})

test_that("Nimble Append fills all barcodes in Seurat but not in Nimble when appending", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  seuratExclusiveBarcode = "12345_TAAAAAAAAAAAAAAA"
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=TRUE)
  nimbleFeatureCounts <- seuratObj@assays$RNA@counts[, seuratExclusiveBarcode][c('E', 'F', 'G', 'H')]
  expect_equal(FALSE %in% (nimbleFeatureCounts == c(0, 0, 0, 0)), FALSE)
})

test_that("Nimble Append output Seurat object is valid when appending", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  expectedBarcodes <- c("12345_AAAAAAAAAAAAAAAA", "12345_CAAAAAAAAAAAAAA", "12345_GAAAAAAAAAAAAAAA", "12345_TAAAAAAAAAAAAAAA")
  expectedFeatures <- c("A", "B", "C", "D", "E", "F", "G", "H")
  expectedValues <- list(c(1, 1, 1, 1, 1, 1, 1, 1), c(1, 1, 1, 1, 2, 2, 2, 2), c(1, 1, 1, 1, 4, 4, 4, 4), c(1, 1, 1, 1, 0, 0, 0, 0))
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=TRUE)
  expect_equal(colnames(seuratObj@assays$RNA@counts), expectedBarcodes)
  expect_equal(rownames(seuratObj@assays$RNA@counts), expectedFeatures)
  expect_equal(FALSE %in% (seuratObj@assays$RNA@counts[, expectedBarcodes[1]] == expectedValues[[1]]), FALSE)
  expect_equal(FALSE %in% (seuratObj@assays$RNA@counts[, expectedBarcodes[2]] == expectedValues[[2]]), FALSE)
  expect_equal(FALSE %in% (seuratObj@assays$RNA@counts[, expectedBarcodes[3]] == expectedValues[[3]]), FALSE)
  expect_equal(FALSE %in% (seuratObj@assays$RNA@counts[, expectedBarcodes[4]] == expectedValues[[4]]), FALSE)
})

test_that("Nimble Append deletes blank feature names when creating new assay", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=FALSE)
  expect_equal('' %in% rownames(seuratObj@assays$Nimble@counts), FALSE)
})

test_that("Nimble Append deletes all barcodes not in Seurat when creating new assay", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  nimbleExclusiveBarcodes = c("12345_CCAGCGAAGTCCGTAT", "12345_CCAGCGAAGTCCGTAC")
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=FALSE)
  expect_equal(nimbleExclusiveBarcodes %in% colnames(seuratObj@assays$Nimble@counts), c(FALSE, FALSE))
})

test_that("Nimble Append fills all barcodes in Seurat but not in Nimble when creating new assay", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  seuratExclusiveBarcode = "12345_TAAAAAAAAAAAAAAA"
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=FALSE)
  nimbleFeatureCounts <- seuratObj@assays$Nimble@counts[, seuratExclusiveBarcode]
  expect_equal(FALSE %in% (nimbleFeatureCounts == c(0, 0, 0, 0)), FALSE)
})

test_that("Nimble Append output Seurat object is valid when creating new assay", {
  seuratObj <- readRDS("../testdata/nimbleTest.rds")
  expectedBarcodes <- c("12345_AAAAAAAAAAAAAAAA", "12345_CAAAAAAAAAAAAAA", "12345_GAAAAAAAAAAAAAAA", "12345_TAAAAAAAAAAAAAAA")
  expectedFeatures <- c("E", "F", "G", "H")
  expectedValues <- list(c(1, 1, 1, 1), c(2, 2, 2, 2), c(4, 4, 4, 4), c(0, 0, 0, 0))
  seuratObj <- AppendNimbleCounts(seuratObj, "../testdata/12345_nimbleCounts.tsv", appendToCounts=FALSE)
  expect_equal(colnames(seuratObj@assays$Nimble@counts), expectedBarcodes)
  expect_equal(rownames(seuratObj@assays$Nimble@counts), expectedFeatures)
  expect_equal(FALSE %in% (seuratObj@assays$Nimble@counts[, expectedBarcodes[1]] == expectedValues[[1]]), FALSE)
  expect_equal(FALSE %in% (seuratObj@assays$Nimble@counts[, expectedBarcodes[2]] == expectedValues[[2]]), FALSE)
  expect_equal(FALSE %in% (seuratObj@assays$Nimble@counts[, expectedBarcodes[3]] == expectedValues[[3]]), FALSE)
  expect_equal(FALSE %in% (seuratObj@assays$Nimble@counts[, expectedBarcodes[4]] == expectedValues[[4]]), FALSE)
})