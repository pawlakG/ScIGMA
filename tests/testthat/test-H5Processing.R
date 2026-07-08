library(testthat)
library(ScIGMA)
library(rhdf5)
library(MultiAssayExperiment)

test_that("loadH5_HDF5_biocond properly loads mocked HDF5 data into MAE", {
  # Setup mock environment
  mock_file <- tempfile(fileext = ".h5")
  create_mock_h5(mock_file, n_cells = 10, n_vars = 5, n_cnv = 3, n_ptn = 2, include_protein = TRUE)
  
  # Ensure cleanup
  on.exit(if (file.exists(mock_file)) file.remove(mock_file), add = TRUE)
  
  # Test the load function
  scigma_obj <- loadH5_HDF5_biocond(mock_file, sample_name = "test_sample", omic_type = "DNA+protein")
  
  # Structure Assertions
  expect_true(inherits(scigma_obj, "ScIGMA_object"))
  expect_true(inherits(scigma_obj$mae, "MultiAssayExperiment"))
  
  # Dimensions Assertions
  expect_equal(nrow(MultiAssayExperiment::colData(scigma_obj$mae)), 10)
  
  # Assays existence
  assay_names <- names(scigma_obj$mae)
  expect_true(all(c("dna_variants", "amplicons", "proteins") %in% assay_names))
  
  # DNA variants check
  expect_equal(nrow(scigma_obj$mae[["dna_variants"]]), 5)
  expect_true("vaf" %in% SummarizedExperiment::assayNames(scigma_obj$mae[["dna_variants"]]))
  
  # CNV check
  expect_equal(nrow(scigma_obj$mae[["amplicons"]]), 3)
  
  # Protein check
  expect_equal(nrow(scigma_obj$mae[["proteins"]]), 2)
  expect_true("counts" %in% SummarizedExperiment::assayNames(scigma_obj$mae[["proteins"]]))
  
  # Lazy Loading validation (DelayedArray check)
  af_assay <- SummarizedExperiment::assay(scigma_obj$mae[["dna_variants"]], "vaf")
  expect_true(inherits(af_assay, "DelayedMatrix"))
})

test_that("loadH5_HDF5_biocond handles DNA only files correctly", {
  mock_file <- tempfile(fileext = ".h5")
  create_mock_h5(mock_file, n_cells = 10, n_vars = 5, n_cnv = 3, n_ptn = 0, include_protein = FALSE)
  on.exit(if (file.exists(mock_file)) file.remove(mock_file), add = TRUE)
  
  # If the user selects DNA only
  scigma_obj <- loadH5_HDF5_biocond(mock_file, sample_name = "test_dna", omic_type = "DNA")
  
  expect_false("proteins" %in% names(scigma_obj$mae))
  expect_equal(scigma_obj$filetype, "DNA")
})

test_that("loadH5_HDF5_biocond fails gracefully on missing files", {
  expect_error(loadH5_HDF5_biocond("non_existent_file.h5", "test", "DNA"))
})
