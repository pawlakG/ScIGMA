library(testthat)
library(ScIGMA)

test_that("filter_variant_ScIGMA properly filters variants based on thresholds", {
  # Generate mock HDF5 with 20 cells and 10 variants
  mock_file <- tempfile(fileext = ".h5")
  create_mock_h5(mock_file, n_cells = 20, n_vars = 10, n_cnv = 5, n_ptn = 2, include_protein = TRUE)
  on.exit(if (file.exists(mock_file)) file.remove(mock_file), add = TRUE)
  
  # Load the object
  scigma_obj <- ScIGMA:::loadH5_HDF5_biocond(mock_file, sample_name = "test_sample", omic_type = "DNA+protein")
  
  # Inject controlled data into the object to guarantee filtering outcomes
  # We will force some variants to fail DP/GQ filters, and some cells to drop
  
  # ScIGMA_object creates DelayedMatrix by default for vaf.mtx, etc. 
  # But the test replaces them to mock specific values. 
  # Let's instantiate them as matrices first so we can modify them directly.
  scigma_obj$dp.mtx <- as.matrix(scigma_obj$dp.mtx)
  scigma_obj$gq.mtx <- as.matrix(scigma_obj$gq.mtx)
  scigma_obj$gt.mtx <- as.matrix(scigma_obj$gt.mtx)
  scigma_obj$vaf.mtx <- as.matrix(scigma_obj$vaf.mtx)
  
  # Force variant 1 to have low DP across all cells (should be filtered out)
  # DP threshold defaults to 10
  scigma_obj$dp.mtx[, 1] <- 5
  
  # Force variant 2 to have low GQ across all cells (should be filtered out)
  # GQ threshold defaults to 30
  scigma_obj$gq.mtx[, 2] <- 20
  
  # For variant 3, ensure it passes by setting good DP, GQ, and a valid mutated GT
  scigma_obj$dp.mtx[, 3] <- 50
  scigma_obj$gq.mtx[, 3] <- 99
  scigma_obj$gt.mtx[, 3] <- 1  # Heterozygous
  scigma_obj$vaf.mtx[, 3] <- 50
  
  # Force cells 1 and 2 to have NA/low quality for all variants (should be filtered out)
  # Cell filtering requires > min.cell.pt (default 50%) variants to pass
  scigma_obj$dp.mtx[1:2, ] <- 5
  
  # Ensure cells 3-20 have good quality for at least variant 3
  scigma_obj$dp.mtx[3:20, 3] <- 50
  
  # Run filter
  filtered_obj <- ScIGMA:::filter_variant_ScIGMA(
    scigma_obj,
    min.dp = 10,
    min.gq = 30,
    vaf.ref = 5,
    vaf.hom = 95,
    vaf.het = 35,
    min.cell.pt = 50,
    min.mut.cell.pt = 1
  )
  
  # Assertions
  # The original object should not be modified
  expect_equal(nrow(scigma_obj$vaf.mtx), 20)
  expect_equal(ncol(scigma_obj$vaf.mtx), 10)
  
  # The filtered object should have dropped cells 1 and 2
  expect_equal(nrow(filtered_obj$vaf.mtx.filtered), 18)
  
  # The filtered object should have dropped variants 1 and 2
  filtered_variant_ids <- filtered_obj$variants.filtered
  original_variant_ids <- scigma_obj$variants
  
  expect_false(original_variant_ids[1] %in% filtered_variant_ids)
  expect_false(original_variant_ids[2] %in% filtered_variant_ids)
  expect_true(original_variant_ids[3] %in% filtered_variant_ids)
  
  # Check cell IDs
  expect_false(scigma_obj$cell.ids[1] %in% filtered_obj$cell.ids.filtered)
  expect_false(scigma_obj$cell.ids[2] %in% filtered_obj$cell.ids.filtered)
  
  # Check status update
  expect_equal(filtered_obj$variant.filter, "filtered")
})
