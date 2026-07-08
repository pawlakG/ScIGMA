library(testthat)
library(ScIGMA)

test_that("CNV processing generates accurate ploidy matrix and clonal profiles", {
  # Setup mock
  mock_file <- tempfile(fileext = ".h5")
  create_mock_h5(mock_file, n_cells = 20, n_vars = 10, n_cnv = 5, n_ptn = 2, include_protein = FALSE)
  
  # Inject amplicon mapping into the mock H5 file for the variants
  rhdf5::h5write(rep(paste0("AMPLICON_", 1:5), each = 2), mock_file, "assays/dna_variants/ca/amplicon")
  
  on.exit(if (file.exists(mock_file)) file.remove(mock_file), add = TRUE)
  
  scigma_obj <- ScIGMA:::loadH5_HDF5_biocond(mock_file, sample_name = "test_sample", omic_type = "DNA")
  
  cell_ids <- colnames(scigma_obj$mae[["dna_variants"]])
  
  # Inject known CNV read counts
  # 5 amplicons x 20 cells
  # Let's say all cells have 100 reads per amplicon, except cells 1-10 which have 50 reads for amplicon 1
  cnv_counts <- matrix(100, nrow = 5, ncol = 20)
  cnv_counts[1, 1:10] <- 50
  
  # Setup DP counts
  dp_counts <- matrix(100, nrow = 10, ncol = 20)
  dp_counts[, 20] <- 5 # depth 5 < amp_readDepth=10
  
  scigma_obj$mae[["dna_variants"]]@assays@data$dp <- as.matrix(dp_counts)
  scigma_obj$mae[["amplicons"]]@assays@data$counts <- as.matrix(cnv_counts)
  
  # Set some clone assignments manually
  scigma_obj$dna.clones <- as.factor(c(rep("Clone 1", 10), rep("WT", 10)))
  names(scigma_obj$dna.clones) <- cell_ids
  
  # 1. Test filtering function
  filtered_data <- ScIGMA:::filter_cnv_matrix_by_completeness(
    scigma_obj,
    amp_completeness = 50,
    amp_readDepth = 10,
    amp_meanCellRead = 10
  )
  
  # Cell 20 should be dropped since its depth is 5 < 10 for all amplicons
  expect_equal(ncol(filtered_data$filtered_cnv_mtx), 19)
  expect_false(cell_ids[20] %in% colnames(filtered_data$filtered_cnv_mtx))
  
  # 2. Test profile generation
  clonal_profile <- ScIGMA:::process_cnv_to_clonal_profile(
    filtered_data = filtered_data,
    dna_variant_clones = scigma_obj$dna.clones,
    diploid_ref = "WT",
    exclude_clone = "small"
  )
  
  # We should have 2 clones (Clone 1 and WT) and 5 amplicons
  expect_equal(ncol(clonal_profile), 2)
  expect_equal(nrow(clonal_profile), 5)
  expect_true(all(c("Clone 1", "WT") %in% colnames(clonal_profile)))
  
  # The normalization should make the WT profile exactly 2 (diploid ref)
  expect_true(all(clonal_profile[, "WT"] == 2))
  
  # Ratio Clone 1 / WT for Amp 1: 50 / 100 = 0.5. 0.5 * 2 = 1.
  expect_equal(clonal_profile[1, "Clone 1"], 1)
  
  # Ratio Clone 1 / WT for Amp 2: 100 / 100 = 1. 1 * 2 = 2.
  expect_equal(clonal_profile[2, "Clone 1"], 2)
})
