library(testthat)
library(ScIGMA)

test_that("extract_variant_genotypes correctly extracts and formats genotypes", {
  # Setup mock
  mock_file <- tempfile(fileext = ".h5")
  create_mock_h5(mock_file, n_cells = 10, n_vars = 5, n_cnv = 3, n_ptn = 2, include_protein = TRUE)
  on.exit(if (file.exists(mock_file)) file.remove(mock_file), add = TRUE)
  
  # Load it
  scigma_obj <- ScIGMA:::loadH5_HDF5_biocond(mock_file, sample_name = "test", omic_type = "DNA+protein")
  
  # Inject specific genotypes for a variant
  gt_matrix <- matrix(3, nrow = 5, ncol = 10)
  gt_matrix[1, ] <- c(0, 1, 2, 3, 0, 1, 2, 3, 0, 1) # WT, HET, HOM, Missing
  rownames(gt_matrix) <- rownames(scigma_obj$mae[["dna_variants"]])
  scigma_obj$mae[["dna_variants"]]@assays@data$gt <- gt_matrix
  
  target_variant <- rownames(gt_matrix)[1]
  
  # Extract
  res <- ScIGMA:::extract_variant_genotypes(scigma_obj$mae, target_variant, use_compass = FALSE)
  
  # Check format
  expect_equal(nrow(res), 10)
  expect_true(all(c("Barcode", "Variant_Genotype") %in% colnames(res)))
  
  # Check decoding
  expect_equal(res$Variant_Genotype[1], "WT")
  expect_equal(res$Variant_Genotype[2], "HET")
  expect_equal(res$Variant_Genotype[3], "HOM")
  expect_equal(res$Variant_Genotype[4], "Missing/ADO")
})

test_that("compute_population_genotype_distribution computes accurate stats", {
  # Setup mock
  mock_file <- tempfile(fileext = ".h5")
  create_mock_h5(mock_file, n_cells = 10, n_vars = 5, n_cnv = 3, n_ptn = 2, include_protein = TRUE)
  on.exit(if (file.exists(mock_file)) file.remove(mock_file), add = TRUE)
  
  # Load it
  scigma_obj <- ScIGMA:::loadH5_HDF5_biocond(mock_file, sample_name = "test", omic_type = "DNA+protein")
  
  # Inject specific genotypes
  gt_matrix <- matrix(3, nrow = 5, ncol = 10)
  gt_matrix[1, ] <- c(0, 0, 1, 1, 2, 2, 3, 3, 0, 0)
  rownames(gt_matrix) <- rownames(scigma_obj$mae[["dna_variants"]])
  scigma_obj$mae[["dna_variants"]]@assays@data$gt <- gt_matrix
  
  target_variant <- rownames(gt_matrix)[1]
  barcodes <- colnames(scigma_obj$mae[["dna_variants"]])
  
  # 1. No cluster provided
  res <- ScIGMA:::compute_population_genotype_distribution(
    scigma_obj$mae,
    target_variant,
    barcodes,
    use_compass = FALSE
  )
  
  expect_equal(nrow(res), 4) # 4 genotype classes
  wt_row <- res[res$Variant_Genotype == "WT", ]
  expect_equal(wt_row$Count, 4)
  expect_equal(wt_row$Percentage, 40)
  
  # 2. Cluster provided
  clusters <- rep(c("C1", "C2"), each = 5)
  names(clusters) <- barcodes
  
  res_clust <- ScIGMA:::compute_population_genotype_distribution(
    scigma_obj$mae,
    target_variant,
    barcodes,
    use_compass = FALSE,
    seurat_cluster = clusters
  )
  
  expect_equal(nrow(res_clust), 6) # C1 has 3 genotypes, C2 has 3 genotypes
  
  c1_wt <- res_clust[res_clust$Cluster == "C1" & res_clust$Variant_Genotype == "WT", ]
  expect_equal(c1_wt$Count, 2)
  expect_equal(c1_wt$Percentage, 40) # 2 out of 5 in C1
})
