# UPDATED: tests/testthat/test_integration_scigma.R
library(testthat)
library(MultiAssayExperiment) # NEW
library(SingleCellExperiment)
library(DelayedArray)
library(SummarizedExperiment)

# Setup path: adjust according to your testing directory structure
h5_test_file <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"

test_that("loadH5_HDF5_biocond_constructs_valid_multiomics_mae", {
    skip_if_not(file.exists(h5_test_file), "Integration h5 file not found.")

    # Execution
    scigma_obj <- loadH5_HDF5_biocond(
        filepath = h5_test_file,
        sample_name = "aml_4_lines",
        omic_type = "DNA+protein"
    )

    # Core architecture validation
    expect_s4_class(scigma_obj$mae, "MultiAssayExperiment")
    expect_equal(metadata(scigma_obj$mae)$name, "aml_4_lines")

    # MAE experiments topology validation (No more altExp)
    exp_names <- names(scigma_obj$mae)
    expect_true(all(c("dna_variants", "amplicons", "proteins") %in% exp_names))

    # DNA main experiment validation (Accessed via list subsetting)
    main_assays <- assayNames(scigma_obj$mae[["dna_variants"]])
    expect_true(all(c("vaf", "gt", "dp", "gq") %in% main_assays))
    expect_s4_class(assay(scigma_obj$mae[["dna_variants"]], "vaf"), "DelayedMatrix")

    # Ensure feature spaces are isolated
    dna_features <- ncol(scigma_obj$mae[["dna_variants"]])
    prot_features <- ncol(scigma_obj$mae[["proteins"]])
    expect_true(dna_features > prot_features)
})

