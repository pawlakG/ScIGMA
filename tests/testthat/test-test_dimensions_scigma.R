library(testthat)
library(ScIGMA)
library(MultiAssayExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(DelayedArray)

# Setup path
h5_test_file <- "/Users/geoffrey/Documents/Pro/projects/dnaPtnApp/inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"
file.exists(h5_test_file)

test_that("loadH5_HDF5_biocond_enforces_strict_features_x_cells_paradigm", {
    skip_if_not(file.exists(h5_test_file), "Integration h5 file not found.")
    # print(getwd())
    # Execution
    scigma <- loadH5_HDF5_biocond(
        filepath = h5_test_file,
        sample_name = "aml_4_lines",
        omic_type = "DNA+protein"
    )

    # ---- 1. Validation de l'encapsulation R6 et de l'archive ----
    expect_true(inherits(scigma, "ScIGMA_object"))
    expect_s4_class(scigma$mae, "MultiAssayExperiment")
    expect_s4_class(scigma$mae_raw, "MultiAssayExperiment")

    # Le MAE actif et le MAE brut doivent être strictement identiques à l'import
    expect_equal(dim(scigma$mae), dim(scigma$mae_raw))

    # ---- 2. Validation de l'architecture globale (Le MAE) ----
    mae_obj <- scigma$mae

    # Le nombre global de cellules se lit strictement sur le colData du MAE parent
    global_cells <- nrow(SummarizedExperiment::colData(mae_obj))

    expect_true(global_cells > 0)
    expect_true("dna_sample_name" %in% colnames(SummarizedExperiment::colData(mae_obj)))

    # ---- 3. Validation de l'expérience principale (Variants ADN) ----
    expect_true("dna_variants" %in% names(mae_obj))
    dna_exp <- mae_obj[["dna_variants"]]

    n_variants <- nrow(dna_exp)
    # L'appel à ncol() est valide ici car dna_exp est un SingleCellExperiment
    n_cells_dna <- ncol(dna_exp)

    # L'expérience doit avoir le même nombre de colonnes que le registre global
    expect_equal(n_cells_dna, global_cells)
    expect_true(n_variants > 0)

    # Le rowData doit contenir les annotations des variants (Lignes)
    expect_equal(nrow(SummarizedExperiment::rowData(dna_exp)), n_variants)
    expect_true("chrom" %in% colnames(SummarizedExperiment::rowData(dna_exp)))

    # Validation dimensionnelle intra-matrice (DelayedMatrix)
    vaf_mtx <- SummarizedExperiment::assay(dna_exp, "vaf")
    expect_s4_class(vaf_mtx, "DelayedMatrix")
    expect_equal(nrow(vaf_mtx), n_variants)
    expect_equal(ncol(vaf_mtx), global_cells)

    # ---- 4. Validation de la topologie Multi-Omique (Amplicons) ----
    if ("amplicons" %in% names(mae_obj)) {
        amp_exp <- mae_obj[["amplicons"]]
        n_amplicons <- nrow(amp_exp)

        # Symétrie absolue des cellules exigée
        expect_equal(ncol(amp_exp), global_cells)

        # Validation du rowData des amplicons
        expect_equal(nrow(rowData(amp_exp)), n_amplicons)

        # Biologiquement, le panel d'amplicons ciblé est toujours inférieur aux variants découverts
        expect_true(n_variants > n_amplicons)
    }
})
