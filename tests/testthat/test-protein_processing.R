test_that("Protein normalization methods return expected assay structures", {
    # Assuming test mock data is available or can be constructed
    # We will use a mock MultiAssayExperiment
    set.seed(42)
    mock_protein_counts <- matrix(rpois(100, lambda = 10), nrow = 5, ncol = 20)
    rownames(mock_protein_counts) <- paste0("CD", 1:5)
    colnames(mock_protein_counts) <- paste0("Cell", 1:20)
    
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(
            proteins = SummarizedExperiment::SummarizedExperiment(
                assays = list(counts = mock_protein_counts)
            )
        )
    )
    
    mock_obj <- list(mae = mae)
    
    # CLR
    res_clr <- normalize_protein_clr(mock_obj)
    expect_true("normalized" %in% SummarizedExperiment::assayNames(res_clr$mae[["proteins"]]))
    
    # ASINH
    res_asinh <- normalize_protein_asinh(mock_obj, cofactor = 5)
    expect_true("normalized" %in% SummarizedExperiment::assayNames(res_asinh$mae[["proteins"]]))
    
    # NSP
    res_nsp <- normalize_protein_nsp(mock_obj)
    expect_true("normalized" %in% SummarizedExperiment::assayNames(res_nsp$mae[["proteins"]]))
})

test_that("UMAP calculation adds projection to colData", {
    set.seed(42)
    mock_protein_counts <- matrix(rnorm(500, mean = 5), nrow = 10, ncol = 50)
    rownames(mock_protein_counts) <- paste0("P", 1:10)
    colnames(mock_protein_counts) <- paste0("Cell", 1:50)
    
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(
            proteins = SummarizedExperiment::SummarizedExperiment(
                assays = list(normalized = mock_protein_counts)
            )
        )
    )
    
    mock_obj <- list(mae = mae)
    
    # UMAP
    res_umap <- run_umap_protein(mock_protein_counts)
    expect_equal(nrow(res_umap$layout), 50)
    expect_equal(ncol(res_umap$layout), 2)
})
