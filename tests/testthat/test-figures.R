test_that("Figure generation returns valid plots", {
    # Generate mock data
    set.seed(42)
    mock_protein_counts <- matrix(rnorm(500, mean = 5), nrow = 10, ncol = 50)
    rownames(mock_protein_counts) <- paste0("P", 1:10)
    colnames(mock_protein_counts) <- paste0("Cell", 1:50)
    
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(
            proteins = SummarizedExperiment::SummarizedExperiment(
                assays = list(counts = mock_protein_counts, normalized = mock_protein_counts)
            )
        )
    )
    
    # Add dummy UMAP to colData
    SummarizedExperiment::colData(mae[["proteins"]])$UMAP_1 <- runif(50)
    SummarizedExperiment::colData(mae[["proteins"]])$UMAP_2 <- runif(50)
    
    mock_obj <- list(mae = mae)
    
    # Test generate_protein_umap
    p_umap <- generate_protein_umap(mock_obj)
    expect_s3_class(p_umap, "plotly")
    
    # Test generate_protein_biplot
    p_biplot <- generate_protein_biplot(mock_obj, xvar = "P1", yvar = "P2")
    expect_s3_class(p_biplot, "plotly")
    
    # Test generate_protein_barplot
    p_bar <- generate_protein_barplot(mock_obj)
    expect_s3_class(p_bar, "plotly")
})
