test_that("Edge cases are handled appropriately", {
    # Empty drops scenarios
    # Mock empty dataset
    mae_empty <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(
            proteins = SummarizedExperiment::SummarizedExperiment(
                assays = list(counts = matrix(0, nrow = 10, ncol = 0))
            )
        )
    )
    mock_obj <- list(mae = mae_empty)
    
    expect_error(normalize_protein_clr(mock_obj), NA) # Should handle gracefully or throw informative error
})
