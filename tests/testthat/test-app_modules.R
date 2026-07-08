test_that("Shiny UI and Server modules initialize properly", {
    # Basic golem tests for UI and Server
    app <- run_app()
    expect_s3_class(app, "shiny.appobj")
})
