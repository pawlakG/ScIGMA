library(ScIGMA)
test_that("App instantiation", {
    # Avoid actually starting the app in tests, just check if the UI/Server objects exist
    expect_true(exists("app_ui"))
    expect_true(exists("app_server"))
})
