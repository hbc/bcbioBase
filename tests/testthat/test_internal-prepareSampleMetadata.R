context("prepareSampleMetadata")

test_that("Missing description column", {
    expect_error(
        .prepareSampleMetadata(mtcars),
        "`description` column is required"
    )
})
