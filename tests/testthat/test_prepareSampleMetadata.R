context("prepareSampleMetadata")

test_that("Missing description column", {
    expect_error(
        .prepareSampleMetadata(mtcars),
        paste(
            "is_subset :",
            "The element 'description' in \"description\" is not in",
            "colnames\\(object\\)."
        )
    )
})
