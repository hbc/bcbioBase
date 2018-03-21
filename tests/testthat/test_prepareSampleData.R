context("prepareDataMetadata")

test_that("Missing description column", {
    expect_error(
        prepareSampleData(mtcars),
        paste(
            "is_subset :",
            "The element 'description' in \"description\" is not in",
            "colnames\\(object\\)."
        )
    )
})
