context("Assert Check Functions")



test_that("assertFormalInterestingGroups", {
    expect_silent(
        assertFormalInterestingGroups(mtcars, colnames(mtcars)[1L:2L])
    )
    expect_error(
        assertFormalInterestingGroups(mtcars, interestingGroups = "XXX"),
        paste(
            "is_subset :",
            "The element 'XXX' in interestingGroups is not in",
            "colnames\\(x\\)."
        )
    )
})



test_that("assertFormalAnnotationCol", {
    x <- assay(rse_dds)
    y <- sampleData(rse_dds)
    expect_silent(assertFormalAnnotationCol(x, y))
    expect_silent(assertFormalAnnotationCol(x, NA))
    expect_silent(assertFormalAnnotationCol(x, NULL))
})
