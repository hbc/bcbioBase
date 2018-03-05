context("Assert Check Functions")

# assertFormalInterestingGroups ================================================
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
