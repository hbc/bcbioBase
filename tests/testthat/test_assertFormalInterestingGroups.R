context("assertFormalInterestingGroups")

test_that("Success", {
    expect_silent(
        assertFormalInterestingGroups(mtcars, colnames(mtcars)[1L:2L])
    )
})

test_that("Failure", {
    expect_error(
        assertFormalInterestingGroups(mtcars, interestingGroups = "XXX"),
        paste(
            "is_subset :",
            "The element 'XXX' in interestingGroups is not in",
            "colnames\\(x\\)."
        )
    )
})
