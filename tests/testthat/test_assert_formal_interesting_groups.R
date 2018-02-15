context("assert_formal_interesting_groups")

test_that("Success", {
    assert_formal_interesting_groups(mtcars, colnames(mtcars)[1L:2L])
})

test_that("Failure", {
    expect_error(
        assert_formal_interesting_groups(mtcars, interestingGroups = "XXX"),
        paste(
            "is_subset :",
            "The element 'XXX' in interestingGroups is not in",
            "colnames\\(object\\)."
        )
    )
})
