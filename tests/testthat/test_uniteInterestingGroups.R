context("uniteInterestingGroups")

test_that("Single interesting group", {
    data <- uniteInterestingGroups(mtcars, interestingGroups = "gear")
    expect_identical(
        data[["interestingGroups"]],
        as.factor(mtcars[["gear"]])
    )
})

test_that("Two interesting groups", {
    data <- uniteInterestingGroups(
        mtcars,
        interestingGroups = c("gear", "carb"))
    expect_identical(
        head(data[["interestingGroups"]]),
        factor(
            c("4:4", "4:4", "4:1", "3:1", "3:2", "3:1"),
            levels = c("3:1", "3:2", "3:3", "3:4",
                       "4:1", "4:2", "4:4",
                       "5:2", "5:4", "5:6", "5:8")
        )
    )
})

test_that("Missing interesting group", {
    expect_error(
        uniteInterestingGroups(mtcars, interestingGroups = c("XXX", "YYY")),
        paste(
            "is_subset :",
            "The elements 'XXX', 'YYY' in interestingGroups are not in",
            "colnames\\(x\\)."
        )
    )
})
