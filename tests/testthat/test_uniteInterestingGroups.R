context("uniteInterestingGroups")

test_that("mtcars", {
    data <- uniteInterestingGroups(
        mtcars,
        interestingGroups = c("gear", "carb")
    )
    expect_identical(
        head(data[["interestingGroups"]]),
        factor(
            c("4:4", "4:4", "4:1", "3:1", "3:2", "3:1"),
            levels = c("3:1", "3:2", "3:3", "3:4",
                       "4:1", "4:2", "4:4",
                       "5:2", "5:4", "5:6", "5:8")
        )
    )
    expect_error(
        uniteInterestingGroups(
            mtcars,
            interestingGroups = c("XXX", "YYY")
        ),
        "Interesting groups not defined in metadata: XXX, YYY"
    )
})
