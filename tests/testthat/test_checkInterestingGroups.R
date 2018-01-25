context("checkInterestingGroups")

test_that("Missing interesting groups", {
    expect_error(
        checkInterestingGroups(mtcars, interestingGroups = "XXX"),
        "Interesting groups not defined in metadata: XXX"
    )
})

test_that("Warn on NULL mode", {
    expect_identical(
        checkInterestingGroups(
            mtcars,
            interestingGroups = NULL,
            warnOnNULL = FALSE),
        "sampleName"
    )
    expect_warning(
        checkInterestingGroups(
            mtcars,
            interestingGroups = NULL,
            warnOnNULL = TRUE),
        "Defaulting to `sampleName`"
    )
})
