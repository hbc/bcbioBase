context("runDate")

test_that("runDate", {
    expect_identical(
        object = runDate(projectDir(uploadDir)),
        expected = as.Date("2018-01-01")
    )
})
