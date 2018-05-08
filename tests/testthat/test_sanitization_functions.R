context("Sanitization Functions")



# sanitizeSampleData ===========================================================
test_that("sanitizeSampleData", {
    x <- sanitizeSampleData(sd)
    expect_is(x, "DataFrame")
    expect_identical(rownames(x), rownames(sd))
    expect_true(all(vapply(
        X = x,
        FUN = is.factor,
        FUN.VALUE = logical(1L)
    )))
})
