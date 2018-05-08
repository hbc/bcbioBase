context("Sanitization Functions")



# sampleData
sd <- DataFrame(
    "genotype" = factor(c("wt", "ko", "wt", "ko")),
    "batch" = factor(c(1L, 1L, 2L, 2L)),
    # not a factor yet
    "day" = c(14L, 14L, 30L, 30L),
    row.names = c("sample_1", "sample_2", "sample_3", "sample_4")
)



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
