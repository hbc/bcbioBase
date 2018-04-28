context("Assert Check Functions")



test_that("assertFormalInterestingGroups", {
    expect_silent(
        assertFormalInterestingGroups(
            x = rse_bcb,
            interestingGroups = c("tissue", "treatment")
        )
    )
    # Must exist as columns in sampleData
    expect_error(
        assertFormalInterestingGroups(
            x = rse_bcb,
            interestingGroups = "XXX"
        ),
        paste(
            "The interesting groups \"XXX\" are not defined"
        )
    )
    # Require interesting groups to be defined as factor columns
    expect_error(
        assertFormalInterestingGroups(
            x = rse_bcb,
            interestingGroups = c("totalReads", "exonicRate")
        ),
        "The interesting groups \"totalReads, exonicRate\" are not factor"
    )
    # Error on blacklisted columns
    expect_error(
        assertFormalInterestingGroups(
            x = rse_bcb,
            interestingGroups = "genomeBuild"
        ),
        "The interesting groups \"genomeBuild\" are blacklisted."
    )
})



test_that("assertFormalAnnotationCol", {
    x <- assay(rse_dds)
    y <- sampleData(rse_dds)
    expect_silent(assertFormalAnnotationCol(x, y))
    expect_silent(assertFormalAnnotationCol(x, NA))
    expect_silent(assertFormalAnnotationCol(x, NULL))
})
