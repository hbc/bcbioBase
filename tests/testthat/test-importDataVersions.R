test_that("importDataVersions", {
    x <- importDataVersions(file.path("cache", "data-versions.csv"))
    expect_s4_class(x, "DataFrame")
    expect_identical(
        object = colnames(x),
        expected = c("genome", "resource", "version")
    )
})

## Allow missing file, since bcbio doesn't always generate this.
## Consider rethinking this approach, and making bcbio stricter?
test_that("Missing file", {
    expect_message(
        object = importDataVersions("XXX.csv"),
        regexp = "Data versions are missing."
    )
    expect_identical(
        object = importDataVersions("XXX.csv"),
        expected = DataFrame()
    )
})
