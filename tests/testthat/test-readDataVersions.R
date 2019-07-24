context("readDataVersions")

test_that("readDataVersions", {
    x <- readDataVersions(file.path("cache", "data-versions.csv"))
    expect_is(x, "DataFrame")
    expect_identical(
        object = colnames(x),
        expected = c("genome", "resource", "version")
    )
})

## Allow missing file, since bcbio doesn't always generate this.
## Consider rethinking this approach, and making bcbio stricter?
test_that("Missing file", {
    expect_message(
        object = readDataVersions("XXX.csv"),
        regexp = "Data versions are missing."
    )
    expect_identical(
        object = readDataVersions("XXX.csv"),
        expected = DataFrame()
    )
})
