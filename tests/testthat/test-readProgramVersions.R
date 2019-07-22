context("readProgramVersions")

test_that("programs.txt", {
    file <- file.path("cache", "programs.txt")
    versions <- readProgramVersions(file)
    expect_is(versions, "DataFrame")
    expect_identical(
        object = colnames(versions),
        expected = c("program", "version")
    )
})

## Allow missing file, since bcbio doesn't always generate this.
test_that("Missing file", {
    expect_message(
        object = readProgramVersions("XXX.csv"),
        regexp = "Program versions are missing"
    )
    expect_identical(
        object = readProgramVersions("XXX.txt"),
        expected = DataFrame()
    )
})
