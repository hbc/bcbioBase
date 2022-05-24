test_that("programs.txt", {
    file <- file.path("cache", "programs.txt")
    versions <- importProgramVersions(file)
    expect_s4_class(versions, "DataFrame")
    expect_identical(
        object = colnames(versions),
        expected = c("program", "version")
    )
})

## Allow missing file, since bcbio doesn't always generate this.
test_that("Missing file", {
    expect_message(
        object = importProgramVersions("XXX.csv"),
        regexp = "Program versions are missing"
    )
    expect_identical(
        object = importProgramVersions("XXX.txt"),
        expected = DataFrame()
    )
})
