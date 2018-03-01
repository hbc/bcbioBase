context("readProgramVersions")

test_that("readProgramVersions", {
    versions <- readProgramVersions(
        paste(
            "http://bcbiobase.seq.cloud",
            "bcbio",
            "programs.txt",
            sep = "/"))
    expect_is(versions, "tbl_df")
    expect_identical(
        colnames(versions),
        c("program", "version")
    )
})

test_that("Missing file", {
    expect_warning(
        readProgramVersions("XXX.txt"),
        "is_existing_file :"
    )
    expect_identical(
        suppressWarnings(
            readProgramVersions("XXX.txt")
        ),
        NULL
    )
})
