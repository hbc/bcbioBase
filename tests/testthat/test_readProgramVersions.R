context("readProgramVersions")

test_that("readProgramVersions", {
    versions <- readProgramVersions(
        file.path("http://bcbiobase.seq.cloud",
                  "bcbio",
                  "programs.txt"),
        quiet = TRUE)
    expect_true(tibble::is_tibble(versions))
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
