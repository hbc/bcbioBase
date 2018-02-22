context("readDataVersions")

test_that("readDataVersions", {
    versions <- readDataVersions(
        paste("http://bcbiobase.seq.cloud",
              "bcbio",
              "data_versions.csv",
              sep = "/"))
    expect_is(versions, "tbl_df")
    expect_identical(
        colnames(versions),
        c("genome", "resource", "version")
    )
})

test_that("Missing file", {
    expect_warning(
        readDataVersions("XXX.csv"),
        "is_existing_file :"
    )
    expect_identical(
        suppressWarnings(
            readDataVersions("XXX.csv")
        ),
        NULL
    )
})
