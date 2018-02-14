context("readDataVersions")

test_that("readDataVersions", {
    versions <- readDataVersions(
        file.path("http://bcbiobase.seq.cloud",
                  "bcbio",
                  "data_versions.csv"),
        quiet = TRUE)
    expect_true(tibble::is_tibble(versions))
    expect_identical(
        colnames(versions),
        c("genome", "resource", "version")
    )
})

test_that("Missing file", {
    expect_warning(
        readDataVersions("XXX.csv"),
        "XXX.csv missing"
    )
    expect_identical(
        suppressWarnings(
            readDataVersions("XXX.csv")
        ),
        NULL
    )
})
