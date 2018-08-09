context("Data Functions")



# minimalSampleData ============================================================
test_that("minimalSampleData", {
    expect_identical(
        minimalSampleData(c("sample 1", "sample 2")),
        data.frame(
            sampleName = c("sample 1", "sample 2"),
            description = c("sample 1", "sample 2"),
            row.names = c("sample_1", "sample_2"),
            stringsAsFactors = TRUE
        )
    )
})



# sampleDirs ===================================================================
test_that("sampleDirs", {
    uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
    expect_identical(
        sampleDirs(uploadDir),
        c(
            "sample_1" = file.path(uploadDir, "sample_1"),
            "sample_2" = file.path(uploadDir, "sample_2")
        )
    )
})
