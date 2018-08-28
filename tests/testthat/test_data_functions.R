context("Data Functions")

uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")



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



# projectDir ===================================================================
test_that("projectDir", {
    expect_identical(
        projectDir(uploadDir),
        file.path(uploadDir, "2018-01-01_bcbio")
    )
})



# sampleDirs ===================================================================
test_that("sampleDirs", {
    expect_identical(
        sampleDirs(uploadDir),
        c(
            "sample_1" = file.path(uploadDir, "sample_1"),
            "sample_2" = file.path(uploadDir, "sample_2")
        )
    )
})
