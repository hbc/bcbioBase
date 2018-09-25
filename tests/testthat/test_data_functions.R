context("Data Functions")

uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")



# projectDir ===================================================================
test_that("projectDir", {
    expect_identical(
        object = projectDir(uploadDir),
        expected = file.path(uploadDir, "2018-01-01_bcbio")
    )
})

test_that("projectDir : Multiple dated directories", {
    uploadDir <- "XXX"
    dir.create(uploadDir)
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    dir.create(file.path(uploadDir, "2018-01-01_rnaseq"))
    dir.create(file.path(uploadDir, "2018-02-01_rnaseq"))
    expect_warning(
        object = projectDir(uploadDir),
        regexp = "Multiple project directories detected"
    )
    object <- suppressWarnings(projectDir(uploadDir))
    expect_identical(
        object = object,
        expected = file.path(uploadDir, "2018-02-01_rnaseq")
    )
    unlink("XXX", recursive = TRUE)
})



# sampleDirs ===================================================================
test_that("sampleDirs", {
    expect_identical(
        object = sampleDirs(uploadDir),
        expected = c(
            "sample_1" = file.path(uploadDir, "sample_1"),
            "sample_2" = file.path(uploadDir, "sample_2")
        )
    )
})
