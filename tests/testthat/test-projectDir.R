context("projectDir")

test_that("Single dated directory (standard)", {
    expect_identical(
        object = projectDir(uploadDir),
        expected = file.path(uploadDir, "2018-01-01_bcbio")
    )
})

test_that("Multiple dated directories", {
    uploadDir <- "XXX"
    unlink(uploadDir, recursive = TRUE)
    dir.create(uploadDir)
    uploadDir <- realpath(uploadDir)
    dir.create(file.path(uploadDir, "2018-01-01_rnaseq"))
    dir.create(file.path(uploadDir, "2018-02-01_rnaseq"))
    expect_message(
        object = projectDir(uploadDir),
        regexp = "Multiple project directories detected"
    )
    object <- suppressWarnings(projectDir(uploadDir))
    expect_identical(
        object = object,
        expected = file.path(uploadDir, "2018-02-01_rnaseq")
    )
    unlink(uploadDir, recursive = TRUE)
})
