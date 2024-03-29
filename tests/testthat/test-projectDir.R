test_that("Single dated directory (standard)", {
    expect_identical(
        object = projectDir(uploadDir),
        expected = file.path(uploadDir, "2018-01-01_bcbio")
    )
})

test_that("Multiple dated directories", {
    uploadDir <- tempdir2()
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
    unlink2(uploadDir)
})
