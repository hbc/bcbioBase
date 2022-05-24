test_that("sampleDirs", {
    expect_identical(
        object = sampleDirs(uploadDir),
        expected = c(
            "sample_1" = file.path(uploadDir, "sample_1"),
            "sample_2" = file.path(uploadDir, "sample_2")
        )
    )
})
