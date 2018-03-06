context("prepareTemplate")

test_that("All default shared files", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "_setup.R",
        "bibliography.bib"
    )
    expect_silent(prepareTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})

test_that("Single file", {
    prepareTemplate("bibliography.bib")
    expect_true(file.exists("bibliography.bib"))
    unlink("bibliography.bib")
})

test_that("Missing source file", {
    expect_error(
        prepareTemplate("XXX.R"),
        "is_existing_file :"
    )
})
