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
    expect_true(all(file_exists(files)))
    file_delete(files)
})

test_that("Single file", {
    prepareTemplate("bibliography.bib")
    expect_true(file_exists("bibliography.bib"))
    file_delete("bibliography.bib")
})

test_that("Missing source file", {
    expect_error(
        prepareTemplate("XXX.R"),
        "is_existing_file :"
    )
})
