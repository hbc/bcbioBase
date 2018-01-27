context("prepareTemplate")

test_that("All default shared files", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "bibliography.bib",
        "setup.R")
    expect_silent(prepareTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})

test_that("Single file", {
    prepareTemplate("setup.R")
    expect_true(file.exists("setup.R"))
    unlink("setup.R")
})

test_that("Missing source file", {
    expect_error(
        prepareTemplate("XXX.R"),
        "Missing source file"
    )
})
