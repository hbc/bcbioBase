context("copyToDropbox")

prepareTemplate("bibliography.bib")
files <- "bibliography.bib"
dropboxDir <- file.path("bcbioBase_examples", "copyToDropbox")
rdsToken <- system.file("token.rds", package = "bcbioBase")

test_that("RDS token", {
    x <- copyToDropbox(
        files = files,
        dir = dropboxDir,
        rdsToken = rdsToken)
    expect_is(x, "list")
    expect_identical(
        lapply(x[[1L]], class),
        list(
            ".tag" = "character",
            "url" = "character",
            "id" = "character",
            "name" = "character",
            "path_lower" = "character",
            "link_permissions" = "list",
            "client_modified" = "character",
            "server_modified" = "character",
            "rev" = "character",
            "size" = "integer"
        )
    )
})

test_that("Invalid parameters", {
    # files
    expect_error(
        copyToDropbox(files = NULL, dir = getwd()),
        "`files` must be a character vector or list"
    )
    expect_error(
        copyToDropbox(files = "XXX.csv.gz", dir = getwd()),
        "Missing local files: XXX.csv.gz"
    )
    # dir
    expect_error(
        copyToDropbox(files = "XXX", dir = NULL),
        "`dir` must be a string"
    )
    # rdsToken
    expect_error(
        copyToDropbox(
            files = files,
            dir = dropboxDir,
            rdsToken = mtcars),
        "`rdsToken` must contain an RDS file or NA"
    )
    expect_error(
        copyToDropbox(
            files = files,
            dir = dropboxDir,
            rdsToken = "XXX.rds"),
        "XXX.rds does not exist"
    )
})

test_that("Shared directory", {
    expect_warning(
        copyToDropbox(
            files = files,
            dir = file.path(dropboxDir, "shared"),
            rdsToken = rdsToken),
        "rdrop2 currently isn't working well with shared directories."
    )
})

unlink("bibliography.bib")
