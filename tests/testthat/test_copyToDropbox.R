context("copyToDropbox")

prepareTemplate("bibliography.bib")
dropboxDir <- file.path("bcbioBase_examples", "copyToDropbox")

test_that("RDS token", {
    x <- copyToDropbox(
        files = "bibliography.bib",
        dir = dropboxDir,
        rdsToken = system.file("token.rds", package = "bcbioBase")
    )
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
            files = "bibliography.bib",
            dir = dropboxDir,
            rdsToken = "XXX.rds"),
        "XXX.rds does not exist"
    )
})

test_that("HBC Team Folder", {
    expect_error(
        copyToDropbox(
            files = "bibliography.bib",
            dir = "HBC Team Folder (1)"),
        "rdrop2 is not detecting directories correctly in the HBC Team Folder"
    )
})

unlink("bibliography.bib")
