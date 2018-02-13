context("copyToDropbox")

prepareTemplate("bibliography.bib")
files <- "bibliography.bib"
dropboxDir <- file.path("bcbioBase_examples", "copyToDropbox")
rdsToken <- system.file("token.rds", package = "bcbioBase")
stopifnot(file.exists(rdsToken))

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
    expect_error(
        copyToDropbox(files = NULL, dir = getwd()),
        paste(
            "is2 :",
            "files is not in any of the classes 'character', 'list'."
        )
    )
    expect_error(
        copyToDropbox(files = "XXX.csv.gz", dir = getwd()),
        paste(
            "is_existing_file :",
            "Some or all of the files specified by files do not exist."
        )
    )
    expect_error(
        copyToDropbox(files = "bibliography.bib", dir = NULL),
        paste(
            "is_a_string :",
            "dir is not of class 'character'; it has class 'NULL'"
        )
    )
    expect_error(
        copyToDropbox(
            files = files,
            dir = dropboxDir,
            rdsToken = mtcars),
        "is2 : rdsToken is not in any of the classes 'character', 'logical'."
    )
    expect_error(
        copyToDropbox(
            files = files,
            dir = dropboxDir,
            rdsToken = "XXX.rds"),
        paste(
            "is_existing_file :",
            "Some or all of the files specified by rdsToken do not exist."
        )
    )
})

test_that("Shared directory", {
    expect_warning(
        copyToDropbox(
            files = files,
            dir = paste0(dropboxDir, "_shared"),
            rdsToken = rdsToken),
        "rdrop2 currently isn't working well with shared directories."
    )
    # Don't unlink this directory, because we won't be able to check if shared
})

unlink("bibliography.bib")
rdrop2::drop_delete(dropboxDir)
