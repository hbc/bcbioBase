context("copyToDropbox")

## Testing locally with Dropbox authentication token.
## This is a pain to set up on Travis CI and AppVeyor.
skip_if_not(file.exists("token.rds"))

files <- c("demultiplexed.csv", "multiplexed.csv")
dropboxDir <- file.path("bcbioBase_examples", "copyToDropbox")

test_that("RDS token enabled", {
    object <- copyToDropbox(
        files = files,
        dir = dropboxDir,
        rdsToken = "token.rds"
    )
    expect_is(object, "list")
    expect_identical(
        object = lapply(object[[1L]], class),
        expected = list(
            ".tag" = "character",
            url = "character",
            id = "character",
            name = "character",
            "path_lower" = "character",
            "link_permissions" = "list",
            "preview_type" = "character",
            "client_modified" = "character",
            "server_modified" = "character",
            rev = "character",
            size = "integer"
        )
    )
})

test_that("Shared Dropbox directory", {
    expect_warning(
        object = copyToDropbox(
            files = files,
            dir = paste0(dropboxDir, "_shared"),
            rdsToken = "token.rds"
        ),
        regexp = "rdrop2 currently isn't working well with shared"
    )
    ## Don't clean up directory, because we won't be able to check if shared.
})

test_that("Invalid parameters", {
    expect_error(
        copyToDropbox(files = NULL, dir = ".")
    )
    expect_error(
        copyToDropbox(files = "XXX.csv.gz", dir = ".")
    )
    expect_error(
        object = copyToDropbox(files = files, dir = NULL)
    )
    expect_error(
        copyToDropbox(files = files, dir = dropboxDir, rdsToken = "XXX.rds")
    )
})
