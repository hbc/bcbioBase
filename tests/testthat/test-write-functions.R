context("Write Functions")



# copyToDropbox ================================================================
if (file.exists("token.rds")) {
    files <- c("demultiplexed.csv", "multiplexed.csv")
    dropboxDir <- file.path("bcbioBase_examples", "copyToDropbox")
    test_that("copyToDropbox : RDS token enabled", {
        x <- copyToDropbox(
            files = files,
            dir = dropboxDir,
            rdsToken = "token.rds"
        )
        expect_is(x, "list")
        expect_identical(
            lapply(x[[1L]], class),
            list(
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

    test_that("copyToDropbox : Shared Dropbox directory", {
        expect_warning(
            copyToDropbox(
                files = files,
                dir = paste0(dropboxDir, "_shared"),
                rdsToken = "token.rds"
            ),
            "rdrop2 currently isn't working well with shared directories."
        )
        # Don't remove directory, because we won't be able to check if shared
    })

    test_that("copyToDropbox : Invalid parameters", {
        expect_error(
            copyToDropbox(files = NULL, dir = "."),
            paste(
                "is2 :",
                "files is not in any of the classes 'character', 'list'."
            )
        )
        expect_error(
            copyToDropbox(files = "XXX.csv.gz", dir = "."),
            paste(
                "is_existing_file :",
                "Some or all of the files specified by files do not exist."
            )
        )
        expect_error(
            copyToDropbox(files = files, dir = NULL),
            paste(
                "is_a_string :",
                "dir is not of class 'character'; it has class 'NULL'"
            )
        )
        expect_error(
            copyToDropbox(
                files = files,
                dir = dropboxDir, rdsToken = "XXX.rds"
            ),
            paste(
                "is_existing_file :",
                "Some or all of the files specified by rdsToken do not exist."
            )
        )
    })
}
