context("copyToDropbox")

test_that("RDS token", {
    prepareTemplate("bibliography.bib")
    x <- copyToDropbox(
        files = "bibliography.bib",
        dir = file.path("bcbioBase_examples", "copyToDropbox"),
        rdsToken = system.file("token.rds")
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
    unlink("bibliography.bib")
})
