## nocov start
## Covering this code locally, since dealing with RDS token is buggy.



#' Copy files to Dropbox
#'
#' @author Michael Steinbaugh, Victor Barerra, John Hutchinson
#' @note Updated 2019-08-05.
#' @export
#'
#' @param files `character`.
#'   Local file paths.
#' @param dir `character(1)`.
#'   Relative path of remote Dropbox directory.
#' @param rdsToken `character(1)` or `NULL`.
#'   RDS file token to use for Dropbox authentication.
#'
#' @return Invisible `list`.
#' rdrop2 output, including file paths.
#'
#' @examples
#' # > copyToDropbox(files = c("raw_counts.csv", "tpm.csv"), dir = "counts")
copyToDropbox <- function(
    files,
    dir,
    rdsToken = NULL
) {
    assert(
        allAreFiles(files),
        isADirectory(dir)
    )
    ## rdrop2 has issues with trailing slash, so sanitize.
    dir <- gsub("/$", "", dir)
    if (is.character(rdsToken)) {
        assert(isAFile(rdsToken))
    } else {
        rdsToken <- NA  # nocov
    }

    ## Ensure user is authenticated with Dropbox.
    drop_auth(rdstoken = rdsToken)

    ## Display account information.
    acc <- drop_acc()

    message(sprintf(
        "Dropbox: %s <%s>",
        acc[["name"]][["display_name"]],
        acc[["email"]]
    ))

    ## Dropbox output directory.
    if (!suppressWarnings(drop_exists(dir))) {
        drop_create(dir)  # nocov
    }

    ## Warn about writes into shared directories.
    metadata <- drop_get_metadata(dir)
    if (any(
        c("parent_shared_folder_id", "sharing_info") %in% names(metadata)
    )) {
        warning(
            "rdrop2 doesn't work well with shared directories.\n",
            "For the time being, please write to an unshared directory.\n",
            "The files can be then moved manually on your Dropbox account ",
            "and the link URLs will be preserved."
        )
    }

    ## Loop across the files in list.
    rdrop <- lapply(files, function(file) {
        dropboxFile <- file.path(dir, basename(file))
        if (suppressWarnings(drop_exists(dropboxFile))) {
            drop_delete(dropboxFile)
        }
        drop_upload(file = file, path = dir)
        drop_share(dropboxFile, requested_visibility = "public")
    })

    invisible(rdrop)
}



## nocov end
