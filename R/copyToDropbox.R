#' Copy Files to Dropbox
#'
#' @family Write Functions
#' @author Michael Steinbaugh, Victor Barerra, John Hutchinson
#'
#' @param files Local file paths.
#' @param dir Relative path of remote Dropbox directory.
#' @param rdsToken RDS file token to use for Dropbox authentication.
#'
#' @return Invisible `list` of rdrop2 output.
#' @export
#'
#' @examples
#' \dontrun{
#' copyToDropbox(files = c("raw_counts.csv", "tpm.csv"), dir = "counts")
#' }
copyToDropbox <- function(
    files,
    dir,
    rdsToken = NULL
) {
    assert_is_any_of(files, c("character", "list"))
    assert_all_are_existing_files(files)
    assert_is_a_string(dir)
    dir <- gsub("/$", "", dir)
    assert_all_are_non_missing_nor_empty_character(dir)
    if (is.character(rdsToken)) {
        assert_is_a_string(rdsToken)
        assert_all_are_existing_files(rdsToken)
    } else {
        rdsToken <- NA  # nocov
    }

    # Ensure user is authenticated with Dropbox
    drop_auth(rdstoken = rdsToken)

    # Display account information
    acc <- drop_acc()

    message(paste(
        "Dropbox:",
        acc[["name"]][["display_name"]],
        paste0("<", acc[["email"]], ">")
    ))

    # Dropbox output directory
    if (!suppressWarnings(drop_exists(dir))) {
        drop_create(dir)  # nocov
    }

    # Warn about writes into shared directories
    metadata <- drop_get_metadata(dir)
    if (any(
        c("parent_shared_folder_id", "sharing_info") %in% names(metadata)
    )) {
        warning(paste(
            "rdrop2 currently isn't working well with shared directories.",
            "For the time being, please write to an unshared directory.",
            "The files can be then moved manually on your Dropbox account",
            "and the link URLs will be preserved."
        ))
    }

    # Loop across the files in list
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
