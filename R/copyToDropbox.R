#' Copy Files to Dropbox
#'
#' @author Michael Steinbaugh, Victor Barerra, John Hutchinson
#'
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#'
#' @param files Local file paths.
#' @param dir Relative path of remote Dropbox directory.
#' @param rdsToken RDS file token to use for Dropbox authentication.
#'
#' @return Invisibly return [list] of rdrop2 output.
#' @export
#' 
#' @examples
#' prepareTemplate("bibliography.bib")
#' copyToDropbox(
#'     files = "bibliography.bib",
#'     dir = file.path("bcbioBase_examples", "copyToDropbox"),
#'     rdsToken = system.file("token.rds", package = "bcbioBase")
#' )
#' unlink("bibliography.bib")
copyToDropbox <- function(
    files,
    dir,
    rdsToken = NA) {
    # files
    if (!(is.character(files) || is.list(files))) {
        abort("`files` must be a character vector or list")
    }
    # dir
    if (!is_string(dir)) {
        abort("`dir` must be a string")
    }
    dir <- gsub("/$", "", dir)
    # rdsToken
    if (!(is_string(rdsToken) || is.na(rdsToken))) {
        abort("`rdsToken` must contain an RDS file or NA")
    }

    # Check that local files exist
    if (!all(vapply(files, file.exists, logical(1L)))) {
        missing <- !vapply(files, file.exists, logical(1L))
        abort(paste(
            "Missing local files:",
            toString(basename(files[missing]))
        ))
    }

    # Ensure user is authenticated with Dropbox
    if (is_string(rdsToken)) {
        if (!file.exists(rdsToken)) {
            abort(paste(rdsToken, "does not exist"))
        }
    }
    drop_auth(rdstoken = rdsToken)

    # Display account information
    acc <- drop_acc()
    inform(paste(
        "Dropbox:",
        acc$name$display_name,
        paste0("<", acc$email, ">")
    ))
    
    # Dropbox output directory
    if (!drop_exists(dir)) {
        drop_create(dir)
    }
    
    # Warn about writes into shared directories
    metadata <- drop_get_metadata(dir)
    if (any(
        c("parent_shared_folder_id", "sharing_info") %in% names(metadata)
    )) {
        warn(paste(
            "rdrop2 currently isn't working well with shared directories.",
            "For the time being, please write to an unshared directory.",
            "The files can be then moved manually on your Dropbox account",
            "and the link URLs will be preserved."
        ))
    }
    
    # Loop across the files in list
    rdrop <- lapply(files, function(file) {
        dropboxFile <- file.path(dir, basename(file))
        if (drop_exists(dropboxFile)) {
            drop_delete(dropboxFile)
        }
        drop_upload(file = file, path = dir)
        drop_share(dropboxFile, requested_visibility = "public")
    })

    invisible(rdrop)
}
