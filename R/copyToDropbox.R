#' Copy Files to Dropbox
#'
#' @author Michael Steinbaugh, Victor Barerra, John Hutchinson
#'
#' @importFrom rdrop2 drop_auth drop_create drop_delete drop_exists drop_share
#'   drop_upload
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
    rdsToken = NULL) {
    # files
    if (!(is.character(files) || is.list(files))) {
        abort("`files` must be a character vector or list")
    }
    # dir
    if (!is_string(dir)) {
        abort("`dir` must be a string")
    }
    # Ensure trailing slash gets stripped
    dir <- gsub("/$", "", dir)
    # rdsToken
    if (!(is_string(rdsToken) || is.null(rdsToken))) {
        abort("`rdsToken` must contain an RDS file or NULL")
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
    } else {
        # Match the default in `drop_auth()`
        rdsToken <- NA
    }
    drop_auth(rdstoken = rdsToken)

    # Create Dropbox directory path
    drop_create(dir)
    if (!drop_exists(dir)) {
        abort("rdrop2 failed to create Dropbox directory")
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
