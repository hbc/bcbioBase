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
    # Error on "HBC Team Folder" detection
    # TODO This should be safe to remove in a future update
    if (grepl("HBC Team Folder", dir)) {
        abort(paste(
            "rdrop2 is not detecting directories correctly",
            "in the HBC Team Folder share at the moment.",
            "Please save the files elsewhere on Dropbox and then",
            "move them manually. The links will still work."
        ))
    }
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
    if (!drop_exists(dir)) {
        drop_create(dir)
        # Double check to ensure creation was successful
        if (!drop_exists(dir)) {
            abort("rdrop2 is failing to create and detect output directory")
        }
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
