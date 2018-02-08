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
#'     rdsToken = system.file("token.rds")
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
    # Ensure trailing slash gets stripped
    dir <- gsub("/$", "", dir)
    # rdsToken
    if (!(is.null(rdsToken) || is_string(rdsToken))) {
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
        rdsToken <- NA
    }
    drop_auth(rdstoken = rdsToken)

    # Check that the desired output directory exists on Dropbox
    recursive <- unlist(strsplit(dir, .Platform[["file.sep"]]))
    # Using `suppressWarnings()` here to avoid:
    # Unknown or uninitialised column: 'path_display'.
    invisible(lapply(
        X = seq_along(recursive),
        FUN = function(a) {
            path <- paste(recursive[1L:a], collapse = .Platform[["file.sep"]])
            if (!suppressWarnings(drop_exists(path))) {
                suppressWarnings(drop_create(path))
            }
        }
    ))
    if (!suppressWarnings(drop_exists(dir))) {
        abort("Failed to create destination directory")
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
