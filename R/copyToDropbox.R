#' Copy Files to Dropbox
#'
#' @author Michael Steinbaugh, Victor Barerra, John Hutchinson
#'
#' @importFrom rdrop2 drop_auth drop_delete drop_exists drop_share drop_upload
#'
#' @param files Local file paths.
#' @param dir Relative path of remote Dropbox directory.
#'
#' @return Invisibly return [list] of rdrop2 output.
#' @export
#' 
#' @examples
#' \dontrun{
#' files <- c("raw_counts.csv.gz", "tpm.csv.gz")
#' copyToDropbox(
#'     files,
#'     dir = file.path("researcher", "project", "results")
#' )
#' }
copyToDropbox <- function(files, dir) {
    if (!(is.character(files) || is.list(files))) {
        abort("`files` must be a character vector or list")
    }
    if (!is_string(dir)) {
        abort("`dir` must be a string")
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
    drop_auth()

    # Loop across the files in list
    rdrop <- lapply(files, function(file) {
        dropboxFile <- file.path(dir, basename(file))
        # Delete file if it exists already. Otherwise `drop_share()`
        # is currently erroring out when we try to re-share an existing
        # file.
        if (drop_exists(dropboxFile)) {
            drop_delete(dropboxFile)
        }
        drop_upload(file = file, path = dir)
        drop_share(dropboxFile, requested_visibility = "public")
    })

    invisible(rdrop)
}
