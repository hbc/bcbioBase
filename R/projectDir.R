#' Project Directory
#'
#' @note This will pick the latest dated directory and warn the user, if bcbio
#' has been run multiple times to the same upload directory.
#'
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `string`. Dated project directory (e.g. "2018-01-01_rnaseq").
#' @export
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' x <- projectDir(uploadDir)
#' basename(x)
projectDir <- function(uploadDir) {
    dir <- list.files(
        path = uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    )
    assert_is_non_empty(dir)
    # Check to see if user has run bcbio multiple times to the same upload
    # directory, and warn when this is detected.
    if (length(dir) > 1L) {
        newest <- tail(dir, n = 1L)
        warning(paste(
            "Multiple project directories detected:",
            printString(dir),
            paste("Using most recent:", newest),
            sep = "\n"
        ), call. = FALSE)
        dir <- newest
    }
    assert_is_a_string(dir)
    message(paste("Dated project directory:", dir))
    normalizePath(
        path = file.path(uploadDir, dir),
        winslash = "/",
        mustWork = TRUE
    )
}
