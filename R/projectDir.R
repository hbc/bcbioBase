#' Project directory
#'
#' @note This will pick the latest dated directory and warn the user, if bcbio
#' has been run multiple times to the same upload directory.
#'
#' @author Michael Steinbaugh
#' @note Updated 2019-08-05.
#' @export
#'
#' @inheritParams acidroxygen::params
#'
#' @return `character(1)`.
#' Dated project directory (e.g. "2018-01-01_rnaseq").
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' x <- projectDir(uploadDir)
#' basename(x)
projectDir <- function(uploadDir) {
    assert(isADirectory(uploadDir))
    dir <- sort(list.files(
        path = uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    ))
    assert(hasLength(dir))
    ## Check to see if user has run bcbio multiple times to the same upload
    ## directory, and warn when this is detected.
    if (length(dir) > 1L) {
        newest <- tail(dir, n = 1L)
        warning(sprintf(
            fmt = paste(
                "Multiple project directories detected:",
                "%s",
                "Using most recent: %s",
                sep = "\n"
            ),
            printString(dir),
            newest
        ))
        dir <- newest
    }
    assert(isString(dir))
    message(sprintf("Dated project directory: %s.", dir))
    realpath(file.path(uploadDir, dir))
}
