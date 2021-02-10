#' Project directory
#'
#' @note This will pick the latest dated directory and warn the user, if bcbio
#' has been run multiple times to the same upload directory.
#'
#' @author Michael Steinbaugh
#' @note Updated 2020-12-03.
#' @export
#'
#' @inheritParams AcidRoxygen::params
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
    uploadDir <- realpath(uploadDir)
    dir <- sort(list.files(
        path = uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    ))
    if (!hasLength(dir)) {
        stop(sprintf(
            "Failed to locate dated bcbio project directory in '%s'.",
            uploadDir
        ))
    }
    ## Check to see if user has run bcbio multiple times to the same upload
    ## directory, and warn when this is detected.
    if (length(dir) > 1L) {
        newest <- tail(dir, n = 1L)
        alertWarning("Multiple project directories detected:")
        ul(dir)
        alert(sprintf("Using most recent: %s", newest))
        dir <- newest
    }
    assert(isString(dir))
    cli_dl(c(projectDir = dir))
    realpath(file.path(uploadDir, dir))
}
