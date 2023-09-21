#' Run date
#'
#' Get the run date from the project directory.
#'
#' Alternatively, can parse YAML data, but this approach is faster and simpler.
#'
#' @author Michael Steinbaugh
#' @note Updated 2023-09-21.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `Date`.
#'
#' @examples
#' runDate(projectDir = "2018-01-01_illumina_rnaseq")
runDate <- function(projectDir) {
    assert(
        isString(projectDir)
    )
    projectDir <- basename(projectDir)
    assert(isMatchingRegex(projectDir, pattern = projectDirPattern))
    match <- strMatch(x = projectDir, pattern = projectDirPattern)
    as.Date(match[[2L]])
}
