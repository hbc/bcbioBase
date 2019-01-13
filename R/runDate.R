# Alternatively, can parse YAML data, but this is faster and simpler.

#' Run date
#'
#' Get the run date from the project directory.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams params
#'
#' @return `Date`.
#'
#' @examples
#' runDate("2018-01-01_illumina_rnaseq")
runDate <- function(projectDir) {
    assert(isString(projectDir))
    projectDir <- basename(projectDir)
    assert(isMatchingRegex(projectDir, pattern = projectDirPattern))
    match <- str_match(string = projectDir, pattern = projectDirPattern)
    as.Date(match[[2L]])
}
