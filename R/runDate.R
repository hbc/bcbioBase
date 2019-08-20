#' Run date
#'
#' Get the run date from the project directory.
#'
#' Alternatively, can parse YAML data, but this approach is faster and simpler.
#'
#' @author Michael Steinbaugh
#' @note Updated 2019-08-20.
#' @export
#'
#' @inheritParams acidroxygen::params
#'
#' @return `Date`.
#'
#' @examples
#' runDate(projectDir = "2018-01-01_illumina_rnaseq")
runDate <- function(projectDir) {
    projectDir <- basename(projectDir)
    assert(
        isString(projectDir),
        allAreMatchingRegex(projectDir, pattern = projectDirPattern)
    )
    match <- str_match(string = projectDir, pattern = projectDirPattern)
    as.Date(match[[2L]])
}
