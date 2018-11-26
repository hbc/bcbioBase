# Alternatively, can parse YAML data, but this is faster and simpler.
#' Run Date
#'
#' Get the run date from the project directory.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams params
#'
#' @examples
#' runDate("2018-01-01_illumina_rnaseq")
runDate <- function(projectDir) {
    assert_is_a_string(projectDir)
    projectDir <- basename(projectDir)
    assert_that(grepl(projectDirPattern, projectDir))
    match <- str_match(
        string = projectDir,
        pattern = projectDirPattern
    )
    as.Date(match[[2L]])
}
