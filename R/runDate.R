#' Run Date
#'
#' Get the run date from the project directory.
#'
#' @family Data Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @param projectDir `string`. Project directory path.
#'
#' @examples
#' runDate("2018-01-01_illumina_rnaseq")
runDate <- function(projectDir) {
    match <- str_match(
        string = basename(projectDir),
        pattern = projectDirPattern
    )
    as.Date(match[[2L]])
}
