#' Read Program Versions
#'
#' @note bcbio doesn't save program versions when run in fast mode.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams basejump::params
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "programs.txt")
#' x <- readProgramVersions(file)
#' print(x)
readProgramVersions <- function(file) {
    assert_is_a_string(file)
    # Program versions are optional
    file <- tryCatch(
        localOrRemoteFile(file),
        error = function(e) {
            message("Program versions are missing.")
            NULL
        }
    )
    if (is.null(file)) {
        return(DataFrame())
    }
    # bcbio outputs `programs.txt`, but the file is comma separated.
    data <- read_csv(
        file,
        col_names = c("program", "version"),
        col_types = "cc"  # character
    )
    as(data, "DataFrame")
}
