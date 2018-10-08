#' Read Program Versions
#'
#' @note bcbio doesn't save program versions when run in fast mode.
#'
#' @family Import/Export
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `tbl_df`.
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
        return(tibble())
    }
    # bcbio outputs `programs.txt`, but the file is comma separated.
    read_csv(
        file,
        col_names = c("program", "version"),
        # `c` denotes character here.
        col_types = "cc"
    )
}
