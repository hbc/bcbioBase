#' Read Program Versions
#'
#' @note bcbio doesn't save program versions when run in fast mode.
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `tbl_df`.
#' @export
#'
#' @examples
#' x <- readProgramVersions("http://bcbiobase.seq.cloud/programs.txt")
#' print(x)
readProgramVersions <- function(file) {
    assert_is_a_string(file)
    # Program versions are optional
    file <- tryCatch(
        localOrRemoteFile(file),
        error = function(e) {
            message("Program versions are missing")
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
