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
#' readProgramVersions("http://bcbiobase.seq.cloud/programs.txt")
readProgramVersions <- function(file) {
    assert_is_a_string(file)
    # Program versions are optional
    file <- tryCatch(
        localOrRemoteFile(file),
        error = function(e) {
            inform("Program versions are missing")
            NULL
        }
    )
    if (is.null(file)) {
        return(tibble())
    }
    # bcbio outputs `programs.txt`, but it's comma separated!
    read_csv(
        file,
        col_names = c("program", "version"),
        # `c` denotes character here
        col_types = "cc"
    )
}
