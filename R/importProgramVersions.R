#' Import program versions
#'
#' @note bcbio doesn't save program versions when run in fast mode.
#'
#' @author Michael Steinbaugh
#' @note Updated 2019-08-27.
#' @export
#'
#' @inheritParams acidroxygen::params
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "programs.txt")
#' x <- importProgramVersions(file)
#' print(x)
importProgramVersions <- function(file) {
    assert(isString(file))
    ## Program versions are optional.
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
    ## bcbio outputs `programs.txt`, but the file is comma separated.
    data <- import(
        file = file,
        format = "csv",
        colnames = c("program", "version")
    )
    as(data, "DataFrame")
}
