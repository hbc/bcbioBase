#' Import program versions
#'
#' @note bcbio doesn't save program versions when run in fast mode.
#'
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "programs.txt")
#' x <- importProgramVersions(file)
#' print(x)
importProgramVersions <- function(file) {
    tryCatch(
        expr = {
            df <- import(
                con = file,
                format = "csv",
                colnames = c("program", "version"),
                engine = "base"
            )
            df <- as(df, "DataFrame")
            df
        },
        error = function(e) {
            alertWarning("Program versions are missing.")
            DataFrame()
        }
    )
}
