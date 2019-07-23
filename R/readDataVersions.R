#' Read data versions
#'
#' @note The `data_versions.csv` file is only generated for special genomes
#' containing additional information (e.g. the built-in "`hg38`" build).
#'
#' @author Michael Steinbaugh
#' @export
#' @inheritParams basejump::params
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "data-versions.csv")
#' x <- readDataVersions(file)
#' print(x)

## Updated 2019-07-23.
readDataVersions <- function(file) {
    assert(isString(file))
    ## Data versions are optional.
    file <- tryCatch(
        localOrRemoteFile(file),
        error = function(e) {
            message("Data versions are missing.")
            NULL
        }
    )
    if (is.null(file)) {
        return(DataFrame())
    }
    data <- import(file)
    as(data, "DataFrame")
}
