#' Import data versions
#'
#' @note The `data_versions.csv` file is only generated for special genomes
#' containing additional information (e.g. the built-in `"hg38"` build).
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
#' file <- file.path(bcbioBaseTestsURL, "data-versions.csv")
#' x <- importDataVersions(file)
#' print(x)
importDataVersions <- function(file) {
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
