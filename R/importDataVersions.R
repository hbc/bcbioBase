## FIXME Take out code using localOrRemoteFile.



#' Import data versions
#'
#' @note The `data_versions.csv` file is only generated for special genomes
#' containing additional information (e.g. the built-in `"hg38"` build).
#'
#' @author Michael Steinbaugh
#' @note Updated 2020-01-17.
#' @export
#'
#' @inheritParams AcidRoxygen::params
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
            alertWarning("Data versions are missing.")
            NULL
        }
    )
    if (is.null(file)) {
        return(DataFrame())
    }
    data <- import(file)
    as(data, "DataFrame")
}
