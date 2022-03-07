#' Import data versions
#'
#' @note The `data_versions.csv` file is only generated for special genomes
#' containing additional information (e.g. the built-in `"hg38"` build).
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
#' file <- file.path(bcbioBaseTestsURL, "data-versions.csv")
#' x <- importDataVersions(file)
#' print(x)
importDataVersions <- function(file) {
    df <- tryCatch(
        expr = {
            df <- import(
                con = file,
                format = "csv",
                engine = "base"
            )
            df <- as(df, "DataFrame")
            df
        },
        error = function(e) {
            alertWarning("Data versions are missing.")
            DataFrame()
        }
    )
    return(df)
}
