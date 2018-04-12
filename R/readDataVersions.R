#' Read Data Versions
#'
#' @note The `data_versions.csv` file is only generated for special genomes
#' containing additional information (e.g. the built-in "hg38" build).
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
#' readDataVersions("http://bcbiobase.seq.cloud/data_versions.csv") %>%
#'     glimpse()
readDataVersions <- function(file) {
    assert_is_a_string(file)
    # Data versions are optional
    file <- tryCatch(
        localOrRemoteFile(file),
        error = function(e) {
            message("Data versions are missing")
            NULL
        }
    )
    if (is.null(file)) {
        return(tibble())
    }
    read_csv(file)
}
