#' Read Data Versions
#'
#' @note The `data_versions.csv` file is only generated for special genomes
#' containing additional information (e.g. the built-in "hg38" build).
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump localOrRemoteFile
#' @importFrom readr read_csv
#' @importFrom tibble tibble
#'
#' @inheritParams readSampleMetadataFile
#'
#' @param file Data versions CSV file.
#'
#' @return `tbl_df`.
#' @export
#'
#' @examples
#' readDataVersions(
#'     "http://bcbiobase.seq.cloud/data_versions.csv"
#' ) %>%
#'     glimpse()
readDataVersions <- function(file) {
    assert_is_a_string(file)
    # Warn if this file is missing
    file <- localOrRemoteFile(file, severity = "warning")
    if (is.null(file)) {
        return(tibble())
    }
    read_csv(file)
}
