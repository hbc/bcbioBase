#' Read Data Versions
#'
#' @note bcbio doesn't save data versions when run in fast mode.
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump localOrRemoteFile
#' @importFrom readr read_csv
#'
#' @inheritParams readSampleMetadataFile
#'
#' @param file Data versions CSV file.
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' readDataVersions(
#'     "http://bcbiobase.seq.cloud/data_versions.csv"
#' ) %>%
#'     glimpse()
readDataVersions <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file, severity = "warning")
    if (is.null(file)) {
        return(NULL)
    }
    read_csv(file)
}
