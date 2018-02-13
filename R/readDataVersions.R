#' Read Data Versions
#'
#' @family File Utilities
#'
#' @importFrom basejump localOrRemoteFile
#' @importFrom readr read_csv
#'
#' @inheritParams readSampleMetadataFile
#'
#' @param file Data versions CSV file.
#'
#' @return [data.frame].
#' @export
#'
#' @examples
#' url <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "data_versions.csv")
#' readDataVersions(url)
readDataVersions <- function(
    file,
    quiet = FALSE) {
    assert_is_a_string(file)
    assert_is_a_bool(quiet)
    file <- localOrRemoteFile(file, quiet = quiet)
    if (is.null(file)) {
        return(invisible())
    }
    read_csv(file, progress = quiet)
}
