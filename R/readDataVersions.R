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
#' url <- paste(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "data_versions.csv",
#'     sep = "/")
#' readDataVersions(url) %>% glimpse()
readDataVersions <- function(file) {
    assert_is_a_string(file)
    # bcbio doesn't save this in fast mode
    file <- localOrRemoteFile(file, severity = "warning")
    if (is.null(file)) {
        return(NULL)
    }
    read_csv(file)
}
