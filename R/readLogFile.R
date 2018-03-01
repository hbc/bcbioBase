#' Read Log File
#'
#' @family File Utilities
#'
#' @importFrom basejump localOrRemoteFile
#' @importFrom readr read_lines
#'
#' @inheritParams readSampleMetadataFile
#'
#' @param file Log file.
#'
#' @return Character vector.
#' @export
#'
#' @examples
#' url <- paste(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "bcbio-nextgen.log",
#'     sep = "/")
#' readLogFile(url) %>% head()
readLogFile <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    assert_all_are_existing_files(file)
    read_lines(file)
}
