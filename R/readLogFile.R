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
#' url <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "bcbio-nextgen.log")
#' readLogFile(url) %>% head()
readLogFile <- function(
    file,
    quiet = FALSE) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file, quiet = quiet)
    if (is.null(file)) {
        return(invisible())
    }
    read_lines(file)
}
