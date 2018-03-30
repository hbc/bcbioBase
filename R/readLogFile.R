#' Read Log File
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return Character vector.
#' @export
#'
#' @examples
#' readLogFile("http://bcbiobase.seq.cloud/bcbio-nextgen.log") %>%
#'     head()
readLogFile <- function(file) {
    assert_is_a_string(file)
    # Log files are always required
    file <- localOrRemoteFile(file)
    assert_all_are_existing_files(file)
    read_lines(file)
}
