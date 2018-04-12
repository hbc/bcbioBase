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
#' readLog("http://bcbiobase.seq.cloud/bcbio-nextgen.log") %>%
#'     head()
readLog <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    read_lines(file)
}
