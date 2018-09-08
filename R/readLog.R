#' Read Log File
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `character`.
#' @export
#'
#' @examples
#' x <- readLog("http://bcbiobase.seq.cloud/bcbio_nextgen.log")
#' head(x)
readLog <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    read_lines(file)
}
