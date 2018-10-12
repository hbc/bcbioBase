#' Read Log File
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `character`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "bcbio_nextgen.log")
#' x <- readLog(file)
#' head(x)
readLog <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    read_lines(file)
}
