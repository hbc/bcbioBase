#' Read Program Versions
#'
#' @rdname readProgramVersions
#' @name readProgramVersions
#' @family bcbio Utilities
#' @keywords internal
#'
#' @importFrom basejump localOrRemoteFile
#' @importFrom readr read_csv
#'
#' @inheritParams readSampleMetadataFile
#'
#' @param file Program versions TXT file.
#'
#' @return [data.frame].
#' @export
#'
#' @examples
#' url <- paste(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "programs.txt",
#'     sep = "/")
#' readProgramVersions(url)
readProgramVersions <- function(file) {
    assert_is_a_string(file)
    # bcbio doesn't save this in fast mode
    file <- localOrRemoteFile(file, severity = "warning")
    if (is.null(file)) {
        return(NULL)
    }
    # bcbio outputs programs.txt, but is comma separated!
    read_csv(
        file,
        col_names = c("program", "version"),
        # `c` denotes character here
        col_types = "cc")
}
