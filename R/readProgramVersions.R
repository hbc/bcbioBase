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
#' url <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "programs.txt")
#' readProgramVersions(url)
readProgramVersions <- function(
    file,
    quiet = FALSE) {
    assert_is_a_string(file)
    assert_is_a_bool(quiet)
    file <- localOrRemoteFile(file, quiet = quiet)
    if (is.null(file)) {
        return(invisible())
    }
    # programs.txt, but is comma separated!
    read_csv(
        file,
        col_names = c("program", "version"),
        # `c` denotes character here
        col_types = "cc",
        progress = quiet)
}
