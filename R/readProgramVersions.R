#' Read Program Versions
#'
#' @note bcbio doesn't save program versions when run in fast mode.
#'
#' @name readProgramVersions
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump localOrRemoteFile
#' @importFrom readr read_csv
#'
#' @inheritParams readSampleMetadataFile
#'
#' @param file Program versions TXT file.
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' readProgramVersions(
#'     "http://bcbiobase.seq.cloud/programs.txt"
#' )
readProgramVersions <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file, severity = "warning")
    if (is.null(file)) {
        return(NULL)
    }
    # bcbio outputs programs.txt, but is comma separated!
    read_csv(
        file,
        col_names = c("program", "version"),
        # `c` denotes character here
        col_types = "cc"
    )
}
