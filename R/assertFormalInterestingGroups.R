#' Interesting Groups Formal Assert Check
#'
#' Prevent unwanted downstream behavior when a missing interesting group
#' is requested by the user.
#'
#' @family Assert Check Functions
#' @author Michael Steinbaugh
#'
#' @inherit assert
#'
#' @param x Object supporting [colnames()], typically a `data.frame`.
#' @param interestingGroups Interesting groups character vector.
#'
#' @return Silent, stop on error.
#' @export
#'
#' @examples
#' data <- readSampleData("http://bcbiobase.seq.cloud/demultiplexed.csv")
#' assertFormalInterestingGroups(data, interestingGroups = "genotype")
assertFormalInterestingGroups <- function(
    x,
    interestingGroups,
    severity = "stop"
) {
    assert_has_colnames(x, severity = severity)
    assert_is_subset(
        x = interestingGroups,
        y = colnames(x),
        severity = severity
    )
}
