#' Interesting Groups Formal Assert Check
#'
#' Prevent unwanted downstream behavior when a missing interesting group
#' is requested by the user.
#'
#' @inherit assert
#'
#' @param object Object supporting [colnames()], typically a [data.frame].
#' @param interestingGroups Interesting groups character vector.
#'
#' @return Valid character of defined interesting groups. Stop on failure.
#' @export
#'
#' @examples
#' demultiplexed <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "sample_metadata",
#'     "demultiplexed.xlsx")
#' meta <- readSampleMetadataFile(demultiplexed)
#' assert_formal_interesting_groups(meta, "genotype")
assert_formal_interesting_groups <- function(
    object,
    interestingGroups,
    severity = "stop") {
    assert_has_colnames(object, severity = severity)
    assert_is_subset(
        x = interestingGroups,
        y = colnames(object),
        severity = severity
    )
}
