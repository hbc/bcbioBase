# nocov start
# TODO Deprecate in future release

#' Check Interesting Groups
#'
#' Prevent unwanted downstream behavior when a missing interesting group
#' is requested by the user.
#'
#' @param object Object supporting [colnames()], typically a [data.frame].
#' @param interestingGroups Interesting groups character vector.
#' @param warnOnNULL Warn the user on `NULL` argument.
#'
#' @return Valid character of defined interesting groups. Stop on failure.
#' @export
#'
#' @examples
#' demultiplexed <- paste(
#'     "http://bcbiobase.seq.cloud",
#'     "sample_metadata",
#'     "demultiplexed.xlsx",
#'     sep = "/")
#' meta <- readSampleMetadataFile(demultiplexed)
#' checkInterestingGroups(object = meta, interestingGroups = "genotype")
checkInterestingGroups <- function(
    object,
    interestingGroups,
    warnOnNULL = FALSE) {
    assert_has_colnames(object)
    assert_is_subset(interestingGroups, colnames(object))

    # Default to `sampleName` if `NULL`
    if (is.null(interestingGroups)) {
        if (isTRUE(warnOnNULL)) {
            warn(paste(
                "`interestingGroups` is NULL.",
                "Defaulting to `sampleName`."
            ))
        }
        interestingGroups <- "sampleName"
    }

    interestingGroups
}

# nocov end
