#' Unite Interesting Groups
#'
#' Create a single interesting groups column (`interestingGroups`) used for
#' coloring in plots. When multiple interesting groups are present, unite into a
#' single column, delimited by `:`.
#'
#' @importFrom tidyr unite
#'
#' @param object Object (e.g. [data.frame]) containing interesting groups
#'   columns.
#' @param interestingGroups Character vector of interesting groups.
#'
#' @return [data.frame].
#' @export
#'
#' @examples
#' demultiplexed <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "sample_metadata",
#'     "demultiplexed.xlsx")
#' meta <- readSampleMetadataFile(demultiplexed)
#' meta <- uniteInterestingGroups(
#'     object = meta,
#'     interestingGroups = c("genotype", "sampleName")
#' )
#' pull(meta, "interestingGroups")
uniteInterestingGroups <- function(object, interestingGroups) {
    assert_has_colnames(object)
    assert_is_character(interestingGroups)
    assert_formal_interesting_groups(object, interestingGroups)

    object[["interestingGroups"]] <- NULL

    if (length(interestingGroups) > 1L) {
        # For multiple groups, unite to a colon-separated string
        object <- unite(
            data = object,
            col = interestingGroups,
            !!!syms(interestingGroups),
            sep = ":",
            remove = FALSE)
    } else {
        object[["interestingGroups"]] <- object[[interestingGroups]]
    }

    # Set the `interestingGroups` column as factor
    object[["interestingGroups"]] <- as.factor(object[["interestingGroups"]])
    object
}
