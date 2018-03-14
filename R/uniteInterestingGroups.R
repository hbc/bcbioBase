#' Unite Interesting Groups
#'
#' Create a single interesting groups column ("`interestingGroups`") used for
#' coloring in plots. When multiple interesting groups are present, unite into a
#' single column, delimited by "`:`".
#'
#' @param object Object (e.g. `data.frame`) containing interesting groups
#'   columns.
#' @param interestingGroups Character vector of interesting groups.
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' meta <- readSampleMetadataFile(
#'     "http://bcbiobase.seq.cloud/demultiplexed.csv"
#' )
#' meta <- uniteInterestingGroups(
#'     object = meta,
#'     interestingGroups = c("genotype", "sampleName")
#' )
#' meta[, "interestingGroups"]
uniteInterestingGroups <- function(object, interestingGroups) {
    assert_has_colnames(object)
    assert_is_character(interestingGroups)
    assertFormalInterestingGroups(object, interestingGroups)
    object[["interestingGroups"]] <- NULL
    object[["interestingGroups"]] <- apply(
        X = as.data.frame(object[, interestingGroups, drop = FALSE]),
        MARGIN = 1L,  # rows
        FUN = paste,
        collapse = ":"
    ) %>%
        # Ensure `interestingGroups` column is factor
        as.factor()
    object
}
