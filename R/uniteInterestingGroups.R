#' Unite Interesting Groups
#'
#' Create a single interesting groups column ("`interestingGroups`") used for
#' coloring in plots. When multiple interesting groups are present, unite into a
#' single column, delimited by "`:`".
#'
#' @importFrom tibble is_tibble
#' @importFrom tidyr unite
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
    if (is_tibble(object)) {
        object <- unite(
            data = object,
            col = interestingGroups,
            !!!syms(interestingGroups),
            sep = ":",
            remove = FALSE
        )
    } else {
        class <- class(object)[[1L]]
        # DataFrame is returning weird V_recycle error, so coerce
        object <- as.data.frame(object)
        object[["interestingGroups"]] <- apply(
            X = object[, interestingGroups, drop = FALSE],
            MARGIN = 1L,  # rows
            FUN = paste,
            collapse = ":"
        )
        object <- as(object, class)
    }
    # Set the `interestingGroups` column as factor
    object[["interestingGroups"]] <- as.factor(object[["interestingGroups"]])
    object
}
