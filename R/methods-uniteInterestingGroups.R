#' Unite Interesting Groups
#'
#' Create a single interesting groups column ("`interestingGroups`") used for
#' coloring in plots. When multiple interesting groups are present, unite into a
#' single column, delimited by "`:`".
#'
#' @name uniteInterestingGroups
#' @author Michael Steinbaugh
#'
#' @param object Object (e.g. `data.frame`) containing interesting groups
#'   columns.
#' @param interestingGroups Character vector of interesting groups.
#'
#' @return Object of same class, containing `interestingGroups` column.
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
NULL



# Constructors =================================================================
.uniteInterestingGroups.base <- function(object, interestingGroups) {
    assert_has_colnames(object)
    assert_is_character(interestingGroups)
    assertFormalInterestingGroups(object, interestingGroups)
    class <- class(object)[[1L]]
        # DataFrame is returning weird V_recycle error, so coerce
        object <- as.data.frame(object)
        object[["interestingGroups"]] <- apply(
            X = object[, interestingGroups, drop = FALSE],
            MARGIN = 1L,  # rows
            FUN = paste,
            collapse = ":"
        ) %>%
            as.factor()
    as(object, class)
}



#' @importFrom tibble is_tibble
#' @importFrom tidyr unite
.uniteInterestingGroups.tidy <- function(object, interestingGroups) {
    assert_is_any_of(object, "tbl_df")
    assert_has_colnames(object)
    assert_is_character(interestingGroups)
    assertFormalInterestingGroups(object, interestingGroups)
    object[["interestingGroups"]] <- NULL
    object <- unite(
        data = object,
        col = interestingGroups,
        !!!syms(interestingGroups),
        sep = ":",
        remove = FALSE
    )
    object[["interestingGroups"]] <- as.factor(object[["interestingGroups"]])
    object
}


#' Methods =====================================================================
#' @rdname uniteInterestingGroups
#' @export
setMethod(
    "uniteInterestingGroups",
    signature("DataFrame"),
    .uniteInterestingGroups.base
)



#' @rdname uniteInterestingGroups
#' @export
setMethod(
    "uniteInterestingGroups",
    signature("data.frame"),
    .uniteInterestingGroups.base
)



#' @rdname uniteInterestingGroups
#' @export
setMethod(
    "uniteInterestingGroups",
    signature("tbl_df"),
    .uniteInterestingGroups.tidy
)
