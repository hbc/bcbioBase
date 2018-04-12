#' Unite Interesting Groups
#'
#' Create a single interesting groups column ("`interestingGroups`") used for
#' coloring in plots. When multiple interesting groups are present, unite into a
#' single column, delimited by "`:`".
#'
#' @name uniteInterestingGroups
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param object Object containing interesting groups in multiple columns.
#'
#' @return Object of same class, containing a united `interestingGroups` column.
#' @export
#'
#' @examples
#' x <- uniteInterestingGroups(
#'     object = colData(rse_bcb),
#'     interestingGroups = c("treatment", "day")
#' )
#' x[, "interestingGroups"]
NULL



#' Methods =====================================================================
#' @rdname uniteInterestingGroups
#' @export
setMethod(
    "uniteInterestingGroups",
    signature("data.frame"),
    function(object, interestingGroups) {
        assert_has_colnames(object)
        assert_is_character(interestingGroups)
        assertFormalInterestingGroups(object, interestingGroups)
        class <- class(object)[[1L]]
        object <- as.data.frame(object)
        intgroup <- apply(
            X = object[, interestingGroups, drop = FALSE],
            MARGIN = 1L,
            FUN = paste,
            collapse = ":"
        )
        object[["interestingGroups"]] <- as.factor(intgroup)
        as(object, class)
    }
)



#' @rdname uniteInterestingGroups
#' @export
setMethod(
    "uniteInterestingGroups",
    signature("DataFrame"),
    getMethod("uniteInterestingGroups", "data.frame")
)



#' @rdname uniteInterestingGroups
#' @export
setMethod(
    "uniteInterestingGroups",
    signature("tbl_df"),
    function(object, interestingGroups) {
        assert_is_any_of(object, "tbl_df")
        assert_has_colnames(object)
        assert_is_character(interestingGroups)
        assertFormalInterestingGroups(object, interestingGroups)
        object[["interestingGroups"]] <- NULL
        object <- unite(
            data = object,
            col = interestingGroups,
            !!interestingGroups,
            sep = ":",
            remove = FALSE
        )
        object[["interestingGroups"]] <-
            as.factor(object[["interestingGroups"]])
        object
    }
)
