#' Interesting Groups
#'
#' @name interestingGroups
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return Character vector.
#'
#' @examples
#' # SummarizedExperiment ====
#' interestingGroups(rse_small)
NULL



# Methods ======================================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("SummarizedExperiment"),
    function(object) {
        validObject(object)
        x <- metadata(object)[["interestingGroups"]]
        ifelse(is.character(x), x, "sampleName")
    }
)



#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "SummarizedExperiment",
        value = "character"
    ),
    function(object, value) {
        assertFormalInterestingGroups(
            x = sampleData(object),
            interestingGroups = value
        )
        metadata(object)[["interestingGroups"]] <- value
        validObject(object)
        object
    }
)
