#' Interesting Groups
#'
#' @name interestingGroups
#'
#' @inheritParams general
#'
#' @return Character vector.
#'
#' @examples
#' \dontrun{
#' load(system.file("extdata/se.rda", package = "bcbioBase"))
#'
#' # SummarizedExperiment
#' interestingGroups(se)
#'
#' # Assignment support
#' interestingGroups(se) <- "sampleID"
#' interestingGroups(se)
#' }
NULL



# Methods ======================================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("SummarizedExperiment"),
    function(object) {
        validObject(object)
        metadata(object)[["interestingGroups"]]
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
            x = colData(object),
            interestingGroups = value
        )
        metadata(object)[["interestingGroups"]] <- value
        validObject(object)
        object
    }
)
