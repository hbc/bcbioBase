#' Interesting Groups
#'
#' @name interestingGroups
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `character`.
#'
#' @examples
#' # SummarizedExperiment ====
#' interestingGroups(rse_dds)
#' colnames(colData(rse_dds))
#' interestingGroups(rse_dds) <- colnames(colData(rse_dds))[[1L]]
#' interestingGroups(rse_dds)
NULL



# Methods ======================================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("SummarizedExperiment"),
    function(object) {
        validObject(object)
        value <- metadata(object)[["interestingGroups"]]
        assertFormalInterestingGroups(
            x = sampleData(object, clean = FALSE, interestingGroups = NULL),
            interestingGroups = value
        )
        value
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
            x = sampleData(object, clean = FALSE, interestingGroups = NULL),
            interestingGroups = value
        )
        metadata(object)[["interestingGroups"]] <- value
        validObject(object)
        object
    }
)
