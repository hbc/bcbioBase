#' Flat Files from S4 Object
#'
#' Extract the slots inside an S4 object for archival storage.
#'
#' @name flatFiles
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @return `list` containing the slots of an S4 object.
NULL



# Methods ======================================================================
#' @rdname deprecated
#' @export
setMethod(
    "flatFiles",
    signature("SummarizedExperiment"),
    function(object) {
        # Need to coerce to RSE prior to SE to keep rowRanges/rowData intact
        if (is(object, "RangedSummarizedExperiment")) {
            object <- as(object, "RangedSummarizedExperiment")
        }
        object <- as(object, "SummarizedExperiment")
        list <- as(object, "list")
        list
    }
)
