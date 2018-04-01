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
        list <- lapply(slotNames(object), function(slot) {
            if (.hasSlot(object, slot)) {
                slot(object, slot)
            } else {
                NULL
            }
        })
        names(list) <- slotNames(object)
        list
    }
)
