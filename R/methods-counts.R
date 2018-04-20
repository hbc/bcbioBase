#' Counts
#'
#' @name counts
#' @family Data Functions
#'
#' @importFrom BiocGenerics counts
#'
#' @examples
#' # SummarizedExperiment ====
#' x <- counts(rse_bcb)
#' glimpse(x)
NULL



# Methods ======================================================================
#' @rdname counts
#' @export
setMethod(
    "counts",
    "SummarizedExperiment",
    function(object) {
        stopifnot("counts" %in% names(assays(object)))
        assays(object)[["counts"]]
    }
)
