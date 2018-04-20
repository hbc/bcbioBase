#' Counts
#'
#' @name counts
#' @family Data Functions
#'
#' @importFrom BiocGenerics counts
#'
#' @inheritParams general
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
