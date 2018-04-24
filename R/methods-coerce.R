#' Methods for Coercing an Object to a Class
#'
#' @name coerce
#' @aliases as
#' @family S4 Class Definition
#' @author Michael Steinbaugh
#'
#' @return Object of new class.
#'
#' @seealso
#' - [methods::as()].
#' - [methods::coerce()].
#'
#' @examples
#' # SummarizedExperiment to list ====
#' x <- as(rse_dds, "list")
#' names(x)
NULL



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-SummarizedExperiment-list
setAs(
    from = "SummarizedExperiment",
    to = "list",
    function(from) {
        flatFiles(from)
    }
)
