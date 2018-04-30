#' Convert Genes to Symbols
#'
#' @name convertGenesToSymbols
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `SummarizedExperiment`, with unique gene symbols as the rownames.
#'
#' @examples
#' # SummarizedExperiment ====
#' x <- convertGenesToSymbols(rse_bcb)
#' show(x)
#' head(rownames(x))
NULL



# Methods ======================================================================
#' @rdname convertGenesToSymbols
#' @export
setMethod(
    "convertGenesToSymbols",
    signature("SummarizedExperiment"),
    function(object) {
        validObject(object)
        gene2symbol <- gene2symbol(object)
        if (is.null(gene2symbol)) {
            return(object)
        }
        symbols <- gene2symbol[, "geneName", drop = TRUE]
        # Note that ".1" will be added here for duplicate gene symbols.
        symbols <- make.unique(symbols)
        rownames(object) <- symbols
        object
    }
)
