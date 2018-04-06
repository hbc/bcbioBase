#' Convert Genes to Symbols
#'
#' @name convertGenesToSymbols
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `SummarizedExperiment`, with unique gene symbols as the rownames.
#'
#' @examples
#' # SummarizedExperiment ====
#' x <- convertGenesToSymbols(rse_small)
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
        assertIsGene2symbol(gene2symbol)
        symbols <- gene2symbol[, "geneName", drop = TRUE]
        symbols <- make.unique(symbols)
        rownames(object) <- symbols
        object
    }
)
