#' Gene to Symbol Mappings
#'
#' @name gene2symbol
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `data.frame` containing gene identifier and gene name (aka symbol)
#'   mappings.
NULL



# Methods ======================================================================
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("SummarizedExperiment"),
    function(object) {
        validObject(object)
        data <- rowData(object)
        assert_is_non_empty(data)
        data <- as.data.frame(data)
        rownames(data) <- rownames(object)
        cols <- c("geneID", "geneName")
        assert_is_subset(cols, colnames(data))
        data <- data[, cols, drop = FALSE]
        assertIsGene2symbol(data)
        data
    }
)
