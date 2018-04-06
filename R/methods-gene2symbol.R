#' Gene to Symbol Mappings
#'
#' @name gene2symbol
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `data.frame` containing gene identifier and gene name (aka symbol)
#'   mappings.
#'
#' @examples
#' # SummarizedExperiment ====
#' gene2symbol(rse_small)
NULL



# Methods ======================================================================
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("SummarizedExperiment"),
    function(object) {
        validObject(object)
        rowData <- rowData(object)
        cols <- c("geneID", "geneName")
        if (!all(cols %in% colnames(rowData))) {
            warn("`rowData(object)` does not contain gene-to-symbol mappings")
            return(NULL)
        }
        assert_is_non_empty(data)
        data <- as.data.frame(data)
        rownames(data) <- rownames(object)
        assert_is_subset(cols, colnames(data))
        data <- data[, cols, drop = FALSE]
        assertIsGene2symbol(data)
        data
    }
)
