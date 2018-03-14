# TODO Improve error message if dataset is transcript level?



#' Gene to Symbol Mappings
#'
#' @name gene2symbol
#' @author Michael Steinbaugh
#'
#' @importFrom basejump gene2symbol
#'
#' @inheritParams general
#'
#' @return `data.frame` containing Ensembl gene identifier and gene name
#'   (aka symbol) mappings.
#'
#' @examples
#' \dontrun{
#' load(system.file("extdata/se.rda", package = "bcbioBase"))
#'
#' # SummarizedExperiment
#' gene2symbol(se) %>% head()
#' }
NULL



# Methods ======================================================================
#' @rdname gene2symbol
#' @importFrom SummarizedExperiment rowData
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
