#' Sample Names
#'
#' Requires that `sampleName` column is defined in [sampleData()].
#'
#' @name sampleNames
#' @family Metadata Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `character` vector of the human readable sample names.
#'
#' @examples
#' # SummarizedExperiment ====
#' sampleNames(rse_bcb)
NULL



# Methods ======================================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleNames",
    signature("SummarizedExperiment"),
    function(object) {
        validObject(object)
        assert_is_subset("sampleName", colnames(sampleData(object)))
        sampleData(object) %>%
            .[, "sampleName", drop = TRUE] %>%
            as.character() %>%
            sort()
    }
)
