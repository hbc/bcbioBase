# nocov start

#' Defunct or Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @return No value.
NULL



# v0.2.0 =======================================================================
# `annotable()` deprecated in basejump v0.4.0
#' @rdname deprecated
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod(
    "annotable",
    signature("SummarizedExperiment"),
    function(object) {
        warn(paste(
            "'annotable' is deprecated.",
            "Use 'rowData' instead.",
            "See help(\"Deprecated\")",
            sep = "\n"
        ))
        rowData(object)
    }
)

# nocov end
