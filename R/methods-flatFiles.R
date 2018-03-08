#' Get Flat Files from S4 Object
#'
#' Prepare an unstructured list of information saved in the S4 object,
#' for improved archival data storage.
#'
#' @name flatFiles
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `list`.
#'
#' @examples
#' \dontrun{
#' load(system.file("extdata/se.rda", package = "bcbioBase"))
#'
#' # SummarizedExperiment
#' flatFiles(se) %>% names()
#' }
NULL



# Methods ======================================================================
#' @rdname flatFiles
#' @export
setMethod(
    "flatFiles",
    signature("SummarizedExperiment"),
    function(object) {
        list(
            assays = assays(object),
            rowData = rowData(object),
            colData = colData(object),
            metadata = metadata(object)
        )
    }
)
