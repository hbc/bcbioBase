#' Feature Data
#'
#' Return the gene or transcript (feature) metadata.
#'
#' @note This is a complement to the standard [rowData()] function.
#'
#' @name featureData
#' @family Metadata Functions
#'
#' @inheritParams sampleData
#'
#' @param return `DataFrame`, `data.frame`, or `AsIs`.
#'
#' @return Gene/transcript (feature) metadata. Note that the features are
#' defined in the rows, similar to [rowData()].
NULL



# Constructors =================================================================
#' @importFrom SummarizedExperiment rowData
.featureData <- function(
    object,
    return = c("DataFrame", "data.frame", "AsIs"),
    ...
) {
    return <- match.arg(return)
    data <- rowData(object)
    if (return != "AsIs") {
        data <- as(data, return)
    }
    data
}



#' @importFrom basejump sanitizeRowData
#' @importFrom SummarizedExperiment rowData<-
`.featureData<-` <- function(object, ..., value) {
    value <- sanitizeRowData(value)
    rowData(object) <- value
    object
}



# Methods ======================================================================
#' @rdname featureData
#' @export
setMethod(
    "featureData",
    signature("ANY"),
    .featureData
)



#' @rdname sampleData
#' @export
setMethod(
    "featureData<-",
    signature("ANY"),
    `.featureData<-`
)
