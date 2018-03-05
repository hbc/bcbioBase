#' Sample Data
#'
#' Returns the sample metadata.
#'
#' @note This is a complement to the standard [colData()] function, but improves
#'   support for accessing sample metadata for datasets where multiple items in
#'   the columns map to a single sample (e.g. cells for a single-cell RNA-seq
#'   experiment).
#'
#' @name sampleData
#'
#' @param return Return as `DataFrame` or `data.frame`.
#'
#' @return `DataFrame` or `data.frame` containing sample metadata. Note that
#'   the samples are defined in the rows, similar to [colData()].
NULL



# Constructors =================================================================
#' @importFrom basejump sanitizeColData
#' @importFrom SummarizedExperiment colData
.sampleData <- function(
    object,
    return = c("DataFrame", "data.frame"),
    ...
) {
    return <- match.arg(return)
    x <- colData(object)
    x <- sanitizeColData(x)
    as(x, return)
}



#' @importFrom SummarizedExperiment colData<-
`.sampleData<-` <- function(object, ..., value) {
    colData(object) <- value
    object
}



# Methods ======================================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("ANY"),
    .sampleData
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature("ANY"),
    `.sampleData<-`
)
