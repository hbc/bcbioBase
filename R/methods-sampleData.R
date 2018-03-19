#' Sample Data
#'
#' Return the sample metadata.
#'
#' @note This is a complement to the standard [colData()] function, but improves
#'   support for accessing sample metadata for datasets where multiple items in
#'   the columns map to a single sample (e.g. cells for a single-cell RNA-seq
#'   experiment).
#'
#' @name sampleData
#' @family Metadata Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param return `data.frame` (default), `DataFrame`, or unmodified (`AsIs`).
#'
#' @return Sample metadata. Note that the samples are defined in the rows,
#'   similar to [colData()].
NULL



# Constructors =================================================================
#' @importFrom basejump sanitizeSampleData
#' @importFrom SummarizedExperiment colData
.sampleData <- function(
    object,
    return = c("data.frame", "DataFrame", "AsIs"),
    ...
) {
    return <- match.arg(return)
    data <- colData(object, ...)
    # Ensure all columns are factors
    data <- sanitizeSampleData(data)
    if (return != "AsIs") {
        data <- as(data, return)
    }
    data
}



#' @importFrom basejump sanitizeSampleData
#' @importFrom SummarizedExperiment colData<-
`.sampleData<-` <- function(object, ..., value) {
    # Ensure all columns are factors
    value <- sanitizeSampleData(value)
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
