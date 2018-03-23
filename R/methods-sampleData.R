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
#' @return Data describing the samples.
NULL



# Constructors =================================================================
.sampleData <- function(
    object,
    return = c("data.frame", "DataFrame", "kable")
) {
    return <- match.arg(return)
    data <- colData(object)
    # Ensure all columns are factors
    data <- sanitizeSampleData(data)
    if (return == "kable") {
        blacklist <- c("description", "fileName", "sampleID")
        data %>%
            as.data.frame() %>%
            .[, setdiff(colnames(.), blacklist), drop = FALSE] %>%
            # Ensure `sampleName` is first
            .[, unique(c("sampleName", colnames(.))), drop = FALSE] %>%
            # Arrange by `sampleName`
            .[order(.[["sampleName"]]), , drop = FALSE] %>%
            kable(row.names = FALSE)
    } else {
        as(data, return)
    }
}



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
