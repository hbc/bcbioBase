#' Sample Data
#'
#' Return the sample metadata. Columns are always sanitized to factor.
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
#'
#' @examples
#' # SummarizedExperiment ====
#' sampleData(rse_dds)
#'
#' # Assignment support
#' x <- rse_dds
#' sampleData(x)[["test"]] <- seq_len(ncol(x))
#' # `test` column should be now defined
#' glimpse(sampleData(x))
NULL



# Methods ======================================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("SummarizedExperiment"),
    function(
        object,
        interestingGroups,
        return = c("DataFrame", "data.frame", "kable")
    ) {
        validObject(object)
        return <- match.arg(return)
        data <- colData(object)

        # Include `interestingGroups` column
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
            if (is.character(interestingGroups)) {
                data <- uniteInterestingGroups(data, interestingGroups)
            }
        }

        # Return
        if (return == "kable") {
            blacklist <- c("description", "fileName", "sampleID")
            data %>%
                as.data.frame() %>%
                .[, setdiff(colnames(.), blacklist), drop = FALSE] %>%
                kable(row.names = TRUE)
        } else {
            as(data, return)
        }
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "SummarizedExperiment",
        value = "DataFrame"
    ),
    function(object, value) {
        value <- as(value, "DataFrame")
        # Ensure all columns are factors
        value <- sanitizeSampleData(value)
        colData(object) <- value
        object
    }
)
