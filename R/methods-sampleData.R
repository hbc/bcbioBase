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
#' @param isFactor Only include `factor` columns. Generally recommended. This
#'   helps easily discard columns in [colData()] that are `numeric` metrics.
#' @param blacklist Drop blacklisted columns defined in [metadataBlacklist].
#'   Only applies when `isFactor = TRUE`.
#'
#' @return Data describing the samples.
#'
#' @seealso [metadataBlacklist].
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
        isFactor = TRUE,
        blacklist = TRUE,
        return = c("DataFrame", "data.frame", "kable")
    ) {
        validObject(object)
        return <- match.arg(return)

        data <- colData(object)

        # Only include factor columns
        if (isTRUE(isFactor)) {
            data <- data[, vapply(data, is.factor, logical(1L)), drop = FALSE]
            # Drop blacklisted factor columns (recommended)
            if (isTRUE(blacklist)) {
                setdiff <- setdiff(colnames(data), metadataBlacklist)
                data <- data[, setdiff, drop = FALSE]
            }
        }

        # Include `interestingGroups` column
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
            if (is.character(interestingGroups)) {
                data <- uniteInterestingGroups(data, interestingGroups)
            }
        }

        # Return
        if (return == "kable") {
            data %>%
                as.data.frame() %>%
                kable(row.names = FALSE)
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
        colData(object) <- value
        object
    }
)
