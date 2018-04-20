#' Select Samples
#'
#' @name selectSamples
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `SummarizedExperiment`.
#'
#' @examples
#' # SummarizedExperiment ====
#' x <- selectSamples(rse_bcb, day = 7L)
#' show(x)
#' colnames(x)
NULL



# Methods ======================================================================
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("SummarizedExperiment"),
    function(object, ...) {
        validObject(object)
        args <- list(...)
        invisible(lapply(args, assert_is_atomic))

        # Match the arguments against the sample metadata
        colData <- colData(object)
        assert_is_subset(names(args), colnames(colData))

        # Obtain the sample identifiers
        list <- mapply(
            col = names(args),
            arg = args,
            MoreArgs = list(data = colData),
            FUN = function(col, arg, data) {
                rownames(data[data[[col]] %in% arg, , drop = FALSE])
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )
        samples <- Reduce(f = intersect, x = list) %>%
            as.character() %>%
            sort()
        assert_is_non_empty(samples)

        object[, samples]
    }
)
