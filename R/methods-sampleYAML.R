#' Sample YAML Metadata Utilities
#'
#' @name sampleYAML
#' @family YAML Functions
#'
#' @inheritParams general
#'
#' @param yaml Project summary YAML list.
#' @param keys Nested operator keys, supplied as a character vector.
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @return [tibble].
#'
#' @examples
#' url <- paste(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "project-summary.yaml",
#'     sep = "/"
#' )
#' yaml <- basejump::readYAML(url)
#' sampleYAML(yaml, "metadata") %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom basejump removeNA
#' @importFrom dplyr arrange bind_rows
#' @importFrom magrittr set_rownames
.sampleYAML <- function(yaml, keys) {
    assert_is_non_empty(yaml)
    yaml <- yaml[["samples"]]
    assert_is_list(yaml)
    assert_is_non_empty(yaml)

    # Currently max 2 keys are supported
    assert_all_are_in_range(length(keys), lower = 1L, upper = 2L)

    # Always check that first key exists (e.g. "metadata")
    assert_is_subset(keys[[1L]], names(yaml[[1L]]))

    # When going deeper (e.g. metrics, check for secondary key presence)
    if (length(keys) == 2L) {
        assert_is_subset(
            keys[[2L]],
            names(yaml[[1L]][[keys[[1L]]]])
        )
    }

    data <- lapply(yaml, function(sample) {
        # Always get the sample description
        assert_is_character(sample[["description"]])
        description <- sample[["description"]]

        # Fix NULL batch and phenotype values found in metadata
        if (rev(keys)[[1L]] == "metadata") {
            if (is.null(sample[[keys]][["batch"]])) {
                sample[[keys]][["batch"]] <- NA
            }
            if (length(sample[[keys]][["phenotype"]])) {
                if (grepl("^$", sample[[keys]][["phenotype"]])) {
                    sample[[keys]][["phenotype"]] <- NA
                }
            }
        }

        unlist <- unlist(sample[[keys]])
        names(unlist) <- camel(names(unlist))
        # Add back the description
        unlist <- c(description = description, unlist)
        unlist
    })

    # Use this method to handle an uneven number of lengths.
    # This can happen for sample metrics values that vary in lengths per sample.
    dflist <- lapply(data, function(x) {
        as_tibble(t(x), stringsAsFactors = FALSE)
    })

    bind_rows(dflist) %>%
        removeNA() %>%
        arrange(.data[["description"]])
}



# Methods ======================================================================
#' @rdname sampleYAML
#' @export
setMethod(
    "sampleYAML",
    signature(
        yaml = "list",
        keys = "character"),
    .sampleYAML
)
