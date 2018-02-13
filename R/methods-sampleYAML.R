#' Sample YAML Metadata Utilities
#'
#' @rdname sampleYAML
#' @name sampleYAML
#' @family YAML Utilities
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
#' @return [data.frame].
#'
#' @examples
#' url <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "project-summary.yaml")
#' yaml <- basejump::readYAML(url)
#' sampleYAML(yaml, "metadata") %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom basejump removeNA
#' @importFrom dplyr arrange bind_rows
#' @importFrom magrittr set_rownames
.sampleYAML <- function(yaml, keys) {
    yaml <- yaml[["samples"]]
    assert_is_list(yaml)
    assert_is_non_empty(yaml)
    assert_is_character(keys)
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
        nested <- sample[[keys]]

        # Always set the description, if present
        if (is.character(sample[["description"]])) {
            nested[["description"]] <- sample[["description"]]
        }

        # Fix NULL batch and phenotype values found in metadata
        if (rev(keys)[[1L]] == "metadata") {
            if (is.null(nested[["batch"]])) {
                nested[["batch"]] <- NA
            }
            if (length(nested[["phenotype"]])) {
                if (grepl("^$", nested[["phenotype"]])) {
                    nested[["phenotype"]] <- NA
                }
            }
        }

        unlist <- unlist(nested)
        names(unlist) <- camel(names(unlist), strict = FALSE)
        unlist
    })

    dflist <- lapply(data, function(x) {
        as.data.frame(t(x), stringsAsFactors = FALSE)
    })

    bind_rows(dflist) %>%
        removeNA() %>%
        arrange(.data[["description"]]) %>%
        set_rownames(.[["description"]])
}



# Methods ======================================================================
#' @rdname sampleYAML
#' @export
setMethod(
    "sampleYAML",
    signature(
        yaml = "list",
        keys = "character"),
    .sampleYAML)
