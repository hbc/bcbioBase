#' Sample YAML Metadata Utilities
#'
#' @rdname sampleYAML
#' @name sampleYAML
#' @family YAML Utilities
#'
#' @inheritParams AllGenerics
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
#' sampleYAML(yaml, "metadata")
NULL



# Constructors =================================================================
#' @importFrom basejump removeNA
#' @importFrom dplyr arrange bind_rows
#' @importFrom magrittr set_rownames
.sampleYAML <- function(yaml, keys) {
    samples <- yaml[["samples"]]
    if (!length(samples)) {
        abort("No sample information in YAML")
    }
    if (!keys[[1L]] %in% names(samples[[1L]])) {
        return(NULL)
    }
    if (length(keys) > 1L) {
        if (!keys[[2L]] %in% names(samples[[1L]][[keys[[1L]]]])) {
            return(NULL)
        }
    }
    data <- lapply(samples, function(sample) {
        nested <- sample[[keys]]
        # Set the description
        nested[["description"]] <- sample[["description"]]
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
    signature(yaml = "list",
              keys = "character"),
    .sampleYAML)
