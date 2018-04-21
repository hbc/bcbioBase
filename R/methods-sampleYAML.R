#' Sample YAML Metadata
#'
#' @name sampleYAML
#' @family YAML Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param yaml Project summary YAML `list`.
#' @param keys Nested operator keys, supplied as a character vector.
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @return `data.frame`.
#'
#' @examples
#' yaml <- basejump::readYAML(
#'     "http://bcbiobase.seq.cloud/project-summary.yaml"
#' )
#'
#' sampleYAML(yaml, keys = "metadata") %>% glimpse()
#' sampleYAMLMetadata(yaml) %>% glimpse()
#' sampleYAMLMetrics(yaml) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname sampleYAML
#' @export
setMethod(
    "sampleYAML",
    signature(
        yaml = "list",
        keys = "character"
    ),
    function(yaml, keys) {
        assert_is_non_empty(yaml)
        assert_is_subset("samples", names(yaml))
        yaml <- yaml[["samples"]]
        assert_is_list(yaml)
        assert_is_non_empty(yaml)

        # Currently max 2 keys are supported (e.g. summary, metrics)
        assert_all_are_in_range(length(keys), lower = 1L, upper = 2L)

        # Check that keys are present and early return on failure
        if (!keys[[1L]] %in% names(yaml[[1L]])) {
            warning(paste(
                deparse(keys[[1L]]),
                "missing in sample YAML")
            )
            return(NULL)
        } else if (
            length(keys) == 2L &&
            !keys[[2L]] %in% names(yaml[[1L]][[keys[[1L]]]])
        ) {
            warning(paste(
                deparse(keys[[2L]]),
                "missing in sample YAML"
            ))
            return(NULL)
        }

        list <- lapply(yaml, function(x) {
            # Always get the sample description
            description <- x[["description"]]
            assert_is_character(description)
            assert_is_non_empty(description)

            x <- x[[keys]]
            assert_is_non_empty(x)
            # Always sanitize names to camel case
            x <- camel(x)
            # Add description column
            x[["description"]] <- description

            # Coerce nested elements to string, if necessary.
            # Consider adding a warning here about this behavior for metadata.
            x <- lapply(x, function(x) {
                if (length(x) > 1L) {
                    toString(x)
                } else {
                    x
                }
            })

            # Remove any `NULL` items
            Filter(Negate(is.null), x)
        })

        # Use this method to coerce a list with uneven lengths
        ldply(list, data.frame, stringsAsFactors = FALSE) %>%
            fixNA() %>%
            removeNA() %>%
            .[, sort(colnames(.))] %>%
            arrange(!!sym("description")) %>%
            set_rownames(.[["description"]])
    }
)



#' @rdname sampleYAML
#' @export
setMethod(
    "sampleYAMLMetadata",
    signature("list"),
    function(yaml) {
        sampleYAML(yaml = yaml, keys = "metadata") %>%
            prepareSampleData()
    }
)



#' @rdname sampleYAML
#' @export
setMethod(
    "sampleYAMLMetrics",
    signature("list"),
    function(yaml) {
        data <- sampleYAML(
            yaml = yaml,
            keys = c("summary", "metrics")
        )

        # Early return on empty metrics
        if (!length(data)) {
            warning("Fast mode detected: No sample metrics were calculated")
            return(NULL)
        }

        # Fix numerics set as characters
        numericAsCharacter <- function(x) {
            any(grepl(x = x, pattern = "^[0-9\\.]+$"))
        }

        data %>%
            rownames_to_column() %>%
            mutate_if(is.factor, as.character) %>%
            mutate_if(numericAsCharacter, as.numeric) %>%
            mutate_if(is.character, as.factor) %>%
            column_to_rownames() %>%
            # Drop any metadata columns. Note we're also dropping the duplicate
            # `name` column present in the metrics YAML.
            .[,
              sort(setdiff(
                  x = colnames(.),
                  y = c(metadataPriorityCols, "name")
              )),
              drop = FALSE]
    }
)
