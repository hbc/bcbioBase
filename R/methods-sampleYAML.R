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
#' @return `tbl_df`.
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
            arrange(!!sym("description"))
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
        # Early return on NULL metrics (fast mode)
        fastMode <- "Fast mode detected: No sample metrics were calculated"

        if (is.null(yaml[["samples"]][[1L]][["summary"]][["metrics"]])) {
            warning(fastMode)
            return(NULL)
        }

        data <- sampleYAML(
            yaml = yaml,
            keys = c("summary", "metrics")
        )
        assert_is_tbl(data)

        if (identical(colnames(data), "description")) {
            warning(fastMode)
            return(NULL)
        }

        # Fix numerics set as characters
        numericAsCharacter <- function(x) {
            any(grepl(x = x, pattern = "^[0-9\\.]+$"))
        }

        data %>%
            mutate_if(is.factor, as.character) %>%
            mutate_if(numericAsCharacter, as.numeric) %>%
            mutate_if(is.character, as.factor) %>%
            as.data.frame %>%
            set_rownames(makeNames(.[["description"]], unique = TRUE)) %>%
            # Drop any metadata columns
            .[,
              sort(setdiff(
                  x = colnames(.),
                  y = c(metadataPriorityCols, "name")
              )),
              drop = FALSE]
    }
)
