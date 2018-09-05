#' Read YAML Sample Data and Metrics
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @name readYAML
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param file `string`. File path to bcbio `project-summary.yaml`.
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- "http://bcbiobase.seq.cloud/project-summary.yaml"
#' readYAMLSampleData(file) %>% glimpse()
#' readYAMLSampleMetrics(file) %>% glimpse()
NULL



.readYAMLSample <- function(file, keys) {
    assert_is_character(keys)
    # Currently max 2 keys are supported (e.g. summary, metrics).
    assert_all_are_in_range(length(keys), lower = 1L, upper = 2L)

    yaml <- readYAML(file)
    assert_are_identical(
        x = names(yaml),
        y = c("date", "upload", "bcbio_system", "samples")
    )

    # Focus on the sample YAML data.
    yaml <- yaml[["samples"]]
    assert_is_list(yaml)
    assert_is_non_empty(yaml)

    # `summary` is only returned for RNA-seq pipeline, not single cell.
    assert_is_subset(
        x = c(
            "description",
            "dirs",
            "genome_build",
            "genome_resources",
            "metadata",
            "sam_ref"
        ),
        y = names(yaml[[1L]])
    )

    # Check that nested keys are present and early return on failure.
    # Return `NULL` here instead of stopping, so we can handle bcbio RNA-seq
    # fast mode runs.
    if (
        length(keys) == 2L &&
        !keys[[2L]] %in% names(yaml[[1L]][[keys[[1L]]]])
    ) {
        return(NULL)
    }

    # Top-level sample metadata in the YAML is relatively easy to parse.
    # Just select for atomic values.
    top <- sapply(
        X = yaml,
        FUN = function(item) {
            x <- sapply(
                X = item,
                FUN = function(item) {
                    if (is.atomic(item)) {
                        item
                    } else {
                        NULL
                    }
                },
                simplify = TRUE,
                USE.NAMES = TRUE
            )
            Filter(Negate(is.null), x)
        },
        simplify = TRUE,
        USE.NAMES = FALSE
    )
    top_df <- as.data.frame(top)
    dim(top_df)
    top_DF <- as(top, "DataFrame")
    dim(top_DF)
    # Coerce to data.frame.
    # TODO Make this a utility function in basejump.
    top <- ldply(
        .data = top,
        .fun = data.frame,
        stringsAsFactors = FALSE
    )
    assert_is_data.frame(top)
    assert_is_non_empty(top)
    top <- as(top, "DataFrame")

    # Handle the nested metadata, defined by the keys.
    # This step is a little tricker but should work consistently.
    nested <- lapply(yaml, function(item) {
        item <- item[[keys]]
        assert_is_non_empty(item)
        # Sanitize names into camel case here, otherwise they'll get modified
        # during the `ldply()` call that coerces `list` to `data.frame`.
        item <- camel(item)
        item <- lapply(item, function(item) {
            if (length(item) > 1L) {
                # Detect and coerce nested metadata back to a string, if
                # necessary. bcbio allows nesting with a semicolon delimiter.
                # Warn the user here about discouraging this with R data.
                # http://bit.ly/2Je1xgO
                warning("Nested sample metadata detected")
                paste(item, collapse = "; ")
            } else {
                item
            }
        })
        # Remove any `NULL` items.
        Filter(Negate(is.null), item)
    })
    # Coerce to data frame.
    # Note that we're using `plyr::ldply()` here because it can coerce a list
    # with uneven lengths. However, it will sanitize column names, so make sure
    # we convert to camel case before running this step.
    nested <- ldply(
        .data = nested,
        .fun = data.frame,
        stringsAsFactors = FALSE
    )
    assert_is_data.frame(nested)
    assert_is_non_empty(nested)
    nested <- as(nested, "DataFrame")

    # Bind the top and nested data frames, coerce to tibble, and return.
    assert_are_disjoint_sets(
        x = colnames(top),
        y = colnames(nested)
    )
    cbind(top, nested) %>%
        as("tbl_df") %>%
        sanitizeNA() %>%
        removeNA() %>%
        camel() %>%
        arrange(!!sym("description")) %>%
        # Order the columns alphabetically.
        .[, sort(colnames(.)), drop = FALSE] %>%
        # Set the rownames.
        mutate(rowname = makeNames(!!sym("description"), unique = TRUE))
}



#' @rdname readYAML
#' @export
readYAMLSampleData <- function(file) {
    message("Reading sample metadata from YAML")
    file %>%
        .readYAMLSample(keys = "metadata") %>%
        .returnSampleData()

}



#' @rdname readYAML
#' @export
readYAMLSampleMetrics <- function(file) {
    message("Reading sample metrics from YAML")
    data <- .readYAMLSample(file, keys = c("summary", "metrics"))

    # Early return on empty metrics (e.g. fast mode).
    if (!length(data)) {
        return(NULL)  # nocov
    }

    # Fix numerics set as characters.
    numericAsCharacter <- function(x) {
        any(grepl(x = x, pattern = "^[0-9\\.]+$"))
    }

    # Drop any metadata columns. Note we're also dropping the duplicate `name`
    # column present in the metrics YAML.
    data %>%
        camel() %>%
        # Coerce any remaining periods (e.g. x5.3Bias).
        set_colnames(gsub("\\.", "x", colnames(.))) %>%
        # Drop blacklisted columns from the return.
        .[, sort(setdiff(colnames(.), metricsBlacklist)), drop = FALSE] %>%
        rownames_to_column() %>%
        # Ensure numerics from YAML are set correctly and not character.
        mutate_if(numericAsCharacter, as.numeric) %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        column_to_rownames()
}
