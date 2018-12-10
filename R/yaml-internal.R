.assertIsSummaryYAML <- function(yaml) {
    assert_is_list(yaml)
    assert_is_non_empty(yaml)
    assert_are_identical(
        x = names(yaml),
        y = c("date", "upload", "bcbio_system", "samples")
    )
}



.sanitizeNumericAsCharacter <- function(x) {
    any(grepl(x = x, pattern = "^[0-9\\.]+$"))
}



# Currently parsing of a maximum of 2 key levels is supported
# (e.g. summary>metrics).
.sampleYAML <- function(yaml, keys) {
    .assertIsSummaryYAML(yaml)
    assert_is_character(keys)
    assert_all_are_in_range(length(keys), lower = 1L, upper = 2L)

    # Focus on the sample YAML data.
    assert_is_subset("samples", names(yaml))
    yaml <- yaml[["samples"]]
    assert_is_list(yaml)
    assert_is_non_empty(yaml)

    # Note that `summary` is only returned for RNA-seq pipeline.
    assert_is_subset(
        x = c("description", "dirs", "metadata"),
        y = names(yaml[[1L]])
    )

    # Check that nested keys are present and early return on failure.
    # Return `NULL` here instead of stopping, so we can handle bcbio RNA-seq
    # fast mode runs.
    if (
        length(keys) == 2L &&
        !keys[[2L]] %in% names(yaml[[1L]][[keys[[1L]]]])
    ) {
        return(NULL)  # nocov
    }

    # Top-level sample metadata in the YAML is relatively easy to parse.
    # Just select for atomic values, otherwise return `NA`.
    top <- vapply(
        X = yaml,
        FUN = function(item) {
            return <- sapply(
                X = item,
                FUN = function(item) {
                    if (is.atomic(item)) {
                        item
                    } else {
                        NA
                    }
                },
                simplify = TRUE,
                USE.NAMES = FALSE
            )
            assert_is_character(return)
            return
        },
        # Require that there's no dimension mismatch when parsing.
        FUN.VALUE = character(length(yaml[[1L]])),
        USE.NAMES = TRUE
    )
    assert_is_matrix(top)
    top <- top %>%
        t() %>%
        as_tibble() %>%
        removeNA()
    assert_are_identical(nrow(top), length(yaml))
    invisible(lapply(top, assert_is_atomic))

    # Handle the nested metadata, defined by the keys.
    # This step is a little tricker but should work consistently.
    nested <- lapply(
        X = yaml,
        FUN = function(item) {
            item <- item[[keys]]
            assert_is_list(item)
            # Remove any entries that are `NULL` (e.g. "batch" in metadata).
            item <- Filter(Negate(is.null), item)
            assert_is_non_empty(item)
            # Sanitize names into camel case here, otherwise they'll get
            # modified during the `ldply()` call that coerces `list` to
            # `data.frame`.
            item <- camel(item)
            lapply(
                X = item,
                FUN = function(item) {
                    assert_is_atomic(item)
                    # Detect and coerce nested metadata back to a string, if
                    # necessary. bcbio allows nesting with a semicolon
                    # delimiter.
                    if (length(item) > 1L) {
                        item <- paste(item, collapse = "; ")
                    }
                    assert_is_scalar(item)
                    item
                }
            )
        }
    )
    # Use `ldply()` method to coerce a list with uneven lengths.
    nested <- ldply(nested, data.frame, stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        removeNA()
    assert_is_non_empty(nested)
    invisible(lapply(nested, assert_is_atomic))

    # Bind the top and nested data frames, coerce to tibble, and return.
    assert_are_identical(x = nrow(top), y = nrow(nested))
    assert_are_disjoint_sets(x = colnames(top), y = colnames(nested))
    cbind(top, nested) %>%
        as_tibble() %>%
        camel() %>%
        # Coerce any periods in colnames to "x"
        # (e.g. x5.3Bias becomes x5x3Bias).
        set_colnames(gsub("\\.", "x", colnames(.))) %>%
        sanitizeNA() %>%
        removeNA() %>%
        arrange(!!sym("description")) %>%
        # Ensure numerics from YAML are set correctly and not character.
        mutate_if(.sanitizeNumericAsCharacter, as.numeric) %>%
        # Order the columns alphabetically.
        .[, sort(colnames(.)), drop = FALSE] %>%
        # Set the rownames.
        mutate(rowname = makeNames(!!sym("description"), unique = TRUE)) %>%
        as("DataFrame")
}
