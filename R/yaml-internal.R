.isSummaryYAML <- function(yaml) {
    ok <- is.list(yaml)
    if (!isTRUE(ok)) {
        return(FALSE)
    }

    ok <- identical(
        x = names(yaml),
        y = c("date", "upload", "bcbio_system", "samples")
    )
    if (!isTRUE(ok)) {
        return(FALSE)
    }

    TRUE
}



.sanitizeNumericAsCharacter <- function(x) {
    any(grepl(x = x, pattern = "^[0-9\\.]+$"))
}



## Currently parsing of a maximum of 2 key levels is supported
## (e.g. summary > metrics).
.sampleYAML <- function(yaml, keys) {
    assert(
        .isSummaryYAML(yaml),
        isCharacter(keys),
        isInRange(length(keys), lower = 1L, upper = 2L),
        isSubset("samples", names(yaml))
    )

    ## Focus on the sample YAML data.
    yaml <- yaml[["samples"]]
    assert(
        is.list(yaml),
        ## Note that `summary` is only returned for RNA-seq pipeline, so we're
        ## not requiring it in the assertion here.
        isSubset(
            x = c("description", "dirs", "metadata"),
            y = names(yaml[[1L]])
        )
    )

    ## Check that nested keys are present and early return on failure. Return
    ## `NULL` here instead of stopping, so we can handle bcbio RNA-seq fast mode.
    if (
        length(keys) == 2L &&
        !keys[[2L]] %in% names(yaml[[1L]][[keys[[1L]]]])
    ) {
        return(NULL)  # nocov
    }

    ## Top-level sample metadata in the YAML is relatively easy to parse. Just
    ## select for atomic values, otherwise return `NA`.
    top <- vapply(
        X = yaml,
        FUN = function(item) {
            ## Note that this will return variable length, so `vapply()` approach
            ## doesn't work here. Better method to use instead?
            return <- sapply(
                X = item,
                FUN = function(item) {
                    if (is.atomic(item)) {
                        item
                    } else {
                        NA_character_
                    }
                },
                simplify = TRUE,
                USE.NAMES = FALSE
            )
            ## Don't use `isCharacter()` assert here because we're allowing NA.
            assert(is.character(return))
            return
        },
        ## Require that there's no dimension mismatch when parsing.
        FUN.VALUE = character(length(yaml[[1L]])),
        USE.NAMES = TRUE
    )
    assert(is.matrix(top))
    top <- top %>%
        t() %>%
        as_tibble() %>%
        removeNA()
    assert(
        identical(nrow(top), length(yaml)),
        allAreAtomic(top)
    )

    ## Handle the nested metadata, defined by the keys. This step is a little
    ## tricker but should work consistently.
    nested <- lapply(
        X = yaml,
        FUN = function(item) {
            item <- item[[keys]]
            assert(is.list(item))
            ## Remove any entries that are `NULL` (e.g. "batch" in metadata).
            item <- Filter(Negate(is.null), item)
            assert(isNonEmpty(item))
            ## Sanitize names into camel case here, otherwise they'll get
            ## modified during the `ldply()` call that coerces `list` to
            ## `data.frame`.
            item <- camel(item)
            lapply(
                X = item,
                FUN = function(item) {
                    assert(is.atomic(item))
                    ## Detect and coerce nested metadata back to a string, if
                    ## necessary. bcbio allows nesting with a semicolon
                    ## delimiter.
                    if (length(item) > 1L) {
                        item <- paste(item, collapse = "; ")
                    }
                    assert(isScalar(item))
                    item
                }
            )
        }
    )
    ## Use `ldply()` method to coerce a list with uneven lengths.
    nested <- ldply(nested, data.frame, stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        removeNA()

    assert(
        isNonEmpty(nested),
        allAreAtomic(nested),
        identical(nrow(top), nrow(nested)),
        areDisjointSets(colnames(top), colnames(nested))
    )

    ## Bind the top and nested data frames, coerce to tibble, and return.
    cbind(top, nested) %>%
        as_tibble() %>%
        camel() %>%
        ## Coerce any periods in colnames to "x"
        ## (e.g. x5.3Bias becomes x5x3Bias).
        set_colnames(gsub("\\.", "x", colnames(.))) %>%
        sanitizeNA() %>%
        removeNA() %>%
        arrange(!!sym("description")) %>%
        ## Ensure numerics from YAML are set correctly and not character.
        mutate_if(.sanitizeNumericAsCharacter, as.numeric) %>%
        ## Order the columns alphabetically.
        .[, sort(colnames(.)), drop = FALSE] %>%
        ## Set the rownames.
        mutate(rowname = makeNames(!!sym("description"), unique = TRUE)) %>%
        as("DataFrame")
}
