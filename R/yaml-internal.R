## Updated 2019-07-23.
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



## Currently parsing of a maximum of 2 key levels is supported.
## (e.g. summary > metrics).
## Updated 2019-08-20.
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
    ## `NULL` here instead of stopping, so we can handle bcbio RNA-seq fast
    ## mode.
    if (
        hasLength(keys, n = 2L) &&
        !isSubset(x = keys[[2L]], y = names(yaml[[1L]][[keys[[1L]]]]))
    ) {
        return(NULL)  # nocov
    }
    ## Top-level sample metadata in the YAML is relatively easy to parse. Just
    ## select for atomic values, otherwise return `NA`.
    top <- vapply(
        X = yaml,
        FUN = function(item) {
            ## Note that this will return variable length, so `vapply()`
            ## approach doesn't work here. Better method to use instead?
            x <- sapply(
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
            assert(is.character(x))
            x
        },
        ## Require that there's no dimension mismatch when parsing.
        FUN.VALUE = character(length(yaml[[1L]])),
        USE.NAMES = TRUE
    )
    assert(is.matrix(top))
    top <- t(top)
    top <- as(top, "DataFrame")
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
            item <- camelCase(item)
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
    ## Using `ldply()` method here to coerce a list with uneven lengths.
    nested <- ldply(
        .data = nested,
        .fun = data.frame,
        stringsAsFactors = FALSE
    )
    nested <- as(nested, "DataFrame")
    assert(
        isNonEmpty(nested),
        allAreAtomic(nested),
        identical(nrow(top), nrow(nested)),
        areDisjointSets(colnames(top), colnames(nested))
    )
    out <- cbind(top, nested)
    out <- camelCase(out)
    ## Coerce any periods in colnames to "x" (e.g. `x5.3Bias` to `x5x3Bias`).
    colnames(out) <- gsub("\\.", "x", colnames(out))
    out <- sanitizeNA(out)
    out <- removeNA(out)
    out <- out[order(out[["description"]]), , drop = FALSE]
    ## Ensure numerics from YAML are set correctly and not character.
    out <- DataFrame(lapply(
        X = out,
        FUN = function(x) {
            if (
                is.character(x) &&
                any(grepl(x = x, pattern = "^[0-9\\.]+$"))
            ) {
                x <- as.numeric(x)
            }
            x
        }
    ))
    ## Order the columns alphabetically.
    out <- out[, sort(colnames(out)), drop = FALSE]
    ## Set the rownames.
    rownames(out) <- makeNames(out[["description"]], unique = TRUE)
    out
}
