#' YAML Parsing Functions
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @name yaml
#' @author Michael Steinbaugh
#' @inheritParams basejump::params
#'
#' @param yaml `list`. Project summary YAML.
#'
#' @return `string` or `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "summary.yaml")
#' yaml <- basejump::import(file)
#' summary(yaml)
#'
#' ## GTF file path.
#' x <- getGTFFileFromYAML(yaml)
#' print(x)
#'
#' ## Sample metrics.
#' x <- getMetricsFromYAML(yaml)
#' summary(x)
#' colnames(x)
#'
#' ## Sample metadata.
#' x <- getSampleDataFromYAML(yaml)
#' summary(x)
#' colnames(x)
NULL



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



# GTF file =====================================================================
#' @describeIn yaml `string`. GTF file path.
#' @export
getGTFFileFromYAML <- function(yaml) {
    .assertIsSummaryYAML(yaml)
    # Assume all samples are using the same GTF file.
    file <- yaml %>%
        .[["samples"]] %>%
        .[[1L]] %>%
        .[["genome_resources"]] %>%
        .[["rnaseq"]] %>%
        .[["transcripts"]]
    assert_is_a_string(file)
    assert_are_identical(basename(file), "ref-transcripts.gtf")
    file
}



# Sample-level YAML information ================================================
# Currently max 2 keys are supported (e.g. summary, metrics).
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
        has_length(keys, n = 2L) &&
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



#' @describeIn yaml `DataFrame`. Quality control summary metrics.
#' @export
getMetricsFromYAML <- function(yaml) {
    message("Getting sample metrics from YAML.")
    data <- .sampleYAML(yaml, keys = c("summary", "metrics"))

    # Early return on empty metrics (e.g. fast mode).
    if (!has_length(data)) {
        # nocov start
        message("No metrics were calculated.")
        return(NULL)
        # nocov end
    }

    # Drop any metadata columns. Note we're also dropping the duplicate `name`
    # column present in the metrics YAML.
    yamlFlatCols <- c("description", "genome_build", "sam_ref")
    blacklist <- c(camel(yamlFlatCols), "name")
    data <- data %>%
        as_tibble(rownames = "rowname") %>%
        # Drop blacklisted columns from the return.
        .[, sort(setdiff(colnames(.), blacklist)), drop = FALSE] %>%
        # Convert all strings to factors.
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        as("DataFrame")
    assertHasRownames(data)
    data
}



#' @describeIn yaml `DataFrame`. Sample metadata.
#' @export
getSampleDataFromYAML <- function(yaml) {
    message("Making sample metadata from YAML.")
    yaml %>%
        .sampleYAML(keys = "metadata") %>%
        .makeSampleData()
}
