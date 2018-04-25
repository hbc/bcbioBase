#' Read Sample Metadata
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param file File path. Supports CSV, TSV, and XLSX file formats.
#' @param lanes Number of lanes used to split the samples into technical
#'   replicates (`_LXXX`) suffix.
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' # Demultiplexed
#' x <- readSampleData("http://bcbiobase.seq.cloud/demultiplexed.csv")
#' glimpse(x)
#'
#' # Multiplexed (e.g. inDrop single-cell RNA-seq)
#' x <- readSampleData("http://bcbiobase.seq.cloud/multiplexed.csv")
#' glimpse(x)
readSampleData <- function(file, lanes = 1L) {
    assert_is_a_string(file)
    assertIsAnImplicitInteger(lanes)
    assert_all_are_positive(lanes)

    # Works with local or remote files.
    data <- readFileByExtension(file) %>%
        camel() %>%
        removeNA()

    # Warn on legacy `samplename` column. This is used in some bcbio
    # documentation examples, and we need to work on improving the consistency.
    if (any(c("samplename", "sampleID") %in% colnames(data))) {
        warning(paste(
            "Invalid metadata columns detected.",
            "Recommended values:",
            "- description: Sample name per file (required)",
            "- sampleName: Multiplexed sample name (only for single-cell)",
            "- fileName: FASTQ file name (optional but recommended)",
            sep = "\n"
        ))
        data[["fileName"]] <- data[["samplename"]]
        data[["samplename"]] <- NULL
        data[["sampleID"]] <- NULL
    }

    # Check for description column (always required)
    requiredCols <- "description"
    assert_is_subset(requiredCols, colnames(data))

    # Valid rows must non-empty description
    data <- data[!is.na(data[["description"]]), , drop = FALSE]
    assert_is_non_empty(data)

    # Determine whether the samples are multiplexed, based on the presence
    # of duplicate values in the `description` column
    if (
        any(duplicated(data[["description"]])) ||
        any(c("index", "sequence") %in% colnames(data))
    ) {
        multiplexed <- TRUE
        message("Multiplexed samples detected")
        requiredCols <- c(requiredCols, "sampleName", "index")
        assert_is_subset(requiredCols, colnames(data))
        assert_has_no_duplicates(data[["sampleName"]])
    } else {
        multiplexed <- FALSE
        message("Demultiplexed samples detected")
        assert_has_no_duplicates(data[["description"]])
        if (!"sampleName" %in% colnames(data)) {
            data[["sampleName"]] <- data[["description"]]
        }
    }
    nameCols <- c("sampleName", "description")

    # Prepare metadata for lane split replicates. This step will expand rows
    # into the number of desired replicates.
    if (lanes > 1L) {
        data <- data %>%
            group_by(!!!syms(nameCols)) %>%
            # Expand by lane (e.g. "L001")
            expand(
                lane = paste0("L", str_pad(1L:lanes, 3L, pad = "0"))
            ) %>%
            left_join(data, by = nameCols) %>%
            ungroup() %>%
            # Ensure lane-split metadata doesn't contain spaces
            mutate_at(
                nameCols,
                funs(
                    paste(
                        makeNames(., unique = FALSE),
                        !!sym("lane"),
                        sep = "_"
                    )
                )
            )
    }

    # This step applies to handling single-cell metadata
    if (isTRUE(multiplexed)) {
        if ("sequence" %in% colnames(data)) {
            # bcbio subdirs (e.g. inDrop): `description`-`revcomp` Note that
            # forward `sequence` is required in metadata file. Index number is
            # also required here for data preservation, but is not used in
            # generation of the sample directory names. Require at least 6
            # nucleotides in the index sequence. inDrop currently uses 8 but
            # SureCell uses 6.
            sequence <- data[["sequence"]]
            assert_all_are_matching_regex(sequence, "^[ACGT]{6,}")
            data[["revcomp"]] <- vapply(
                X = sequence,
                FUN = function(x) {
                    x %>%
                        as("character") %>%
                        as("DNAStringSet") %>%
                        reverseComplement() %>%
                        as("character")
                },
                FUN.VALUE = character(1L)
            )
            # Match the sample directories exactly here, using the hyphen
            data[["description"]] <- paste(
                data[["description"]],
                data[["revcomp"]],
                sep = "-"
            )
        } else if ("index" %in% colnames(data)) {
            # CellRanger: `description`-`index`
            data[["description"]] <- paste(
                data[["description"]],
                data[["index"]],
                sep = "-"
            )
        }
    }

    .returnSampleData(data)
}



# Consistent sanitization for YAML and external file
.returnSampleData <- function(data) {
    assert_has_dimnames(data)
    assert_is_subset("description", colnames(data))
    assert_are_disjoint_sets("sampleID", colnames(data))

    # Set sampleName from description, if necessary
    if (!"sampleName" %in% colnames(data)) {
        data[["sampleName"]] <- data[["description"]]
    }

    data <- data %>%
        # Ensure `sampleID` has valid names. This allows for input of samples
        # beginning with numbers or containing hyphens for example, which aren't
        # valid names in R.
        mutate(rowname = makeNames(!!sym("description"), unique = TRUE)) %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        arrange(!!sym("rowname")) %>%
        select(!!sym("sampleName"), everything()) %>%
        as.data.frame() %>%
        column_to_rownames()

    data
}
