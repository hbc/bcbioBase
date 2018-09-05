#' Read Sample Metadata
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param file `string`. File path. Supports CSV, TSV, and XLSX file formats.
#' @param lanes `scalar integer`. Number of lanes used to split the samples into
#'   technical replicates (`_LXXX`) suffix.
#'
#' @return `DataFrame`.
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
    # Ensure coercion to tibble here, for consistent handling.
    data <- readFileByExtension(file) %>%
        as("tbl_df") %>%
        camel() %>%
        removeNA()

    # Stop on input of blacklisted columns.
    intersect <- intersect(.sampleDataBlacklist, colnames(data))
    if (length(intersect)) {
        stop(paste0(
            paste("Invalid columns:", toString(intersect)), "\n",
            "Recommended values:\n",
            "- description: Sample name per file (required).\n",
            "- sampleName: Unique sample name",
            " (required for multiplexed samples).\n",
            "- fileName: FASTQ file name (optional, but recommended)."
        ))
    }

    # Check for required columns.
    # The `description` column is always required.
    required <- "description"
    assert_is_subset(required, colnames(data))

    # Valid rows must contain a non-empty description.
    data <- data[!is.na(data[["description"]]), , drop = FALSE]
    assertAllAreValidNames(data[["description"]])

    # Determine whether the samples are multiplexed, based on the presence
    # of duplicate values in the `description` column.
    if (
        any(duplicated(data[["description"]])) ||
        any(c("index", "sequence") %in% colnames(data))
    ) {
        multiplexed <- TRUE
        message("Multiplexed samples detected")
        required <- c(required, "sampleName", "index")
        assert_is_subset(required, colnames(data))
        assert_has_no_duplicates(data[["sampleName"]])
    } else {
        multiplexed <- FALSE
        message("Demultiplexed samples detected")
        assert_has_no_duplicates(data[["description"]])
        # Note that `sampleName` column isn't required for demultiplexed
        # samples. We can assign from the bcbio `description` automatically.
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

    # This step applies to handling single-cell metadata.
    # - bcbio subdirs (e.g. inDrops): `description`-`revcomp`.
    # - Note that forward `sequence` is required in metadata file.
    # - Index number is also required here for data preservation, but is not
    #   used in generation of the sample directory names.
    # - Require at least 6 nucleotides in the index sequence.
    # - inDrops currently uses 8 but SureCell uses 6.
    if (isTRUE(multiplexed)) {
        if ("sequence" %in% colnames(data)) {
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
    assert_is_tbl_df(data)
    assert_is_subset(
        x = "description",
        y = colnames(data)
    )
    assert_are_disjoint_sets(
        x = .sampleDataBlacklist,
        y = colnames(data)
    )

    # Set sampleName from description, if necessary
    if (!"sampleName" %in% colnames(data)) {
        data[["sampleName"]] <- data[["description"]]
    }

    # Ensure `sampleID` has valid names here. This allows for input of samples
    # beginning with numbers or containing hyphens for example, which aren't
    # valid names in R.
    data <- data %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        mutate(rowname = makeNames(!!sym("description"), unique = TRUE)) %>%
        arrange(!!sym("rowname")) %>%
        select(!!sym("sampleName"), everything()) %>%
        as("DataFrame")

    assertHasRownames(data)
    data
}



# Consider adding "rowname" here.
.sampleDataBlacklist <- c(
    "filename",  # note case: use "fileName" instead.
    "interestingGroups",
    "samplename",  # note case: use "description" instead.
    "sampleID"
)
