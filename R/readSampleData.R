#' Read sample metadata
#'
#' This function reads user-defined sample metadata saved in a spreadsheet.
#' The "`description`" column is always required, and must match the bcbio
#' per sample directory names exactly. Inclusion of the "`fileName`" column
#' isn't required but is recommended for data provenance.
#'
#' @note Some bcbio examples on readthedocs use "`samplename`" instead of
#'   "`fileName`". This function checks for that and will error out
#'   intentionally, since we're using the `sampleName` column (note case) to
#'   define unique sample names, in the event that bcbio has processed
#'   multiplexed samples.
#'
#' @section Demultiplexed samples:
#'
#' This applies to bulk RNA-seq samples.
#'
#' Normally when loading a bcbio run of demultiplexed samples, the sample
#' metadata will be imported automatically from the `project-summary.yaml` file
#' in the final upload directory using the `readYAMLSampleData` function. If
#' you notice any typos in your metadata after completing the run, these can be
#' corrected by editing the YAML file. Alternatively, you can pass in a
#' spreadsheet with the `readSampleData` function.
#'
#' The samples in the bcbio run must map to the `description` column. The values
#' provided in `description` for demultiplexed samples must be unique. They must
#' also be *syntactically valid*, meaning that they cannot contain illegal
#' characters (e.g. spaces, non-alphanumerics, *dashes*) or *begin with a
#' number*. Consult the documentation in `help(topic = "make.names")` for more
#' information on valid names in R.
#'
#' @section Multiplexed samples:
#'
#' This applies to some single-cell RNA-seq formats, including inDrops. In this
#' case, bcbio will output per-sample directories with this this structure:
#' "`description`-`revcomp`".
#'
#' `readSampleData` checks to see if the `description` column is unique. If
#' the values are duplicated, the function assumes that bcbio processed
#' multiplexed FASTQs, where multiple samples of interest are barcoded inside a
#' single FASTQ. This this case, you must supply additional "`index`",
#' "`sequence`", and "`sampleName`" columns.
#'
#' Note that bcbio currently outputs the reverse complement index sequence in
#' the sample directory names (e.g. "`sample-ATAGAGAG`"). Define the forward
#' index barcode in the `sequence` column here, not the reverse complement. The
#' reverse complement will be calculated automatically and added as the
#' `revcomp` column in the sample metadata.
#'
#' @author Michael Steinbaugh
#' @inheritParams basejump::params
#' @export
#'
#' @param file `character(1)`.
#'   File path. Supports CSV, TSV, and XLSX file formats.
#' @param lanes `integer(1)`.
#'   Number of lanes used to split the samples into technical replicates
#'   suffix (i.e. `_LXXX`).
#'
#' @return `DataFrame`.
#'
#' @examples
#' ## Demultiplexed
#' file <- file.path(bcbioBaseCacheURL, "demultiplexed.csv")
#' x <- readSampleData(file)
#' print(x)
#'
#' ## Multiplexed
#' file <- file.path(bcbioBaseCacheURL, "multiplexed_indrops.csv")
#' x <- readSampleData(file)
#' print(x)
readSampleData <- function(file, lanes = 0L) {
    # Coerce detectLanes empty integer return to 0.
    if (length(lanes) == 0L) {
        lanes <- 0L
    }
    # Note that we're allowing import from URL here (primarily for unit tests).
    assert(
        isAFile(file) || containsAURL(file),
        isInt(lanes),
        isNonNegative(lanes)
    )
    lanes <- as.integer(lanes)

    # Convert to a sequence, if necessary.
    if (
        length(lanes) == 1L &&
        lanes > 1L
    ) {
        lanes <- seq_len(lanes)
    }

    # Works with local or remote files.
    # Ensure coercion to tibble here, for consistent handling.
    data <- import(file) %>%
        as_tibble(rownames = NULL) %>%
        camel() %>%
        removeNA()

    # Check to ensure that columns are valid, before proceeding.
    assert(.isSampleData(data))

    # Check for required columns. The `description` column is always required.
    required <- "description"
    assert(isSubset(required, colnames(data)))

    # Valid rows must contain a non-empty description.
    data <- data[!is.na(data[["description"]]), , drop = FALSE]

    # Determine whether the samples are multiplexed, based on the presence
    # of duplicate values in the `description` column.
    if (
        any(duplicated(data[["description"]])) ||
        any(c("index", "sequence") %in% colnames(data))
    ) {
        multiplexed <- TRUE
        message("Multiplexed samples detected.")
        required <- c(required, "sampleName", "index")
        assert(isSubset(required, colnames(data)))
        # Note that `description` column is expected to have duplicates.
        assert(
            validNames(unique(data[["description"]])),
            hasNoDuplicates(data[["sampleName"]])
        )
    } else {
        multiplexed <- FALSE
        message("Demultiplexed samples detected.")
        assert(
            hasNoDuplicates(data[["description"]]),
            validNames(data[["description"]])
        )

        # Note that `sampleName` column isn't required for demultiplexed
        # samples. We can assign from the bcbio `description` automatically.
        if (!"sampleName" %in% colnames(data)) {
            data[["sampleName"]] <- data[["description"]]
        }
    }
    nameCols <- c("sampleName", "description")

    # Prepare metadata for lane split replicates. This step will expand rows
    # into the number of desired replicates.
    if (length(lanes) > 1L) {
        data <- data %>%
            group_by(!!!syms(nameCols)) %>%
            # Expand by lane (e.g. "L001").
            expand(
                lane = paste0(
                    "L", str_pad(string = lanes, width = 3L, pad = "0")
                )
            ) %>%
            left_join(data, by = nameCols) %>%
            ungroup() %>%
            # Ensure lane-split metadata doesn't contain spaces.
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
        assert(isSubset(c("index", "sequence"), colnames(data)))
        sequence <- data[["sequence"]]
        assert(allAreMatchingRegex(sequence, pattern = "^[ACGT]{6,}"))
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
        # Match the sample directories exactly here, using the hyphen.
        data[["description"]] <- paste(
            data[["description"]],
            data[["revcomp"]],
            sep = "-"
        )
    }

    .makeSampleData(data)
}
