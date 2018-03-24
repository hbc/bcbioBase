#' Read Sample Metadata File
#'
#' @name readSampleMetadataFile
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param file Metadata file. Supports CSV, TSV, and XLSX file formats.
#' @param lanes *Optional*. Number of lanes used to split the samples into
#'   technical replicates (`_LXXX`) suffix.
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' # Demultiplexed
#' readSampleMetadataFile("http://bcbiobase.seq.cloud/demultiplexed.csv") %>%
#'     glimpse()
#'
#' # Multiplexed (e.g. inDrop single-cell RNA-seq)
#' readSampleMetadataFile("http://bcbiobase.seq.cloud/multiplexed.csv") %>%
#'     glimpse()
readSampleMetadataFile <- function(file, lanes = 1L) {
    assert_is_a_string(file)
    assertIsAnImplicitInteger(lanes)

    # Works with local or remote files.
    data <- readFileByExtension(file) %>%
        camel() %>%
        removeNA()

    # Don't allow the user to manually define `sampleID` column
    assert_are_disjoint_sets("sampleID", colnames(data))

    # Warn on legacy `samplename` column. This is used in some bcbio
    # documentation examples, and we need to work on improving the consistency.
    if ("samplename" %in% colnames(data)) {
        warn(paste(
            "`samplename` (note case) is used in some bcbio examples for",
            "FASTQ file names and `description` for sample names.",
            "Here we are requiring `fileName` for FASTQ file names",
            "(e.g. `control_replicate_1.fastq.gz`),",
            "`description` for demultiplexed per file sample names",
            "(e.g. `control_replicate_1`, and `sampleName` for multiplexed",
            "sample names (i.e. inDrop barcoded samples)."
        ))
        data[["fileName"]] <- data[["samplename"]]
        data[["samplename"]] <- NULL
    }

    # Check for basic required columns
    requiredCols <- c("fileName", "description")
    assert_is_subset(requiredCols, colnames(data))

    # Valid rows must contain both `fileName` and `description`
    data <- data %>%
        .[!is.na(.[["fileName"]]), , drop = FALSE] %>%
        .[!is.na(.[["description"]]), , drop = FALSE]

    # Determine whether the samples are multiplexed, based on the presence
    # of duplicate values in the `description` column
    if (
        any(duplicated(data[["fileName"]])) ||
        any(c("index", "sequence") %in% colnames(data))
    ) {
        multiplexed <- TRUE
        inform("Multiplexed samples detected")
        requiredCols <- c(requiredCols, "sampleName", "index")
        assert_is_subset(requiredCols, colnames(data))
    } else {
        multiplexed <- FALSE
        inform("Demultiplexed samples detected")
        assert_has_no_duplicates(data[["description"]])
        # Set `sampleName` column as `description` if unset
        if (!"sampleName" %in% colnames(data)) {
            data[["sampleName"]] <- data[["description"]]
        }
    }

    # Prepare metadata for lane split replicates. This step will expand rows
    # into the number of desired replicates.
    if (lanes > 1L) {
        data <- data %>%
            group_by(!!sym("description")) %>%
            # Expand by lane (e.g. "L001")
            tidyr::expand(
                lane = paste0("L", str_pad(1L:lanes, 3L, pad = "0"))
            ) %>%
            left_join(data, by = "description") %>%
            ungroup() %>%
            # Ensure `description` and `sampleName` don't contain spaces
            # upon lane expansion
            mutate(
                description = paste(
                    makeNames(!!sym("description"), unique = FALSE),
                    !!sym("lane"),
                    sep = "_"
                ),
                sampleName = paste(
                    makeNames(!!sym("sampleName"), unique = FALSE),
                    !!sym("lane"),
                    sep = "_"
                )
            )
    }

    # Ensure that `sampleName` is unique
    assert_has_no_duplicates(data[["sampleName"]])

    # Demultiplexed samples ====================================================
    if (!isTRUE(multiplexed)) {
        # Always set sampleID by the description column (simple)
        data[["sampleID"]] <- data[["description"]]
    }

    # Multiplexed samples ======================================================
    # This step is more complicated and applies to single-cell RNA-seq data.
    #
    # loadSingleCell()
    # bcbio subdirs (e.g. inDrop): `description`-`revcomp`
    # Note that forward `sequence` is required in metadata file
    # Index number is also required here for data preservation, but is not used
    # in generation of the sample directory names.
    #
    # loadCellRanger()
    # cellranger matrix: `description`-`index`
    if (isTRUE(multiplexed)) {
        # bcbio pipeline (default)
        if ("sequence" %in% colnames(data)) {
            sequence <- data[["sequence"]]
            # Require at least 6 nucleotides in the index sequence.
            # inDrop currently uses 8 but SureCell uses 6.
            assert_all_are_matching_regex(
                x = sequence,
                pattern = "^[ACGT]{6,}"
            )
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
            data[["sampleID"]] <- paste(
                data[["description"]],
                data[["revcomp"]],
                sep = "-"
            )
        } else if ("index" %in% colnames(data)) {
            # cellranger pipeline
            data[["sampleID"]] <- paste(
                data[["description"]],
                data[["index"]],
                sep = "-"
            )
        }
    }

    prepareSampleData(data)
}
