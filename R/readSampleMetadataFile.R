#' Read Sample Metadata File
#'
#' @name readSampleMetadataFile
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump camel readFileByExtension removeNA
#' @importFrom Biostrings reverseComplement
#' @importFrom dplyr group_by left_join mutate mutate_all ungroup
#' @importFrom stringr str_pad
#'
#' @inheritParams general
#'
#' @param file Metadata file. Supports CSV, TSV, and XLSX file formats.
#' @param lanes *Optional*. Number of lanes used to split the samples into
#'   technical replicates (`_LXXX`) suffix.
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' # Demultiplexed
#' readSampleMetadataFile(
#'     "http://bcbiobase.seq.cloud/demultiplexed.csv"
#' ) %>%
#'     glimpse()
#'
#' # Multiplexed (e.g. inDrop single-cell RNA-seq)
#' readSampleMetadataFile(
#'     "http://bcbiobase.seq.cloud/multiplexed.csv"
#' ) %>%
#'     glimpse()
readSampleMetadataFile <- function(file, lanes = 1L) {
    assert_is_a_string(file)
    assert_is_integer(lanes)

    # Works with local or remote files.
    data <- readFileByExtension(file)

    # Don't allow the user to manually define `sampleID` column
    assert_are_disjoint_sets("sampleID", colnames(data))

    # Warn on legacy `samplename` column. We need to work on improving the
    # consistency in examples or the internal handlng of file and sample
    # names in a future update.
    if ("samplename" %in% colnames(data)) {
        warn(paste(
            "`samplename` (note case) is used in some bcbio examples for",
            "FASTQ file names and `description` for sample names.",
            "Here we are requiring `fileName for FASTQ file names",
            "(e.g. `control_replicate_1.fastq.gz`),",
            "`description` for multiplexed per file sample names",
            "(e.g. `control replicate 1`, and `sampleName` for multiplexed",
            "sample names (i.e. inDrop barcoded samples)."
        ))
        data <- dplyr::rename(data, fileName = .data[["samplename"]])
    }

    # Check for basic required columns
    requiredCols <- c("fileName", "description")
    assert_is_subset(requiredCols, colnames(data))

    # Determine whether the samples are multiplexed, based on the presence
    # of duplicate values in the `description` column
    if (any(duplicated(data[["fileName"]])) || "index" %in% colnames(data)) {
        multiplexed <- TRUE
    } else {
        multiplexed <- FALSE
    }

    if (isTRUE(multiplexed)) {
        requiredCols <- c(requiredCols, "sampleName", "index")
        assert_is_subset(requiredCols, colnames(data))
        assert_has_no_duplicates(data[["sampleName"]])
    } else {
        assert_has_no_duplicates(data[["description"]])
        assert_are_disjoint_sets("sampleName", colnames(data))
        data[["sampleName"]] <- data[["description"]]
    }

    data <- data %>%
        # Valid rows must contain `description` and `sampleName`. Imported Excel
        # files can contain empty rows, so this helps correct that problem.
        .[!is.na(.[["description"]]), , drop = FALSE] %>%
        .[!is.na(.[["sampleName"]]), , drop = FALSE] %>%
        # Strip all NA rows and columns
        removeNA() %>%
        camel()

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
            mutate(
                description = paste(
                    .data[["description"]], .data[["lane"]], sep = "_"),
                sampleName = paste(
                    .data[["sampleName"]], .data[["lane"]], sep = "_")
            )
    }

    # This code is only applicable to multiplexed files used for single-cell
    # RNA-seq analysis. For bcbio single-cell RNA-seq, the multiplexed per
    # sample directories are created by combining the `sampleName` column
    # with the reverse complement (`revcomp`) of the index barcode sequence
    # (`sequence`). This is the current behavior for the inDrop pipeline.
    # Let's check for an ACGT sequence and use the revcomp if there's a
    # match. Otherwise just return the `sampleName` as the `sampleID`.
    if (
        isTRUE(multiplexed) &&
        is.character(data[["sequence"]])
    ) {
        detectSequence <- all(grepl("^[ACGT]{6,}", data[["sequence"]]))
        if (isTRUE(detectSequence)) {
            data[["revcomp"]] <- vapply(
                X = data[["sequence"]],
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
            # We'll sanitize into valid names using `make.names()` in
            # the final return chain.
            data[["sampleID"]] <- paste(
                data[["description"]],
                data[["revcomp"]],
                sep = "-"
            )
        }
    }

    # Default to sanitized `sampleName` column for `sampleID`
    if (!"sampleID" %in% colnames(data)) {
        data[["sampleID"]] <- data[["sampleName"]]
    }

    data %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        prepareSampleMetadata()
}
