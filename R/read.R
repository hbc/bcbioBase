#' Read Data Versions
#'
#' @note The `data_versions.csv` file is only generated for special genomes
#' containing additional information (e.g. the built-in "hg38" build).
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "data_versions.csv")
#' x <- readDataVersions(file)
#' print(x)
readDataVersions <- function(file) {
    assert_is_a_string(file)
    # Data versions are optional.
    file <- tryCatch(
        localOrRemoteFile(file),
        error = function(e) {
            message("Data versions are missing.")
            NULL
        }
    )
    if (is.null(file)) {
        return(DataFrame())
    }
    import(file)
}



#' Read Program Versions
#'
#' @note bcbio doesn't save program versions when run in fast mode.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "programs.txt")
#' x <- readProgramVersions(file)
#' print(x)
readProgramVersions <- function(file) {
    assert_is_a_string(file)
    # Program versions are optional
    file <- tryCatch(
        localOrRemoteFile(file),
        error = function(e) {
            message("Program versions are missing.")
            NULL
        }
    )
    if (is.null(file)) {
        return(DataFrame())
    }
    # bcbio outputs `programs.txt`, but the file is comma separated.
    data <- read_csv(
        file,
        col_names = c("program", "version"),
        col_types = "cc"  # character
    )
    as(data, "DataFrame")
}



#' Read Sample Metadata
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
#' in the final upload directory using the [readYAMLSampleData()] function. If
#' you notice any typos in your metadata after completing the run, these can be
#' corrected by editing the YAML file. Alternatively, you can pass in a
#' spreadsheet with the [readSampleData()] function.
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
#' [readSampleData()] checks to see if the `description` column is unique. If
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
#' @export
#'
#' @inheritParams general
#' @param file `string`. File path. Supports CSV, TSV, and XLSX file formats.
#' @param lanes `scalar integer`. Number of lanes used to split the samples into
#'   technical replicates (`_LXXX`) suffix.
#'
#' @return `DataFrame`.
#'
#' @examples
#' ## Demultiplexed
#' file <- file.path(bcbioBaseCacheURL, "demultiplexed.csv")
#' basejump::import(file)
#' x <- readSampleData(file)
#' print(x)
#'
#' ## Multiplexed
#' file <- file.path(bcbioBaseCacheURL, "multiplexed_indrops.csv")
#' basejump::import(file)
#' x <- readSampleData(file)
#' print(x)
readSampleData <- function(file, lanes = 0L) {
    assert_is_a_string(file)
    assertIsImplicitInteger(lanes)
    assert_all_are_non_negative(lanes)
    lanes <- as.integer(lanes)

    # Convert to a sequence, if necessary.
    if (has_length(lanes, n = 1L) && lanes > 1L) {
        lanes <- seq_len(lanes)
    }

    # Works with local or remote files.
    # Ensure coercion to tibble here, for consistent handling.
    data <- import(file) %>%
        as_tibble(rownames = NULL) %>%
        camel() %>%
        removeNA()

    # Stop on input of blacklisted columns.
    intersect <- intersect(.sampleDataBlacklist, colnames(data))
    if (has_length(intersect)) {
        stop(paste0(
            paste("Invalid columns:", toString(intersect)), "\n",
            "Recommended columns:\n",
            "  - fileName: FASTQ file name (optional, but recommended).\n",
            "  - description: Sample description per file (required).\n",
            "  - sampleName: Unique sample name",
            " (multiplexed samples only).\n",
            "The `description` column is sanitized into the sample ID ",
            "for demultiplexed samples."
        ))
    }

    # Check for required columns.
    # The `description` column is always required.
    required <- "description"
    assert_is_subset(required, colnames(data))

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
        assert_is_subset(required, colnames(data))
        # Note that `description` column is expected to have duplicates.
        assertAllAreValidNames(unique(data[["description"]]))
        assert_has_no_duplicates(data[["sampleName"]])
    } else {
        multiplexed <- FALSE
        message("Demultiplexed samples detected.")
        assert_has_no_duplicates(data[["description"]])
        assertAllAreValidNames(data[["description"]])
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



#' Read Transcript-to-Gene Annotations
#'
#' Generates a `Tx2Gene` object containing `transcriptID` and `geneID` columns.
#'
#' @note Doesn't attempt to strip transcript versions.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param organism `string` or `NULL`. Full Latin organism name
#'   (e.g. `"Homo sapiens"`).
#' @param genomeBuild `string` or `NULL`. Genome build assembly name
#'   (e.g. `"GRCh38"`).
#' @param ensemblRelease `scalar integer` or `NULL`. Ensembl release version
#'   (e.g. `90`).
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "tx2gene.csv")
#' x <- readTx2Gene(
#'     file = file,
#'     organism = "Mus musculus",
#'     genomeBuild = "GRCm38",
#'     ensemblRelease = 90L
#' )
#' print(x)
readTx2Gene <- function(
    file,
    organism = NULL,
    genomeBuild = NULL,
    ensemblRelease = NULL
) {
    data <- read_csv(
        file = file,
        col_names = c("transcriptID", "geneID"),
        col_types = "cc"  # character
    )
    data <- as(data, "DataFrame")
    metadata(data) <- list(
        organism = as.character(organism),
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease)
    )
    Tx2Gene(data)
}
