context("readSampleMetadataFile")

test_that("Demultiplexed FASTQ", {
    file <- file.path(
        "http://bcbiobase.seq.cloud",
        "sample_metadata",
        "demultiplexed.xlsx")
    meta <- readSampleMetadataFile(file, quiet = TRUE)

    # Check that names are sanitized correctly
    rows <- c("sample_1", "sample_2", "sample_3", "sample_4")
    expect_identical(as.character(meta[["sampleID"]]), rows)
    expect_identical(rownames(meta), rows)

    # Check that column names get set correctly
    expect_identical(
        colnames(meta),
        c("sampleID", "sampleName", "description", "fileName", "genotype")
    )

    # Lane-split technical replicate support
    meta <- readSampleMetadataFile(file, lanes = 4L, quiet = TRUE)
    expect_identical(
        rownames(meta)[1L:8L],
        c(
            "sample_1_L001",
            "sample_1_L002",
            "sample_1_L003",
            "sample_1_L004",
            "sample_2_L001",
            "sample_2_L002",
            "sample_2_L003",
            "sample_2_L004")
    )
    expect_identical(
        meta[1L, metadataPriorityCols],
        data.frame(
            sampleID = factor(
                "sample_1_L001",
                levels = levels(meta[["sampleID"]])),
            sampleName = factor(
                "sample 1_L001",
                levels = levels(meta[["sampleName"]])),
            description = factor(
                "sample 1_L001",
                levels = levels(meta[["description"]])),
            row.names = "sample_1_L001")
    )

    # Error on file containing redundant `description` and `sampleName` columns
    expect_error(
        readSampleMetadataFile(
            file.path(
                "http://bcbiobase.seq.cloud",
                "sample_metadata",
                "demultiplexed_with_sampleName.csv"),
            quiet = TRUE),
        paste(
            "are_disjoint_sets :",
            "\"sampleName\" and colnames\\(data\\) have common elements:",
            "sampleName."
        )
    )

    # Required column check failure
    expect_error(
        readSampleMetadataFile(
            file.path(
                "http://bcbiobase.seq.cloud",
                "sample_metadata",
                "demultiplexed_missing_cols.csv"),
            quiet = TRUE),
        paste(
            "is_subset :",
            "The element 'description' in requiredCols is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicated description
    expect_error(
        readSampleMetadataFile(
            file.path(
                "http://bcbiobase.seq.cloud",
                "sample_metadata",
                "demultiplexed_duplicated_description.csv"),
            quiet = TRUE),
        paste(
            "has_no_duplicates :",
            "data\\[\\[\"description\"\\]\\] has a duplicate at position 2."
        )
    )
})

test_that("Multiplexed FASTQ", {
    file <- file.path(
        "http://bcbiobase.seq.cloud",
        "sample_metadata",
        "multiplexed.xlsx")
    meta <- readSampleMetadataFile(file, quiet = TRUE)

    expect_identical(
        rownames(meta),
        c(
            "run_1_CAGTTATG",
            "run_1_TTACCTCC",
            "run_2_ATAGCCTT",
            "run_2_CTTAATAG",
            "run_2_TAAGGCTC",
            "run_2_TCGCATAA",
            "run_2_TCTTACGC")
    )

    # Lane-split technical replicate support
    meta <- readSampleMetadataFile(file, lanes = 4L, quiet = TRUE)
    expect_identical(
        rownames(meta),
        c(
            "run_1_L001_CAGTTATG",
            "run_1_L001_TTACCTCC",
            "run_1_L002_CAGTTATG",
            "run_1_L002_TTACCTCC",
            "run_1_L003_CAGTTATG",
            "run_1_L003_TTACCTCC",
            "run_1_L004_CAGTTATG",
            "run_1_L004_TTACCTCC",
            "run_2_L001_ATAGCCTT",
            "run_2_L001_CTTAATAG",
            "run_2_L001_TAAGGCTC",
            "run_2_L001_TCGCATAA",
            "run_2_L001_TCTTACGC",
            "run_2_L002_ATAGCCTT",
            "run_2_L002_CTTAATAG",
            "run_2_L002_TAAGGCTC",
            "run_2_L002_TCGCATAA",
            "run_2_L002_TCTTACGC",
            "run_2_L003_ATAGCCTT",
            "run_2_L003_CTTAATAG",
            "run_2_L003_TAAGGCTC",
            "run_2_L003_TCGCATAA",
            "run_2_L003_TCTTACGC",
            "run_2_L004_ATAGCCTT",
            "run_2_L004_CTTAATAG",
            "run_2_L004_TAAGGCTC",
            "run_2_L004_TCGCATAA",
            "run_2_L004_TCTTACGC")
    )

    # Required column check failure
    expect_error(
        readSampleMetadataFile(
            file.path(
                "http://bcbiobase.seq.cloud",
                "sample_metadata",
                "multiplexed_missing_cols.csv"),
            quiet = TRUE),
        paste(
            "is_subset :",
            "The element 'index' in requiredCols is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicate rows in `sampleName` column
    expect_error(
        readSampleMetadataFile(
            file.path("http://bcbiobase.seq.cloud",
                      "sample_metadata",
                      "multiplexed_duplicated_sampleName.csv"),
            quiet = TRUE),
        paste(
            "has_no_duplicates :",
            "data\\[\\[\"sampleName\"\\]\\] has duplicates at positions 2, 4."
        )
    )
})

test_that("Legacy bcbio samplename column", {
    file <- file.path(
        "http://bcbiobase.seq.cloud",
        "sample_metadata",
        "bcbio_legacy_samplename.csv")
    meta <- suppressWarnings(readSampleMetadataFile(file, quiet = TRUE))
    expect_identical(
        meta,
        data.frame(
            # sanitized
            sampleID = factor(
                "sample_1",
                levels = "sample_1"),
            # matches description
            sampleName = factor(
                "sample-1",
                levels = "sample-1"),
            # unmodified
            description = factor(
                "sample-1",
                levels = "sample-1"),
            # renamed `samplename`
            fileName = factor(
                "sample-1.fastq.gz",
                levels = "sample-1.fastq.gz"),
            # sanitized
            row.names = "sample_1"
        )
    )
    expect_warning(
        readSampleMetadataFile(file, quiet = TRUE),
        "`samplename` \\(note case\\) is used in some bcbio examples"
    )
})

test_that("`sampleID` already defined by the user", {
    file <- file.path(
        "http://bcbiobase.seq.cloud",
        "sample_metadata",
        "sampleID_column_defined.xlsx")
    expect_error(
        readSampleMetadataFile(file, quiet = TRUE),
        paste(
            "are_disjoint_sets :",
            "\"sampleID\" and colnames\\(data\\) have common elements:",
            "sampleID."
        )
    )
})
