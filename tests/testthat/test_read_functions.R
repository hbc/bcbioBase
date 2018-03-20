context("Read Functions")

# readDataVersions =============================================================
test_that("readDataVersions", {
    versions <- readDataVersions("data_versions.csv")
    expect_is(versions, "tbl_df")
    expect_identical(
        colnames(versions),
        c("genome", "resource", "version")
    )
})

test_that("readDataVersions : Missing file", {
    expect_warning(
        readDataVersions("XXX.csv"),
        "is_existing_file :"
    )
    expect_identical(
        suppressWarnings(
            readDataVersions("XXX.csv")
        ),
        tibble::tibble()
    )
})



# readLogFile ==================================================================
test_that("readLogFile", {
    log <- readLogFile("bcbio-nextgen.log")
    expect_true(is.character(log))
    expect_identical(
        log[[1L]],
        paste(
            "[2017-08-15T14:53Z]",
            "compute-a-16-44.o2.rc.hms.harvard.edu:",
            "System YAML configuration:",
            "/n/app/bcbio/dev/galaxy/bcbio_system.yaml"
        )
    )
})

test_that("readLogFile : Missing file", {
    expect_error(
        readLogFile("XXX.log"),
        "is_existing_file :"
    )
})



# readProgramVersions ==========================================================
test_that("readProgramVersions", {
    versions <- readProgramVersions("programs.txt")
    expect_is(versions, "tbl_df")
    expect_identical(
        colnames(versions),
        c("program", "version")
    )
})

test_that("readProgramVersions : Missing file", {
    expect_warning(
        readProgramVersions("XXX.txt"),
        "is_existing_file :"
    )
    expect_identical(
        suppressWarnings(
            readProgramVersions("XXX.txt")
        ),
        tibble::tibble()
    )
})



# readSampleMetadataFile =======================================================
test_that("readSampleMetadataFile : Demultiplexed FASTQ", {
    file <- "demultiplexed.csv"
    meta <- readSampleMetadataFile(file)

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
    meta <- readSampleMetadataFile(file, lanes = 4L)
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
            "sample_2_L004"
        )
    )
    expect_identical(
        meta[1L, metadataPriorityCols],
        data.frame(
            sampleID = factor(
                "sample_1_L001",
                levels = levels(meta[["sampleID"]])
            ),
            sampleName = factor(
                "sample 1_L001",
                levels = levels(meta[["sampleName"]])
            ),
            description = factor(
                "sample 1_L001",
                levels = levels(meta[["description"]])
            ),
            row.names = "sample_1_L001"
        )
    )

    # Error on file containing redundant `description` and `sampleName` columns
    expect_error(
        readSampleMetadataFile("demultiplexed_with_sampleName.csv"),
        paste(
            "are_disjoint_sets :",
            "\"sampleName\" and colnames\\(data\\) have common elements:",
            "sampleName."
        )
    )

    # Required column check failure
    expect_error(
        readSampleMetadataFile("demultiplexed_missing_cols.csv"),
        paste(
            "is_subset :",
            "The element 'description' in requiredCols is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicated description
    expect_error(
        readSampleMetadataFile("demultiplexed_duplicated_description.csv"),
        paste(
            "has_no_duplicates :",
            "data\\[\\[\"description\"\\]\\] has a duplicate at position 2."
        )
    )
})

test_that("readSampleMetadataFile : Multiplexed FASTQ", {
    file <- "multiplexed.csv"
    meta <- readSampleMetadataFile(file)
    expect_identical(
        rownames(meta),
        c(
            "run_1_CAGTTATG",
            "run_1_TTACCTCC",
            "run_2_ATAGCCTT",
            "run_2_CTTAATAG",
            "run_2_TAAGGCTC",
            "run_2_TCGCATAA",
            "run_2_TCTTACGC"
        )
    )

    # Lane-split technical replicate support
    meta <- readSampleMetadataFile(file, lanes = 4L)
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
            "run_2_L004_TCTTACGC"
        )
    )

    # Required column check failure
    expect_error(
        readSampleMetadataFile("multiplexed_missing_cols.csv"),
        paste(
            "is_subset :",
            "The element 'index' in requiredCols is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicate rows in `sampleName` column
    expect_error(
        readSampleMetadataFile("multiplexed_duplicated_sampleName.csv"),
        paste(
            "has_no_duplicates :",
            "data\\[\\[\"sampleName\"\\]\\] has duplicates at positions 2, 4."
        )
    )
})

test_that("readSampleMetadataFile : Legacy bcbio samplename column", {
    file <- "bcbio_legacy_samplename.csv"
    meta <- suppressWarnings(readSampleMetadataFile(file))
    expect_identical(
        meta,
        data.frame(
            # sanitized
            sampleID = factor(
                "sample_1",
                levels = "sample_1"
            ),
            # matches description
            sampleName = factor(
                "sample-1",
                levels = "sample-1"
            ),
            # unmodified
            description = factor(
                "sample-1",
                levels = "sample-1"
            ),
            # renamed `samplename`
            fileName = factor(
                "sample-1.fastq.gz",
                levels = "sample-1.fastq.gz"
            ),
            # sanitized
            row.names = "sample_1"
        )
    )
    expect_warning(
        readSampleMetadataFile(file),
        "`samplename` \\(note case\\) is used in some bcbio examples"
    )
})

test_that("readSampleMetadataFile : sampleID defined by user", {
    expect_error(
        readSampleMetadataFile("sampleID_column_defined.csv"),
        paste(
            "are_disjoint_sets :",
            "\"sampleID\" and colnames\\(data\\) have common elements:",
            "sampleID."
        )
    )
})

test_that("readSampleMetadataFile : Missing file", {
    # Always stop on missing
    expect_error(
        readSampleMetadataFile("XXX.csv"),
        "is_existing_file :"
    )
})
