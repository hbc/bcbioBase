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

test_that("readDataVersions : Silent on missing file", {
    expect_identical(
        readDataVersions("XXX.csv"),
        tibble::tibble()
    )
})



# readLog ==================================================================
test_that("readLog", {
    log <- readLog("bcbio-nextgen.log")
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

test_that("readLog : Missing file", {
    expect_error(
        readLog("XXX.log"),
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

test_that("readProgramVersions : Silent on missing file", {
    expect_identical(
        readProgramVersions("XXX.txt"),
        tibble::tibble()
    )
})



# readSampleData =======================================================
test_that("readSampleData : Demultiplexed FASTQ", {
    file <- "demultiplexed.csv"
    x <- readSampleData(file)

    # Check that names are sanitized correctly
    rows <- c("sample_1", "sample_2", "sample_3", "sample_4")
    expect_identical(as.character(x[["sampleID"]]), rows)
    expect_identical(rownames(x), rows)

    # Check that column names get set correctly
    expect_identical(
        colnames(x),
        c("sampleID", "sampleName", "description", "fileName", "genotype")
    )

    # Lane-split technical replicate support
    x <- readSampleData(file, lanes = 4L)
    expect_identical(
        rownames(x)[1L:8L],
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
    # Check that makeNames is applied correctly
    expect_identical(
        x[1L, metadataPriorityCols],
        data.frame(
            sampleID = factor(
                "sample_1_L001",
                levels = levels(x[["sampleID"]])
            ),
            sampleName = factor(
                "sample_1_L001",
                levels = levels(x[["sampleName"]])
            ),
            description = factor(
                "sample_1_L001",
                levels = levels(x[["description"]])
            ),
            row.names = "sample_1_L001"
        )
    )

    # Required column check failure
    expect_error(
        readSampleData("demultiplexed_missing_cols.csv"),
        paste(
            "is_subset :",
            "The element 'description' in requiredCols is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicated description
    expect_error(
        readSampleData("demultiplexed_duplicated_description.csv"),
        paste(
            "has_no_duplicates :",
            "data\\[\\[\"description\"\\]\\] has a duplicate at position 2."
        )
    )
})

test_that("readSampleData : Multiplexed FASTQ", {
    file <- "multiplexed.csv"
    meta <- readSampleData(file)
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
    meta <- readSampleData(file, lanes = 4L)
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
        readSampleData("multiplexed_missing_cols.csv"),
        paste(
            "is_subset :",
            "The element 'index' in requiredCols is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicate rows in `sampleName` column
    expect_error(
        readSampleData("multiplexed_duplicated_sampleName.csv"),
        paste(
            "has_no_duplicates :",
            "data\\[\\[\"sampleName\"\\]\\] has duplicates at positions 2, 4."
        )
    )
})

test_that("readSampleData : CellRanger pipeline", {
    x <- readSampleData("cellranger_metadata.csv")
    y <- data.frame(
        "sampleID" = factor(c("aggregation_1", "aggregation_2")),
        "sampleName" = factor(c("proximal", "distal")),
        "description" = factor("aggregation"),
        "fileName" = factor("aggregation.fastq.gz"),
        "index" = factor(c(1L, 2L)),
        row.names = c("aggregation_1", "aggregation_2")
    )
    assert_are_identical(x, y)
})

test_that("readSampleData : Legacy bcbio samplename column", {
    file <- "bcbio_legacy_samplename.csv"
    meta <- suppressWarnings(readSampleData(file))
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
        readSampleData(file),
        "`samplename` \\(note case\\) is used in some bcbio examples"
    )
})

test_that("readSampleData : sampleID defined by user", {
    expect_error(
        readSampleData("sampleID_column_defined.csv"),
        paste(
            "are_disjoint_sets :",
            "\"sampleID\" and colnames\\(data\\) have common elements:",
            "sampleID."
        )
    )
})

test_that("readSampleData : Missing file", {
    # Always stop on missing
    expect_error(
        readSampleData("XXX.csv"),
        "is_existing_file :"
    )
})



# readTx2gene ==================================================================
test_that("readTx2gene", {
    x <- readTx2gene("tx2gene.csv")
    expect_is(x, "data.frame")
    expect_identical(
        colnames(x),
        c("txID", "geneID")
    )
})
