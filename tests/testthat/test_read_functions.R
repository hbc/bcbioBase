context("Read Functions")



# readDataVersions =============================================================
test_that("readDataVersions", {
    versions <- readDataVersions("data_versions.csv")
    expect_is(versions, "tbl_df")
    expect_identical(
        object = colnames(versions),
        expected = c("genome", "resource", "version")
    )
})

test_that("readDataVersions : Silent on missing file", {
    expect_identical(
        object = readDataVersions("XXX.csv"),
        expected = tibble::tibble()
    )
})



# readLog ==================================================================
test_that("readLog", {
    log <- readLog("bcbio-nextgen.log")
    expect_true(is.character(log))
    expect_identical(
        object = log[[1L]],
        expected = paste(
            "[2017-08-15T14:53Z]",
            "compute-a-16-44.o2.rc.hms.harvard.edu:",
            "System YAML configuration:",
            "/n/app/bcbio/dev/galaxy/bcbio_system.yaml"
        )
    )
})

test_that("readLog : Missing file", {
    expect_error(
        object = readLog("XXX.log"),
        regexp = "is_existing_file :"
    )
})



# readProgramVersions ==========================================================
test_that("readProgramVersions", {
    versions <- readProgramVersions("programs.txt")
    expect_is(versions, "tbl_df")
    expect_identical(
        object = colnames(versions),
        expected = c("program", "version")
    )
})

test_that("readProgramVersions : Silent on missing file", {
    expect_identical(
        object = readProgramVersions("XXX.txt"),
        expected = tibble::tibble()
    )
})



# readSampleData =======================================================
test_that("readSampleData : Demultiplexed FASTQ", {
    file <- "demultiplexed.csv"
    object <- readSampleData(file)

    # Check that names are sanitized correctly.
    expect_identical(
        object = rownames(object),
        expected = paste0("sample", seq_len(4L))
    )

    # Check that column names get set correctly.
    expect_identical(
        object = colnames(object),
        expected = c("sampleName", "fileName", "description", "genotype")
    )

    # Lane-split technical replicate support.
    object <- readSampleData(file, lanes = 4L)
    expect_identical(
        object = colnames(object),
        expected = c("sampleName", "description", "lane", "fileName", "genotype")
    )
    expect_identical(
        object = rownames(object)[1L:8L],
        expected = c(
            "sample1_L001",
            "sample1_L002",
            "sample1_L003",
            "sample1_L004",
            "sample2_L001",
            "sample2_L002",
            "sample2_L003",
            "sample2_L004"
        )
    )

    # Required column check failure.
    expect_error(
        object = readSampleData("demultiplexed_missing_cols.csv"),
        regexp = paste(
            "is_subset :",
            "The element 'description'"
        )
    )

    # Duplicated description.
    expect_error(
        object = readSampleData("demultiplexed_duplicated_description.csv"),
        regexp = paste(
            "is_subset :",
            "The elements 'sampleName', 'index'"
        )
    )
})

test_that("readSampleData : Multiplexed FASTQ", {
    file <- "multiplexed.csv"

    object <- readSampleData(file)
    expect_identical(
        object = rownames(object),
        expected = c(
            "run1_CAGTTATG",
            "run1_TTACCTCC",
            "run2_ATAGCCTT",
            "run2_CTTAATAG",
            "run2_TAAGGCTC",
            "run2_TCGCATAA",
            "run2_TCTTACGC"
        )
    )

    # Lane-split technical replicate support.
    object <- readSampleData(file, lanes = 4L)
    expect_identical(
        object = rownames(object),
        expected = c(
            "run1_L001_CAGTTATG",
            "run1_L001_TTACCTCC",
            "run1_L002_CAGTTATG",
            "run1_L002_TTACCTCC",
            "run1_L003_CAGTTATG",
            "run1_L003_TTACCTCC",
            "run1_L004_CAGTTATG",
            "run1_L004_TTACCTCC",
            "run2_L001_ATAGCCTT",
            "run2_L001_CTTAATAG",
            "run2_L001_TAAGGCTC",
            "run2_L001_TCGCATAA",
            "run2_L001_TCTTACGC",
            "run2_L002_ATAGCCTT",
            "run2_L002_CTTAATAG",
            "run2_L002_TAAGGCTC",
            "run2_L002_TCGCATAA",
            "run2_L002_TCTTACGC",
            "run2_L003_ATAGCCTT",
            "run2_L003_CTTAATAG",
            "run2_L003_TAAGGCTC",
            "run2_L003_TCGCATAA",
            "run2_L003_TCTTACGC",
            "run2_L004_ATAGCCTT",
            "run2_L004_CTTAATAG",
            "run2_L004_TAAGGCTC",
            "run2_L004_TCGCATAA",
            "run2_L004_TCTTACGC"
        )
    )

    # Required column check failure.
    expect_error(
        object = readSampleData("multiplexed_missing_cols.csv"),
        expected = paste(
            "is_subset :",
            "The element 'index' in required is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicate rows in `sampleName` column.
    expect_error(
        object = readSampleData("multiplexed_duplicated_sampleName.csv"),
        regexp = paste(
            "has_no_duplicates :",
            "data\\[\\[\"sampleName\"\\]\\] has duplicates at positions 2, 4."
        )
    )
})

test_that("readSampleData : Multiplexed CellRanger data", {
    object <- readSampleData("cellranger_metadata.csv")
    expected <- DataFrame(
        sampleName = factor(c("proximal", "distal")),
        fileName = factor("aggregation.fastq.gz"),
        description = factor(c("aggregation-1", "aggregation-2")),
        index = factor(c("1", "2")),
        row.names = factor(c("aggregation_1", "aggregation_2"))
    )
    expect_identical(object, expected)
})

test_that("readSampleData : Legacy bcbio `samplename` column", {
    expect_error(
        object = readSampleData("bcbio_legacy_samplename.csv"),
        regexp = "Invalid columns: samplename"
    )
})

test_that("readSampleData : sampleID defined by user", {
    expect_error(
        object = readSampleData("sampleID_column_defined.csv"),
        regexp = "Invalid columns: sampleID"
    )
})

test_that("readSampleData : Missing file", {
    expect_error(
        object = readSampleData("XXX.csv"),
        regexp = "is_existing_file :"
    )
})



# readTx2gene ==================================================================
test_that("readTx2gene", {
    object <- readTx2gene("tx2gene.csv")
    expect_is(object, "tx2gene")
    expect_identical(
        object = colnames(object),
        expected = c("transcriptID", "geneID")
    )
})



# readYAMLSampleData ===========================================================
test_that("readYAMLSampleData", {
    object <- readYAMLSampleData("project-summary.yaml")
    samples <- c("group1_1", "group1_2", "group2_1", "group2_2")
    expect_identical(
        object = object,
        DataFrame(
            sampleName = factor(samples),
            description = factor(samples),
            genomeBuild = factor("mm10"),
            group = factor(c("ctrl", "ctrl", "ko", "ko")),
            samRef = factor(paste(
                "",
                "groups",
                "bcbio",
                "bcbio_dev",
                "genomes",
                "Mmusculus",
                "mm10",
                "seq",
                "mm10.fa",
                sep = "/"
            )),
            row.names = samples
        )
    )
})

# Testing against Kayleigh's example here.
test_that("readYAMLSampleData : nested metadata", {
    object <- suppressWarnings(
        readYAMLSampleData("project-summary-nested-metadata.yaml")
    )
    expect_is(object, "DataFrame")
})



# readYAMLSampleMetrics ========================================================
test_that("readYAMLSampleMetrics", {
    expected <- list(
        averageInsertSize = "numeric",
        duplicates = "numeric",
        duplicationRateOfMapped = "numeric",
        exonicRate = "numeric",
        intergenicRate = "numeric",
        intronicRate = "numeric",
        mappedPairedReads = "numeric",
        mappedReads = "numeric",
        percentGC = "numeric",
        qualityFormat = "factor",
        rrna = "numeric",
        rrnaRate = "numeric",
        sequenceLength = "factor",
        sequencesFlaggedAsPoorQuality = "numeric",
        totalReads = "numeric",
        x5x3Bias = "numeric"
    )

    object <- readYAMLSampleMetrics("project-summary.yaml")
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )

    # Check for proper handling of metrics with mismatched number of values
    object <- readYAMLSampleMetrics("project-summary-metrics-mismatch.yaml")
    expected[["sequenceLength"]] <- "numeric"
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )
})
