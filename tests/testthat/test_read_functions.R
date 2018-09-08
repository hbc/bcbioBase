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

test_that("readDataVersions : Missing file", {
    expect_message(
        object = readDataVersions("XXX.csv"),
        regexp = "Data versions are missing"
    )
    expect_identical(
        object = readDataVersions("XXX.csv"),
        expected = tibble::tibble()
    )
})



# readLog ==================================================================
test_that("readLog", {
    log <- readLog("bcbio_nextgen.log")
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

test_that("readProgramVersions : Missing file", {
    expect_message(
        object = readProgramVersions("XXX.csv"),
        regexp = "Program versions are missing"
    )
    expect_identical(
        object = readProgramVersions("XXX.txt"),
        expected = tibble::tibble()
    )
})



# readSampleData =======================================================
test_that("readSampleData : Demultiplexed FASTQ", {
    file <- "demultiplexed.csv"
    object <- readSampleData(file)

    # Check for message.
    expect_message(
        object = readSampleData(file),
        regexp = "Demultiplexed samples detected"
    )

    # Check DataFrame return.
    expect_identical(
        object = object,
        expected = DataFrame(
            sampleName = factor(paste0("sample", seq_len(4L))),
            fileName = factor(paste0("sample", seq_len(4L), "_R1.fastq.gz")),
            description = factor(paste0("sample", seq_len(4L))),
            genotype = factor(rep(c("wildtype", "knockout"), times = 2L)),
            row.names = paste0("sample", seq_len(4L))
        )
    )

    # Lane-split technical replicate support.
    object <- readSampleData(file, lanes = 4L)
    expect_true("lane" %in% colnames(object))
    expect_identical(
        object = rownames(object)[1L:8L],
        expected = c(
            paste0("sample1_L00", seq_len(4L)),
            paste0("sample2_L00", seq_len(4L))
        )
    )

    # Required column check failure.
    expect_error(
        object = readSampleData("demultiplexed_missing_columns.csv"),
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

    # Check for message.
    expect_message(
        object = readSampleData(file),
        regexp = "Multiplexed samples detected"
    )

    # Note that we're expecting this to sort by the rownames, and not by the
    # `sampleName` column.
    expect_identical(
        object = object,
        expected = DataFrame(
            sampleName = factor(
                c(
                    "sample2_1",
                    "sample1_1",
                    "sample3_1",
                    "sample4_1",
                    "sample2_2",
                    "sample1_2",
                    "sample3_2",
                    "sample4_2"
                ),
                levels = c(
                    "sample1_1",
                    "sample1_2",
                    "sample2_1",
                    "sample2_2",
                    "sample3_1",
                    "sample3_2",
                    "sample4_1",
                    "sample4_2"
                )
            ),
            fileName = factor(c(
                rep("indrops1_R1.fastq.gz", times = 4L),
                rep("indrops2_R1.fastq.gz", times = 4L)
            )),
            # Valid rownames (sampleIDs) are generated from this column.
            # Note that we're sorting the sample metadata by this column.
            description = factor(c(
                "indrops1-AGAGGATA",
                "indrops1-ATAGAGAG",
                "indrops1-CTCCTTAC",
                "indrops1-TATGCAGT",
                "indrops2-AGAGGATA",
                "indrops2-ATAGAGAG",
                "indrops2-CTCCTTAC",
                "indrops2-TATGCAGT"
            )),
            index = factor(
                c(2L, 1L, 3L, 4L, 2L, 1L, 3L, 4L),
                levels = seq_len(4L)
            ),
            sequence = factor(
                c(
                    "TATCCTCT",
                    "CTCTCTAT",
                    "GTAAGGAG",
                    "ACTGCATA",
                    "TATCCTCT",
                    "CTCTCTAT",
                    "GTAAGGAG",
                    "ACTGCATA"
                ),
                levels = c("ACTGCATA", "CTCTCTAT", "GTAAGGAG", "TATCCTCT")
            ),
            aggregate = factor(
                paste0("sample", c(2L, 1L, 3L, 4L, 2L, 1L, 3L, 4L)),
                levels = paste0("sample", seq_len(4L))
            ),
            revcomp = factor(
                c(
                    "AGAGGATA",
                    "ATAGAGAG",
                    "CTCCTTAC",
                    "TATGCAGT",
                    "AGAGGATA",
                    "ATAGAGAG",
                    "CTCCTTAC",
                    "TATGCAGT"
                ),
                levels = c("AGAGGATA", "ATAGAGAG", "CTCCTTAC", "TATGCAGT")
            ),
            row.names = c(
                "indrops1_AGAGGATA",
                "indrops1_ATAGAGAG",
                "indrops1_CTCCTTAC",
                "indrops1_TATGCAGT",
                "indrops2_AGAGGATA",
                "indrops2_ATAGAGAG",
                "indrops2_CTCCTTAC",
                "indrops2_TATGCAGT"
            )
        )
    )

    # Lane-split technical replicate support.
    object <- readSampleData(file, lanes = 4L)
    expect_identical(
        object = rownames(object),
        expected = c(
            "indrops1_L001_AGAGGATA",
            "indrops1_L001_ATAGAGAG",
            "indrops1_L001_CTCCTTAC",
            "indrops1_L001_TATGCAGT",
            "indrops1_L002_AGAGGATA",
            "indrops1_L002_ATAGAGAG",
            "indrops1_L002_CTCCTTAC",
            "indrops1_L002_TATGCAGT",
            "indrops1_L003_AGAGGATA",
            "indrops1_L003_ATAGAGAG",
            "indrops1_L003_CTCCTTAC",
            "indrops1_L003_TATGCAGT",
            "indrops1_L004_AGAGGATA",
            "indrops1_L004_ATAGAGAG",
            "indrops1_L004_CTCCTTAC",
            "indrops1_L004_TATGCAGT",
            "indrops2_L001_AGAGGATA",
            "indrops2_L001_ATAGAGAG",
            "indrops2_L001_CTCCTTAC",
            "indrops2_L001_TATGCAGT",
            "indrops2_L002_AGAGGATA",
            "indrops2_L002_ATAGAGAG",
            "indrops2_L002_CTCCTTAC",
            "indrops2_L002_TATGCAGT",
            "indrops2_L003_AGAGGATA",
            "indrops2_L003_ATAGAGAG",
            "indrops2_L003_CTCCTTAC",
            "indrops2_L003_TATGCAGT",
            "indrops2_L004_AGAGGATA",
            "indrops2_L004_ATAGAGAG",
            "indrops2_L004_CTCCTTAC",
            "indrops2_L004_TATGCAGT"
        )
    )

    # Required column check failure.
    expect_error(
        object = readSampleData("multiplexed_missing_columns.csv"),
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
    object <- readSampleData("multiplexed_cellranger.csv")
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
        object = readSampleData("multiplexed_sampleID_column_defined.csv"),
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
    object <- readYAMLSampleData("project_summary.yaml")
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
        readYAMLSampleData("project_summary_nested_metadata.yaml")
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

    object <- readYAMLSampleMetrics("project_summary.yaml")
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )

    # Check for proper handling of metrics with mismatched number of values
    object <- readYAMLSampleMetrics("project_summary_metrics_mismatch.yaml")
    expected[["sequenceLength"]] <- "numeric"
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )
})
