context("Import")



# readDataVersions =============================================================
test_that("readDataVersions", {
    x <- readDataVersions("data_versions.csv")
    expect_is(x, "DataFrame")
    expect_identical(
        object = colnames(x),
        expected = c("genome", "resource", "version")
    )
})

# Allow missing file, since bcbio doesn't always generate this.
test_that("readDataVersions : Missing file", {
    expect_message(
        object = readDataVersions("XXX.csv"),
        regexp = "Data versions are missing."
    )
    expect_identical(
        object = readDataVersions("XXX.csv"),
        expected = DataFrame()
    )
})



# readProgramVersions ==========================================================
test_that("readProgramVersions", {
    versions <- readProgramVersions("programs.txt")
    expect_is(versions, "DataFrame")
    expect_identical(
        object = colnames(versions),
        expected = c("program", "version")
    )
})

# Allow missing file, since bcbio doesn't always generate this.
test_that("readProgramVersions : Missing file", {
    expect_message(
        object = readProgramVersions("XXX.csv"),
        regexp = "Program versions are missing"
    )
    expect_identical(
        object = readProgramVersions("XXX.txt"),
        expected = DataFrame()
    )
})



# readSampleData =======================================================
test_that("readSampleData : Demultiplexed", {
    file <- "demultiplexed.csv"

    # Check DataFrame return.
    expect_identical(
        object = readSampleData(file),
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
        object = readSampleData("demultiplexed_invalid_missing.csv"),
        regexp = "isSampleData"
    )

    # Duplicated description.
    expect_error(
        object = readSampleData("demultiplexed_invalid_duplicated.csv"),
        regexp = "isSubset"
    )
})

test_that("readSampleData : Multiplexed", {
    file <- "multiplexed_indrops.csv"

    # Note that we're expecting this to sort by the rownames (`description`),
    # and not by the `sampleName` column.
    expect_identical(
        object = readSampleData(file),
        expected = DataFrame(
            sampleName = factor(c(
                "sample1_1",
                "sample2_1",
                "sample3_1",
                "sample4_1",
                "sample1_2",
                "sample2_2",
                "sample3_2",
                "sample4_2"
            )),
            fileName = factor(c(
                rep("indrops1_R1.fastq.gz", times = 4L),
                rep("indrops2_R1.fastq.gz", times = 4L)
            )),
            # Valid rownames (sampleIDs) are generated from this column.
            # Note that we're sorting the sample metadata by this column.
            description = factor(c(
                "indrops1-ATAGAGAG",
                "indrops1-AGAGGATA",
                "indrops1-CTCCTTAC",
                "indrops1-TATGCAGT",
                "indrops2-ATAGAGAG",
                "indrops2-AGAGGATA",
                "indrops2-CTCCTTAC",
                "indrops2-TATGCAGT"
            )),
            index = factor(rep(seq_len(4L), times = 2L)),
            sequence = factor(
                rep(c(
                    "CTCTCTAT",
                    "TATCCTCT",
                    "GTAAGGAG",
                    "ACTGCATA"
                ), times = 2L)
            ),
            aggregate = factor(
                paste0("sample", rep(seq_len(4L), times = 2L)),
                levels = paste0("sample", seq_len(4L))
            ),
            genotype = factor(
                rep(c("wildtype", "knockout"), times = 4L),
                # Note that the order should be alphabetical here.
                levels = c("knockout", "wildtype")
            ),
            revcomp = factor(
                rep(c(
                    "ATAGAGAG",
                    "AGAGGATA",
                    "CTCCTTAC",
                    "TATGCAGT"
                ), times = 2L)
            ),
            row.names = c(
                "indrops1_ATAGAGAG",
                "indrops1_AGAGGATA",
                "indrops1_CTCCTTAC",
                "indrops1_TATGCAGT",
                "indrops2_ATAGAGAG",
                "indrops2_AGAGGATA",
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
            "indrops1_L001_ATAGAGAG",
            "indrops1_L002_ATAGAGAG",
            "indrops1_L003_ATAGAGAG",
            "indrops1_L004_ATAGAGAG",
            "indrops2_L001_ATAGAGAG",
            "indrops2_L002_ATAGAGAG",
            "indrops2_L003_ATAGAGAG",
            "indrops2_L004_ATAGAGAG",
            "indrops1_L001_AGAGGATA",
            "indrops1_L002_AGAGGATA",
            "indrops1_L003_AGAGGATA",
            "indrops1_L004_AGAGGATA",
            "indrops2_L001_AGAGGATA",
            "indrops2_L002_AGAGGATA",
            "indrops2_L003_AGAGGATA",
            "indrops2_L004_AGAGGATA",
            "indrops1_L001_CTCCTTAC",
            "indrops1_L002_CTCCTTAC",
            "indrops1_L003_CTCCTTAC",
            "indrops1_L004_CTCCTTAC",
            "indrops2_L001_CTCCTTAC",
            "indrops2_L002_CTCCTTAC",
            "indrops2_L003_CTCCTTAC",
            "indrops2_L004_CTCCTTAC",
            "indrops1_L001_TATGCAGT",
            "indrops1_L002_TATGCAGT",
            "indrops1_L003_TATGCAGT",
            "indrops1_L004_TATGCAGT",
            "indrops2_L001_TATGCAGT",
            "indrops2_L002_TATGCAGT",
            "indrops2_L003_TATGCAGT",
            "indrops2_L004_TATGCAGT"
        )
    )
})

test_that("readSampleData : Multiplexed : Invalid", {
    # Required column check failure.
    expect_error(
        object = readSampleData("multiplexed_invalid_missing.csv"),
        expected = paste(
            "is_subset :",
            "The element 'index' in required is not in",
            "colnames\\(data\\)."
        )
    )

    # Duplicate rows in `sampleName` column.
    expect_error(
        object = readSampleData("multiplexed_invalid_duplicated.csv"),
        regexp = "hasNoDuplicates"
    )

    # Legacy bcbio `samplename` column.
    expect_error(
        object = readSampleData("demultiplexed_invalid_legacy_samplename.csv"),
        regexp = "Invalid columns: samplename"
    )

    # sampleID defined by user.
    expect_error(
        object = readSampleData("demultiplexed_invalid_sample_id.csv"),
        regexp = "Invalid columns: sampleID"
    )

    # Missing file.
    expect_error(
        object = readSampleData("XXX.csv"),
        regexp = "isAFile"
    )
})



# readTx2Gene ==================================================================
test_that("readTx2Gene", {
    object <- readTx2Gene(
        file = "tx2gene.csv",
        organism = "Mus musculus",
        genomeBuild = "GRCm38",
        ensemblRelease = 90L
    )
    expect_is(object, "Tx2Gene")
    expect_identical(
        object = colnames(object),
        expected = c("transcriptID", "geneID")
    )
})
