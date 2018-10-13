context("Import/Export")



# copyToDropbox ================================================================
if (file.exists("token.rds")) {
    files <- c("demultiplexed.csv", "multiplexed.csv")
    dropboxDir <- file.path("bcbioBase_examples", "copyToDropbox")
    test_that("copyToDropbox : RDS token enabled", {
        object <- copyToDropbox(
            files = files,
            dir = dropboxDir,
            rdsToken = "token.rds"
        )
        expect_is(object, "list")
        expect_identical(
            object = lapply(object[[1L]], class),
            expected = list(
                ".tag" = "character",
                url = "character",
                id = "character",
                name = "character",
                "path_lower" = "character",
                "link_permissions" = "list",
                "preview_type" = "character",
                "client_modified" = "character",
                "server_modified" = "character",
                rev = "character",
                size = "integer"
            )
        )
    })

    test_that("copyToDropbox : Shared Dropbox directory", {
        expect_warning(
            object = copyToDropbox(
                files = files,
                dir = paste0(dropboxDir, "_shared"),
                rdsToken = "token.rds"
            ),
            regexp = "rdrop2 currently isn't working well with shared"
        )
        # Don't clean up directory, because we won't be able to check if shared.
    })

    test_that("copyToDropbox : Invalid parameters", {
        expect_error(
            object = copyToDropbox(files = NULL, dir = "."),
            regexp = paste(
                "is2 :",
                "files is not in any of the classes 'character', 'list'."
            )
        )
        expect_error(
            object = copyToDropbox(files = "XXX.csv.gz", dir = "."),
            regexp = paste(
                "is_existing_file :",
                "Some or all of the files specified by files do not exist."
            )
        )
        expect_error(
            object = copyToDropbox(files = files, dir = NULL),
            regexp = paste(
                "is_a_string :",
                "dir is not of class 'character'; it has class 'NULL'"
            )
        )
        expect_error(
            object = copyToDropbox(
                files = files,
                dir = dropboxDir, rdsToken = "XXX.rds"
            ),
            regexp = paste(
                "is_existing_file :",
                "Some or all of the files specified by rdsToken do not exist."
            )
        )
    })
}



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
})

test_that("readSampleData : Demultiplexed : Invalid", {
    # Required column check failure.
    expect_error(
        object = readSampleData("demultiplexed_invalid_missing.csv"),
        regexp = "is_subset : The element 'description'"
    )

    # Duplicated description.
    expect_error(
        object = readSampleData("demultiplexed_invalid_duplicated.csv"),
        regexp = "is_subset : The elements 'sampleName', 'index'"
    )
})

test_that("readSampleData : Multiplexed : bcbioSingleCell", {
    file <- "multiplexed_indrops.csv"

    # Note that we're expecting this to sort by the rownames (`description`),
    # and not by the `sampleName` column.
    expect_identical(
        object = readSampleData(file),
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
            genotype = factor(
                c(
                    "knockout",
                    "wildtype",
                    "wildtype",
                    "knockout",
                    "knockout",
                    "wildtype",
                    "wildtype",
                    "knockout"
                ),
                levels = c("knockout", "wildtype")  # Note order here.
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
})

test_that("readSampleData : Multiplexed : Cell Ranger", {
    expect_identical(
        object = readSampleData("multiplexed_cellranger.csv"),
        expected = DataFrame(
            sampleName = factor(c("proximal", "distal")),
            fileName = factor("aggr_R1.fastq.gz"),
            description = factor(c("aggr-1", "aggr-2")),
            index = factor(c("1", "2")),
            row.names = factor(c("aggr_1", "aggr_2"))
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
        regexp = paste(
            "has_no_duplicates :",
            "data\\[\\[\"sampleName\"\\]\\] has duplicates at positions 2, 4."
        )
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
        regexp = "is_existing_file :"
    )
})



# readTx2Gene ==================================================================
test_that("readTx2Gene", {
    object <- readTx2Gene("tx2gene.csv")
    expect_is(object, "Tx2Gene")
    expect_identical(
        object = colnames(object),
        expected = c("transcriptID", "geneID")
    )
})
