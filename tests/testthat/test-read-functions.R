context("Read Functions")



# readDataVersions =============================================================
test_that("readDataVersions", {
    versions <- readDataVersions(file.path("cache", "data-versions.csv"))
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
    log <- readLog(file.path("cache", "bcbio-nextgen.log"))
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
        "XXX.log"
    )
})



# readProgramVersions ==========================================================
test_that("readProgramVersions", {
    versions <- readProgramVersions(file.path("cache", "programs.txt"))
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
    file <- file.path("cache", "demultiplexed.csv")
    x <- readSampleData(file)

    # Check that names are sanitized correctly
    expect_identical(
        rownames(x),
        c("sample1", "sample2", "sample3", "sample4")
    )

    # Check that column names get set correctly
    expect_identical(
        colnames(x),
        c("sampleName", "fileName", "description", "genotype")
    )

    # Lane-split technical replicate support
    x <- readSampleData(file, lanes = 4L)
    expect_identical(
        colnames(x),
        c("sampleName", "description", "lane", "fileName", "genotype")
    )
    expect_identical(
        rownames(x)[1L:8L],
        c(
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

    # Required column check failure
    expect_error(
        readSampleData(file.path("cache", "demultiplexed-missing-cols.csv")),
        "description"
    )

    # Duplicated description
    expect_error(
        readSampleData(
            file.path("cache", "demultiplexed-duplicated-description.csv")
        ),
        "'sampleName', 'index'"
    )
})

test_that("readSampleData : Multiplexed FASTQ", {
    file <- file.path("cache", "multiplexed.csv")

    x <- readSampleData(file)
    expect_identical(
        rownames(x),
        c(
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

    # Lane-split technical replicate support
    x <- readSampleData(file, lanes = 4L)
    expect_identical(
        head(rownames(x), n = 4L),
        c(
            "indrops1_L001_AGAGGATA",
            "indrops1_L001_ATAGAGAG",
            "indrops1_L001_CTCCTTAC",
            "indrops1_L001_TATGCAGT"
        )
    )

    # Required column check failure
    expect_error(
        readSampleData(file.path("cache", "multiplexed-missing-cols.csv")),
        "index"
    )

    # Duplicate rows in `sampleName` column
    expect_error(
        readSampleData(
            file.path("cache", "multiplexed-duplicated-sampleName.csv")
        ),
        "sampleName"
    )
})

test_that("readSampleData : Multiplexed CellRanger data", {
    x <- readSampleData(file.path("cache", "cellranger-metadata.csv"))
    y <- data.frame(
        sampleName = c("proximal", "distal"),
        fileName = "aggregation.fastq.gz",
        description = c("aggregation-1", "aggregation-2"),
        index = c("1", "2"),
        row.names = c("aggregation_1", "aggregation_2"),
        stringsAsFactors = TRUE
    )
    expect_identical(x, y)
})

test_that("readSampleData : Legacy bcbio samplename column", {
    file <- file.path("cache", "bcbio-legacy-samplename.csv")
    # Warn on `samplename`
    expect_warning(
        readSampleData(file),
        "Invalid metadata columns detected."
    )
    expect_identical(
        suppressWarnings(readSampleData(file)),
        data.frame(
            sampleName = "sample1",
            description = "sample1",
            fileName = "sample1.fastq.gz",
            row.names = "sample1",
            stringsAsFactors = TRUE
        )
    )
})

test_that("readSampleData : sampleID defined by user", {
    expect_warning(
        readSampleData(file.path("cache", "sampleID-column-defined.csv")),
        "Invalid metadata columns detected."
    )
})

test_that("readSampleData : Missing file", {
    # Always stop on missing
    expect_error(
        readSampleData("XXX.csv"),
        "XXX.csv"
    )
})



# readTx2gene ==================================================================
test_that("readTx2gene", {
    x <- readTx2gene(file.path("cache", "tx2gene.csv"))
    expect_is(x, "data.frame")
    expect_identical(
        colnames(x),
        c("transcriptID", "geneID")
    )
})



# readYAMLSampleData ===========================================================
test_that("readYAMLSampleData", {
    x <- readYAMLSampleData(file.path("cache", "project-summary.yaml"))
    samples <- c("group1_1", "group1_2", "group2_1", "group2_2")
    expect_identical(
        x,
        data.frame(
            sampleName = samples,
            description = samples,
            genomeBuild = "mm10",
            group = c("ctrl", "ctrl", "ko", "ko"),
            samRef = paste(
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
            ),
            row.names = samples,
            stringsAsFactors = TRUE
        )
    )
})

test_that("readYAMLSampleData : nested metadata", {
    # Testing against Kayleigh's example
    x <- suppressWarnings(
        readYAMLSampleData(
            file.path("cache", "project-summary-nested-metadata.yaml")
        )
    )
    expect_is(x, "data.frame")
})



# readYAMLSampleMetrics ========================================================
test_that("readYAMLSampleMetrics", {
    y <- list(
        averageInsertSize = "numeric",
        duplicates = "numeric",
        duplicationRateOfMapped = "numeric",
        exonicRate = "numeric",
        intergenicRate = "numeric",
        intronicRate = "numeric",
        mappedPairedReads = "numeric",
        mappedReads = "numeric",
        qualityFormat = "factor",
        rrna = "numeric",
        rrnaRate = "numeric",
        sequenceLength = "factor",
        sequencesFlaggedAsPoorQuality = "numeric",
        totalReads = "numeric",
        x5x3Bias = "numeric",
        xGC = "numeric"
    )

    x <- readYAMLSampleMetrics(file.path("cache", "project-summary.yaml"))
    expect_identical(lapply(x, class), y)

    # Check for proper handling of metrics with mismatched number of values
    x <- readYAMLSampleMetrics(
        file.path("cache", "project-summary-metrics-mismatch.yaml")
    )
    y[["sequenceLength"]] <- "numeric"
    expect_identical(lapply(x, class), y)
})
