yaml <- import(file.path("cache", "summary.yaml"))

test_that("getGtfFileFromYaml", {
    expect_message(
        object = getGtfFileFromYaml(yaml),
        regexp = "ref-transcripts.gtf"
    )
    expect_null(getGtfFileFromYaml(yaml))
})

test_that("getSampleDataFromYaml", {
    object <- getSampleDataFromYaml(yaml)
    samples <- c("group1_1", "group1_2", "group2_1", "group2_2")
    expected <- DataFrame(
        "sampleName" = factor(samples),
        "description" = factor(samples),
        "genomeBuild" = factor("mm10"),
        "group" = factor(c("ctrl", "ctrl", "ko", "ko")),
        "samRef" = factor(paste(
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
    expect_identical(object, expected)
})

## Testing against Kayleigh's nested example here.
test_that("Nested metadata", {
    ## Expecting warnings about integer range here.
    suppressWarnings({
        yaml <- import(file.path("cache", "summary-nested-metadata.yaml"))
    })
    object <- getSampleDataFromYaml(yaml)
    expect_s4_class(object, "DataFrame")
    expect_identical(dim(object), c(218L, 29L))
    expect_identical(
        object = rownames(object)[[1L]],
        expected = "SRR1022936"
    )
})

expected <- list(
    "averageInsertSize" = "integer",
    "duplicates" = "integer",
    "duplicationRateOfMapped" = "numeric",
    "exonicRate" = "numeric",
    "intergenicRate" = "numeric",
    "intronicRate" = "numeric",
    "mappedPairedReads" = "integer",
    "mappedReads" = "integer",
    "percentGc" = "integer",
    "qualityFormat" = "factor",
    "rrna" = "numeric",
    "rrnaRate" = "numeric",
    "sequenceLength" = "factor",
    "sequencesFlaggedAsPoorQuality" = "integer",
    "totalReads" = "integer",
    "x5x3Bias" = "numeric"
)

test_that("getMetricsFromYaml", {
    object <- getMetricsFromYaml(yaml)
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )
})

## Check for proper handling of metrics with mismatched number of values.
test_that("Mismatched values", {
    file <- file.path("cache", "summary-invalid-metrics-mismatch.yaml")
    yaml <- import(file)
    object <- getMetricsFromYaml(yaml)
    expected[["sequenceLength"]] <- "integer"
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )
})
