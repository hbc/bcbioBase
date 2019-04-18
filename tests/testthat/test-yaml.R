yaml <- import(file.path("cache", "summary.yaml"))



context("getGTFFileFromYAML")

test_that("getGTFFileFromYAML", {
    expect_message(
        object = getGTFFileFromYAML(yaml),
        regexp = "ref-transcripts.gtf"
    )
    expect_null(getGTFFileFromYAML(yaml))
})



context("getSampleDataFromYAML")

test_that("getSampleDataFromYAML", {
    object <- getSampleDataFromYAML(yaml)
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

# Testing against Kayleigh's nested example here.
test_that("Nested metadata", {
    # Expecting warnings about integer range here.
    object <- suppressWarnings(
        getSampleDataFromYAML(
            yaml = import(file.path("cache", "summary-nested-metadata.yaml"))
        )
    )
    expect_is(object, "DataFrame")
})



context("getMetricsFromYAML")

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

test_that("getMetricsFromYAML", {
    object <- getMetricsFromYAML(yaml)
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )
})

# Check for proper handling of metrics with mismatched number of values.
test_that("Mismatched values", {
    yaml <- import(file.path("summary-invalid-metrics-mismatch.yaml"))
    object <- getMetricsFromYAML(yaml)
    expected[["sequenceLength"]] <- "numeric"
    expect_identical(
        object = lapply(object, class),
        expected = expected
    )
})
