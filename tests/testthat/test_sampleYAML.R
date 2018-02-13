context("sampleYAML")

yaml <- readYAML(
    file.path(
        "http://bcbiobase.seq.cloud",
        "bcbio",
        "project-summary.yaml"),
    quiet = TRUE)

test_that("sampleYAML", {
    expect_identical(
        sampleYAML(yaml, "metadata"),
        data.frame(
            "group" = c(
                "ctrl",
                "ctrl",
                "ko",
                "ko"),
            "description" = c(
                "group1_1",
                "group1_2",
                "group2_1",
                "group2_2"),
            row.names = c(
                "group1_1",
                "group1_2",
                "group2_1",
                "group2_2"),
            stringsAsFactors = FALSE
        )
    )
})

test_that("Invalid YAML input", {
    expect_error(
        sampleYAML(yaml = list(), keys = "metadata"),
        "is_non_empty : samples has length 0."
    )
    # Primary key failure
    expect_error(
        sampleYAML(yaml, keys = "XXX"),
        paste(
            "is_subset :",
            "The element 'XXX' in keys\\[\\[1L\\]\\] is not in",
            "names\\(samples\\[\\[1L\\]\\]\\)."
        )
    )
    # Secondary key failure
    expect_error(
        sampleYAML(yaml, keys = c("summary", "XXX")),
        paste(
            "is_subset :",
            "The element 'XXX' in keys\\[\\[2L\\]\\] is not in",
            "names\\(samples\\[\\[1L\\]\\]\\[\\[keys\\[\\[1L\\]\\]\\]\\]\\)."
        )
    )
})

test_that("sampleYAMLMetadata", {
    samples <- c("group1_1", "group1_2", "group2_1", "group2_2")
    expect_identical(
        sampleYAMLMetadata(yaml),
        data.frame(
            sampleID = factor(samples, levels = samples),
            sampleName = factor(samples, levels = samples),
            description = factor(samples, levels = samples),
            group = factor(
                c("ctrl", "ctrl", "ko", "ko"),
                levels = c("ctrl", "ko")
            ),
            row.names = samples
        )
    )
})

test_that("sampleYAMLMetrics", {
    metrics <- sampleYAMLMetrics(yaml)
    expect_identical(
        vapply(
            X = metrics,
            FUN = class,
            FUN.VALUE = "character"),
        c(
            sampleID = "factor",
            sampleName = "factor",
            description = "factor",
            xGC = "numeric",
            x5x3Bias = "numeric",  # 5'3 now sanitized to 5x3 in camel
            averageInsertSize = "numeric",
            duplicates = "numeric",
            duplicationRateOfMapped = "numeric",
            exonicRate = "numeric",
            intergenicRate = "numeric",
            intronicRate = "numeric",
            mappedPairedReads = "numeric",
            mappedReads = "numeric",
            name = "factor",
            qualityFormat = "factor",
            sequenceLength = "factor",
            sequencesFlaggedAsPoorQuality = "numeric",
            totalReads = "numeric",
            rrna = "numeric",
            rrnaRate = "numeric")
    )

    # Check for proper handling of metrics with mismatched number of values
    yaml2 <- readYAML(
        file.path(
            "http://bcbiobase.seq.cloud",
            "bcbio",
            "project-summary-metrics-mismatch.yaml"),
        quiet = TRUE)
    metrics2 <- sampleYAMLMetrics(yaml2)
    expect_identical(
        vapply(metrics2, class, FUN.VALUE = "character"),
        c(
            sampleID = "factor",
            sampleName = "factor",
            description = "factor",
            xGC = "numeric",
            x5x3Bias = "numeric",
            averageInsertSize = "numeric",
            duplicates = "numeric",
            duplicationRateOfMapped = "numeric",
            exonicRate = "numeric",
            intergenicRate = "numeric",
            intronicRate = "numeric",
            mappedPairedReads = "numeric",
            mappedReads = "numeric",
            name = "factor",
            qualityFormat = "factor",
            sequenceLength = "numeric",  # factor in the main example
            sequencesFlaggedAsPoorQuality = "numeric",
            totalReads = "numeric",
            rrna = "numeric",
            rrnaRate = "numeric")
    )
})

test_that("No sample metrics", {
    nometrics <- yaml
    # Subset to only include the first sample
    nometrics[["samples"]] <- nometrics[["samples"]][[1]]
    nometrics[["samples"]][["summary"]][["metrics"]] <- NULL
    expect_error(
        sampleYAMLMetrics(nometrics),
        paste(
            "is_subset :",
            "The element 'summary' in keys\\[\\[1L\\]\\] is not in",
            "names\\(samples\\[\\[1L\\]\\]\\)."
        )
    )
})
