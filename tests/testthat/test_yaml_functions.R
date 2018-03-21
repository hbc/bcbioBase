context("YAML Functions")

yaml <- readYAML("project-summary.yaml")



# sampleYAML ===================================================================
test_that("sampleYAML", {
    expect_identical(
        sampleYAML(yaml, "metadata"),
        tibble(
            "description" = c(
                "group1_1",
                "group1_2",
                "group2_1",
                "group2_2"
            ),
            "group" = c(
                "ctrl",
                "ctrl",
                "ko",
                "ko"
            )
        )
    )
})



# sampleYAMLMetadata ===========================================================
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



# sampleYAMLMetrics ============================================================
test_that("sampleYAMLMetrics", {
    classChecks <- list(
        "averageInsertSize" = "numeric",
        "duplicates" = "numeric",
        "duplicationRateOfMapped" = "numeric",
        "exonicRate" = "numeric",
        "intergenicRate" = "numeric",
        "intronicRate" = "numeric",
        "mappedPairedReads" = "numeric",
        "mappedReads" = "numeric",
        "qualityFormat" = "factor",
        "rrna" = "numeric",
        "rrnaRate" = "numeric",
        "sequenceLength" = "factor",
        "sequencesFlaggedAsPoorQuality" = "numeric",
        "totalReads" = "numeric",
        "x5x3Bias" = "numeric",  # 5'3 now sanitized to 5x3 in camel
        "xGC" = "numeric"
    )

    x <- sampleYAMLMetrics(yaml)
    expect_identical(lapply(x, class), classChecks)

    # Check for proper handling of metrics with mismatched number of values
    yaml <- readYAML("project-summary-metrics-mismatch.yaml")
    x <- sampleYAMLMetrics(yaml)
    classChecks[["sequenceLength"]] <- "numeric"
    expect_identical(lapply(x, class), classChecks)
})

test_that("sampleYAMLMetrics : Fast mode", {
    fastmode <- "Fast mode detected: No sample metrics were calculated"

    # Subset to only include the first sample
    single <- yaml
    single[["samples"]] <- single[["samples"]][1L]

    # NULL metrics
    nullmetrics <- single
    nullmetrics[["samples"]][[1L]][["summary"]][["metrics"]] <- NULL
    expect_warning(
        sampleYAMLMetrics(nullmetrics),
        fastmode
    )
    expect_identical(
        suppressWarnings(sampleYAMLMetrics(nullmetrics)),
        NULL
    )

    # Empty metrics
    emptymetrics <- single
    emptymetrics[["samples"]][[1L]][["summary"]][["metrics"]] <- list()
    expect_warning(
        sampleYAMLMetrics(emptymetrics),
        fastmode
    )
    expect_identical(
        suppressWarnings(sampleYAMLMetrics(emptymetrics)),
        NULL
    )
})



# General ======================================================================
test_that("Invalid parameters", {
    expect_error(
        sampleYAML(yaml = list(), keys = "metadata"),
        "is_non_empty : yaml has length 0."
    )
    # Primary key failure
    expect_error(
        sampleYAML(yaml, keys = "XXX"),
        paste(
            "is_subset :",
            "The element 'XXX' in keys\\[\\[1L\\]\\] is not in",
            "names\\(yaml\\[\\[1L\\]\\]\\)."
        )
    )
    # Secondary key failure
    expect_error(
        sampleYAML(yaml, keys = c("summary", "XXX")),
        paste(
            "is_subset :",
            "The element 'XXX' in keys\\[\\[2L\\]\\] is not in",
            "names\\(yaml\\[\\[1L\\]\\]\\[\\[keys\\[\\[1L\\]\\]\\]\\]\\)."
        )
    )
})
