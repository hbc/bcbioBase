context("Data Functions")

yaml <- readYAML("project-summary.yaml")



# convertGenesToSymbols ========================================================
test_that("convertGenesToSymbols", {
    x <- convertGenesToSymbols(rse_bcb)
    expect_identical(
        head(rownames(x)),
        c("Cox5a", "Comt", "Dazap2", "Rpl13", "Calm1", "Ddt")
    )
})

test_that("convertGenesToSymbols : unmodified return", {
    x <- convertGenesToSymbols(rse_dds)
    expect_identical(rownames(x), rownames(rse_dds))
})



# counts =======================================================================
test_that("counts", {
    x <- counts(rse_bcb)
    expect_is(x, "matrix")
})



# gene2symbol ==================================================================
test_that("gene2symbol", {
    x <- gene2symbol(rse_bcb)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), c("geneID", "geneName"))
    expect_true(tibble::has_rownames(x))
})

test_that("gene2symbol : NULL return", {
    expect_identical(
        gene2symbol(rse_dds),
        NULL
    )
})



# flatFiles ====================================================================
test_that("flatFiles : SummarizedExperiment", {
    x <- flatFiles(rse_bcb)
    expect_is(x, "list")
    expect_identical(
        names(x),
        c(
            "rowRanges",
            "colData",
            "assays",
            "NAMES",
            "elementMetadata",
            "metadata"
        )
    )
})



# interestingGroups ============================================================
test_that("interestingGroups : SummarizedExperiment", {
    expect_identical(
        interestingGroups(rse_bcb),
        "treatment"
    )
})

test_that("interestingGroups : Assignment method", {
    x <- rse_bcb
    interestingGroups(x) <- "sampleName"
    expect_identical(
        interestingGroups(x),
        "sampleName"
    )
    expect_error(
        interestingGroups(x) <- "XXX",
        "is_subset : The element 'XXX' in interestingGroups"
    )
})



# sampleData ===================================================================
test_that("sampleData : SummarizedExperiment", {
    # Check output of `return` parameter
    return <- methodFormals(
        f = "sampleData",
        signature = "SummarizedExperiment"
    ) %>%
        .[["return"]] %>%
        as.character() %>%
        .[-1L]
    list <- lapply(return, function(x) {
        sampleData(rse_bcb, return = x)
    })
    expect_identical(
        lapply(list, class),
        list(
            structure("DataFrame", package = "S4Vectors"),
            "data.frame",
            "knitr_kable"
        )
    )
    # Check dimensions
    expect_identical(
        lapply(list, dim),
        list(
            c(6L, 9L),
            c(6L, 9L),
            NULL
        )
    )
    # Check rownames
    expect_identical(
        lapply(list, rownames),
        list(
            colnames(rse_bcb),
            colnames(rse_bcb),
            NULL
        )
    )
    # Ensure all columns are factor
    invisible(lapply(list[[1L]], function(x) {
        expect_is(x, "factor")
    }))
})

test_that("sampleData : Assignment method", {
    x <- rse_bcb
    sampleData(x)[["test"]] <- seq_len(ncol(x))
    expect_is(sampleData(x)[["test"]], "factor")
})



# sampleDirs ===================================================================
test_that("sampleDirs", {
    uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
    expect_identical(
        sampleDirs(uploadDir),
        c(
            "sample_1" = file.path(uploadDir, "sample_1"),
            "sample_2" = file.path(uploadDir, "sample_2")
        )
    )
})



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

test_that("sampleYAML : Invalid parameters", {
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

test_that("sampleYAML : nested metadata", {
    # Using Kayleigh's bcbio example
    yaml <- suppressWarnings(readYAML("project-summary-nested-metadata.yaml"))
    x <- sampleYAMLMetadata(yaml)
    expect_is(x, "data.frame")
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
    single[["samples"]] <- head(single[["samples"]], 1L)

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



# selectSamples ================================================================
test_that("selectSamples : SummarizedExperiment", {
    x <- selectSamples(rse_bcb, treatment = "folic_acid")
    expect_identical(dim(x), c(500L, 3L))
    expect_identical(names(assays(x)), "counts")
})



# uniteInterestingGroups =======================================================
test_that("uniteInterestingGroups : Single interesting group", {
    x <- uniteInterestingGroups(mtcars, interestingGroups = "gear")
    expect_identical(
        x[["interestingGroups"]],
        as.factor(mtcars[["gear"]])
    )
})

test_that("uniteInterestingGroups : tidy (tibble) mode", {
    x <- uniteInterestingGroups(
        object = dplyr::starwars,
        interestingGroups = c("hair_color", "skin_color")
    )
    expect_is(x, "tbl_df")
    expect_is(x[["interestingGroups"]], "factor")
    expect_identical(
        x[["interestingGroups"]] %>%
            as.character() %>%
            head(2L),
        c("blond:fair", "NA:gold")
    )
})

test_that("uniteInterestingGroups : Two interesting groups", {
    x <- uniteInterestingGroups(
        mtcars,
        interestingGroups = c("gear", "carb")
    )
    expect_identical(
        head(x[["interestingGroups"]]),
        factor(
            c("4:4", "4:4", "4:1", "3:1", "3:2", "3:1"),
            levels = c(
                "3:1", "3:2", "3:3", "3:4",
                "4:1", "4:2", "4:4",
                "5:2", "5:4", "5:6", "5:8"
            )
        )
    )
})

test_that("uniteInterestingGroups : Missing groups", {
    expect_error(
        uniteInterestingGroups(mtcars, interestingGroups = c("XXX", "YYY")),
        paste(
            "is_subset :",
            "The elements 'XXX', 'YYY' in interestingGroups are not in",
            "colnames\\(x\\)."
        )
    )
})
