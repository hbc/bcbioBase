context("Data Functions")



# convertGenesToSymbols ========================================================
test_that("convertGenesToSymbols", {
    x <- convertGenesToSymbols(rse_bcb)
    expect_identical(
        head(rownames(x)),
        c("Cox5a", "Comt", "Dazap2", "Rpl13", "Calm1", "Ddt")
    )
})

test_that("convertGenesToSymbols : unmodified return", {
    expect_warning(
        convertGenesToSymbols(rse_dds),
        "Object does not contain gene-to-symbol mappings"
    )
    x <- suppressWarnings(convertGenesToSymbols(rse_dds))
    expect_identical(rownames(x), rownames(rse_dds))
})



# counts =======================================================================
test_that("counts", {
    x <- counts(rse_dds)
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
    expect_warning(
        gene2symbol(rse_dds),
        "Object does not contain gene-to-symbol mappings"
    )
    expect_identical(
        suppressWarnings(gene2symbol(rse_dds)),
        NULL
    )
})



# flatFiles ====================================================================
test_that("flatFiles : SummarizedExperiment", {
    x <- flatFiles(rse_dds)
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
    expect_identical(
        interestingGroups(rse_dds),
        NULL
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
        "The interesting groups \"XXX\" are not defined"
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
        sampleData(
            rse_bcb,
            clean = TRUE,
            return = x
        )
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
            c(6L, 7L),
            c(6L, 7L),
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
    sampleData(x)[["test"]] <- as.factor(seq_len(ncol(x)))
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



# selectSamples ================================================================
test_that("selectSamples : SummarizedExperiment", {
    x <- selectSamples(rse_dds, condition = "A")
    expect_identical(dim(x), c(1000L, 6L))
    expect_identical(colnames(x), paste0("sample", seq(6L)))
})



# uniteInterestingGroups =======================================================
test_that("uniteInterestingGroups : Single interesting group", {
    x <- uniteInterestingGroups(
        object = datasets::mtcars,
        interestingGroups = c("vs", "am", "gear")
    )
    expect_identical(
        levels(x[["interestingGroups"]]),
        c("0:0:3", "0:1:4", "0:1:5", "1:0:3", "1:0:4", "1:1:4", "1:1:5")
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
        object = datasets::mtcars,
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
        uniteInterestingGroups(
            object = datasets::mtcars,
            interestingGroups = c("XXX", "YYY")
        ),
        "is_subset : The elements 'XXX', 'YYY' in interestingGroups"
    )
})
