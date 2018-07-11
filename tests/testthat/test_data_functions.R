context("Data Functions")



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
    # S4 coercion to list method support
    y <- as(rse_dds, "list")
    expect_identical(x, y)
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



# minimalSampleData ============================================================
test_that("minimalSampleData", {
    expect_identical(
        minimalSampleData(c("sample 1", "sample 2")),
        data.frame(
            sampleName = c("sample 1", "sample 2"),
            description = c("sample 1", "sample 2"),
            row.names = c("sample_1", "sample_2"),
            stringsAsFactors = TRUE
        )
    )
})



# sampleData ===================================================================
# Check output of `return` parameter
return <- methodFormals(
    f = "sampleData",
    signature = "SummarizedExperiment"
) %>%
    .[["return"]] %>%
    as.character() %>%
    .[-1L]

test_that("sampleData: Verbose mode", {
    list <- lapply(return, function(x) {
        sampleData(
            rse_bcb,
            clean = FALSE,
            return = x
        )
    })

    # Check returns
    expect_identical(
        lapply(list, class),
        list(
            structure("DataFrame", package = "S4Vectors"),
            "data.frame",
            "knitr_kable"
        )
    )

    # Check dimnames
    expected <- list(
        colnames(rse_bcb),
        c(colnames(colData(rse_bcb)), "interestingGroups")
    )
    expect_identical(
        lapply(list, dimnames),
        list(
            expected,
            expected,
            NULL
        )
    )
})

test_that("sampleData : Clean mode", {
    list <- lapply(return, function(x) {
        sampleData(
            rse_bcb,
            clean = TRUE,
            return = x
        )
    })

    # Check returns
    expect_identical(
        lapply(list, class),
        list(
            structure("DataFrame", package = "S4Vectors"),
            "data.frame",
            "knitr_kable"
        )
    )

    # Check dimnames
    expected <- list(
        colnames(rse_bcb),
        c(
            "sampleName",
            "day",
            "replicate",
            "strain",
            "tissue",
            "treatment"
        )
    )
    expect_identical(
        lapply(list, dimnames),
        list(
            expected,
            expected,
            NULL
        )
    )

    # Ensure all columns are factor
    invisible(lapply(list[[1L]], function(x) {
        expect_is(x, "factor")
    }))

    # Interesting groups
    x <- sampleData(rse_bcb, clean = FALSE, interestingGroups = NULL)
    expect_identical(
        x[["interestingGruops"]],
        NULL
    )
    x <- sampleData(rse_bcb, clean = FALSE, interestingGroups = "day")
    expect_identical(
        levels(x[["interestingGroups"]]),
        c("0", "7")
    )
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



# sampleNames ==================================================================
test_that("sampleNames", {
    x <- sampleNames(rse_bcb)
    expect_identical(
        x[seq_len(2L)],
        c(
            control_rep1 = "control_rep1",
            control_rep2 = "control_rep2"
        )
    )

    x <- sampleNames(rse_dds)
    expect_identical(
        x[seq_len(2L)],
        c(
            sample1 = "sample1",
            sample10 = "sample10"
        )
    )
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
