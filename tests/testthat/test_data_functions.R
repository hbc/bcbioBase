context("Data Functions")



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
        sampleData(rse_small, return = x)
    })
    expect_identical(
        lapply(list, class),
        list(
            "data.frame",
            structure("DataFrame", package = "S4Vectors"),
            "knitr_kable"
        )
    )
    # Check dimensions
    expect_identical(
        lapply(list, dim),
        list(
            c(6L, 8L),
            c(6L, 8L),
            NULL
        )
    )
    # Check rownames
    expect_identical(
        lapply(list, rownames),
        list(
            colnames(rse_small),
            colnames(rse_small),
            NULL
        )
    )
    # Ensure all columns are factor
    invisible(lapply(list[[1L]], function(x) {
        expect_is(x, "factor")
    }))
})

test_that("sampleData : Assignment method", {
    x <- rse_small
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



# uniteInterestingGroups =======================================================
test_that("uniteInterestingGroups : Single interesting group", {
    data <- uniteInterestingGroups(mtcars, interestingGroups = "gear")
    expect_identical(
        data[["interestingGroups"]],
        as.factor(mtcars[["gear"]])
    )
})

test_that("uniteInterestingGroups : Two interesting groups", {
    data <- uniteInterestingGroups(
        mtcars,
        interestingGroups = c("gear", "carb")
    )
    expect_identical(
        head(data[["interestingGroups"]]),
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
