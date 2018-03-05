context("Generics")

withMethods <- c(
    "prepareSummarizedExperiment",
    "prepareTemplate",
    "sampleYAML",
    "sampleYAMLMetadata",
    "sampleYAMLMetrics"
)
withoutMethods <- c(
    "bcbio",
    "bcbio<-",
    "flatFiles",
    "interestingGroups",
    "interestingGroups<-",
    "metrics",
    "plotDot",
    "plotGene",
    "plotQC",
    "plotViolin",
    "sampleMetadata",
    "sampleMetadata<-",
    "selectSamples",
    "tpm"
)

test_that("S4 generics", {
    generics <- lapply(
        X = c(withMethods, withoutMethods),
        FUN = get)
    expect_true(all(
        vapply(
            X = generics,
            FUN = isS4,
            FUN.VALUE = logical(1L))
    ))
})

test_that("No methods defined", {
    generics <- lapply(
        X = withoutMethods,
        FUN = get)

    methods <- vapply(
        X = withoutMethods,
        FUN = function(x) {
            showMethods(x, printTo = FALSE) %>%
                .[[2L]]
        },
        FUN.VALUE = "character")
    expect_true(all(grepl("<No methods>", methods)))

    invisible(lapply(
        X = generics,
        FUN = function(x) {
            expect_error(
                x(),
                "unable to find an inherited method for function"
            )
        }
    ))
})
