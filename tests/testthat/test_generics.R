context("bcbioGenerics")

names <- c(
    "bcbio",
    "bcbio<-",
    "flatFiles",
    "interestingGroups",
    "metrics",
    "plotDot",
    "plotGene",
    "plotQC",
    "plotViolin",
    "sampleMetadata",
    "selectSamples",
    "tpm"
)
generics <- lapply(names, get)

test_that("Exported generics", {
    classes <- vapply(
        X = generics,
        FUN = class,
        FUN.VALUE = "character")
    expect_true(all(classes == "nonstandardGenericFunction"))
})

test_that("No methods defined", {
    methods <- vapply(
        X = names,
        FUN = function(x) {
            showMethods(x, printTo = FALSE) %>%
                .[[2L]]
        },
        FUN.VALUE = "character"
    )
    names(methods) <- names
    expect_true(all(grepl("<No methods>", methods)))
})
