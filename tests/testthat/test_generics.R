context("S4 Generics")

test_that("No methods defined", {
    s4 <- list(
        bcbio,
        `bcbio<-`,
        metrics,
        plotGene,
        plotQC
    )
    methods <- vapply(
        X = s4,
        FUN = function(x) {
            showMethods(x, printTo = FALSE) %>%
                .[[2L]]
        },
        FUN.VALUE = "character"
    )
    expect_true(all(grepl("<No methods>", methods)))
    invisible(lapply(
        X = s4,
        FUN = function(x) {
            expect_error(
                x(),
                "unable to find an inherited method for function"
            )
        }
    ))
})
