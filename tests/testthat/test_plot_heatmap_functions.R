context("Plot Heatmap Functions")

fxns <- c(
    "plotCorrelationHeatmap",
    "plotHeatmap",
    "plotQuantileHeatmap"
)

test_that("SummarizedExperiment", {
    invisible(lapply(fxns, function(f) {
        object <- rse_bcb
        f <- get(f)
        p <- f(object)

        # Expect pheatmap return
        expect_is(p, "list")
        expect_identical(names(p), heatmapList)

        # Test that plots contain annotation data
        gtable <- p[["gtable"]]
        expect_true("annotation_legend" %in% gtable[["layout"]][["name"]])

        # Test color and title support
        expect_is(
            f(
                object = object,
                color = NULL,
                legendColor = NULL,
                title = NULL
            ),
            "list"
        )
        expect_is(
            f(
                object = object,
                color = viridis,
                legendColor = viridis
            ),
            "list"
        )
    }))
})

test_that("matrix", {
    invisible(lapply(fxns, function(f) {
        object = as.matrix(datasets::USArrests)
        f <- get(f)
        p <- f(object)

        # Expect pheatmap return
        expect_is(p, "list")
        expect_identical(names(p), heatmapList)

        # Test that plots do not contain annotation data
        gtable <- p[["gtable"]]
        expect_false("annotation_legend" %in% gtable[["layout"]][["name"]])
    }))
})

test_that("Invalid pheatmap passthrough", {
    expect_error(
        plotHeatmap(rse_bcb, show_colnames = FALSE),
        "Define formalArgs in camel case: show_colnames"
    )
})
