context("Plot Heatmap Functions")



# All heatmap functions ========================================================
fxns <- c(
    "plotCorrelationHeatmap",
    "plotHeatmap",
    "plotQuantileHeatmap"
)

test_that("SummarizedExperiment", {
    invisible(lapply(fxns, function(f) {
        object <- rse_small
        f <- get(f)
        p <- f(object)

        # Expect pheatmap return
        expect_is(p, "list")
        expect_identical(names(p), heatmapList)

        # Test that plots contain annotation data
        gtable <- p[["gtable"]]
        expect_true("annotation_legend" %in% gtable[["layout"]][["name"]])

        # Test color palette support
        expect_silent(f(object, color = NULL, legendColor = NULL))
    }))
})



# plotHeatmap ==================================================================
test_that("plotHeatmap : Turn off column names for many samples", {
    p <- plotHeatmap(
        object = matrix(seq(1L:1000L), ncol = 100L)
    )
    layout <- p[["gtable"]][["layout"]]
    expect_is(layout, "data.frame")
    expect_false("col_names" %in% layout[["name"]])
})
