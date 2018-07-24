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
        expect_is(p, "pheatmap")
        expect_identical(names(p), pheatmapList)

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
            "pheatmap"
        )
        # Hexadecimal color functions (e.g. viridis)
        expect_is(
            f(
                object = object,
                color = viridis,
                legendColor = viridis
            ),
            "pheatmap"
        )
        # Hexadecimal color palettes (e.g. RColorBrewer)
        color <- colorRampPalette(brewer.pal(n = 11, name = "PuOr"))(256)
        expect_is(
            f(
                object = object,
                color = color
            ),
            "pheatmap"
        )
        # Disable interesting groups
        expect_is(
            f(
                object = object,
                interestingGroups = NULL
            ),
            "pheatmap"
        )
    }))
})

test_that("Invalid pheatmap passthrough", {
    expect_error(
        plotHeatmap(rse_bcb, show_colnames = FALSE),
        "Define formalArgs in camel case: show_colnames"
    )
})
