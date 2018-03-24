context("Plot Functions")



# plotHeatmap ==================================================================
test_that("plotHeatmap : matrix", {
    p <- plotHeatmap(mat)
    expect_is(p, "list")
    expect_identical(names(p), heatmapList)
})

test_that("Annotation columns support", {
    annotationCol <- data.frame(
        genotype = c(
            "wildtype",
            "wildtype",
            "mutant",
            "mutant"
        ),
        row.names = colnames(mat)
    )
    p <- plotHeatmap(mat, annotationCol = annotationCol)
    expect_identical(names(p), heatmapList)
    expect_length(p[["gtable"]], 10L)
})

test_that("plotHeatmap : Default color palette", {
    expect_silent(plotHeatmap(mat, color = NULL))
})

test_that("plotHeatmap : Turn off columns for many samples", {
    mat <- matrix(seq(1L:1000L), ncol = 100L)
    p <- plotHeatmap(mat)
    expect_identical(names(p), heatmapList)
    expect_length(p[["gtable"]], 5L)
})

test_that("plotHeatmap : Matrix dimensions are too small", {
    expect_error(
        plotHeatmap(matrix(seq(1L:10L), ncol = 1L)),
        "is_greater_than : ncol"
    )
    expect_error(
        plotHeatmap(matrix(seq(1L:10L), nrow = 1L)),
        "is_greater_than : nrow"
    )
})



# plotQuantileHeatmap ==========================================================
test_that("plotQuantileHeatmap : matrix", {
    mat <- matrix(seq_len(10000L), nrow = 100L, ncol = 100L)
    p <- plotQuantileHeatmap(mat)
    expect_is(p, "list")
    expect_identical(names(p), c("quantiles", "plot"))
    expect_identical(names(p[["plot"]]), heatmapList)
    expect_length(p[["plot"]][["gtable"]], 5L)
})

test_that("plotQuantileHeatmap : Annotation columns support", {
    annotationCol <- data.frame(
        genotype = c(
            "wildtype",
            "wildtype",
            "mutant",
            "mutant"
        ),
        row.names = colnames(mat)
    )
    p <- plotQuantileHeatmap(mat, annotationCol = annotationCol)
    expect_identical(names(p[["plot"]]), heatmapList)
    expect_length(p[["plot"]][["gtable"]], 10L)
})

test_that("plotQuantileHeatmap : Default color palette", {
    p <- plotQuantileHeatmap(mat, color = NULL)
    expect_is(p, "list")
})

test_that("plotQuantileHeatmap : Matrix dimensions are too small", {
    expect_error(
        plotQuantileHeatmap(matrix(seq(1L:10L), ncol = 1L)),
        "is_greater_than : ncol"
    )
    expect_error(
        plotQuantileHeatmap(matrix(seq(1L:10L), nrow = 1L)),
        "is_greater_than : nrow"
    )
})
