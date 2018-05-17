context("ggplot2")



test_that("bcbio_geom_abline", {
    g <- bcbio_geom_abline(xintercept = 1L)
    expect_is(g, "Layer")

    g <- bcbio_geom_abline(yintercept = 1L)
    expect_is(g, "Layer")

    expect_error(
        bcbio_geom_abline(xintercept = 1L, yintercept = 1L),
        "Specify only `xintercept` or `yintercept` but not both"
    )
})



test_that("bcbio_geom_label", {
    g <- bcbio_geom_label()
    expect_is(g, "Layer")
})



test_that("bcbio_geom_label_average", {
    data <- data.frame(
        sampleName = c("sample1", "sample1", "sample1", "sample1"),
        counts = c(1, 2, 3, 4)
    )
    g <- bcbio_geom_label_average(data, col = "counts")
    expect_is(g, "Layer")
})



test_that("bcbio_geom_label_repel", {
    g <- bcbio_geom_label_repel()
    expect_is(g, "Layer")
})
