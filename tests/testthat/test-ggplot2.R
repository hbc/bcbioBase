context("ggplot2")



test_that("bcbio_geom_abline", {
    g <- bcbio_geom_abline(xintercept = 1L)
    expect_is(g, "Layer")

    g <- bcbio_geom_abline(yintercept = 1L)
    expect_is(g, "Layer")

    # Require single xintercept or yintercept
    e <- "Either `xintercept` or `yintercept` is required"
    expect_error(bcbio_geom_abline(), e)
    expect_error(bcbio_geom_abline(xintercept = 1L, yintercept = 1L), e)
})



test_that("bcbio_geom_label", {
    g <- bcbio_geom_label()
    expect_is(g, "Layer")
})



test_that("bcbio_geom_label_average", {
    # Normal mode
    data <- data.frame(
        sampleName = c("sample1", "sample2"),
        counts = seq_len(8L)
    )
    g <- bcbio_geom_label_average(data, col = "counts")
    expect_is(g, "Layer")

    # Aggregate mode, for facet wrapping
    data <- data.frame(
        sampleName = c("sample1", "sample2"),
        aggregate = "sample",
        counts = seq_len(8L)
    )
    g <- bcbio_geom_label_average(data, col = "counts")
    expect_is(g, "Layer")
})



test_that("bcbio_geom_label_repel", {
    g <- bcbio_geom_label_repel()
    expect_is(g, "Layer")

    # Single color mode
    g <- bcbio_geom_label_repel(color = "orange")
    expect_identical(
        g[["aes_params"]][["colour"]],
        "orange"
    )
})
