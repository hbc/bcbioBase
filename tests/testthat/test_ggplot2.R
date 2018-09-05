context("ggplot2")



test_that("bcbio_geom_abline", {
    object <- bcbio_geom_abline(xintercept = 1L)
    expect_is(object, "Layer")

    object <- bcbio_geom_abline(yintercept = 1L)
    expect_is(object, "Layer")

    # Require single xintercept or yintercept
    error <- "Either `xintercept` or `yintercept` is required"
    expect_error(
        object = bcbio_geom_abline(),
        regexp = error
    )
    expect_error(
        object = bcbio_geom_abline(xintercept = 1L, yintercept = 1L),
        regexp = error
    )
})



test_that("bcbio_geom_label", {
    object <- bcbio_geom_label()
    expect_is(object, "Layer")
})



test_that("bcbio_geom_label_average", {
    # Normal mode.
    data <- data.frame(
        sampleName = c("sample1", "sample2"),
        counts = seq_len(8L)
    )
    object <- bcbio_geom_label_average(data, col = "counts")
    expect_is(object, "Layer")

    # Aggregate mode, for facet wrapping.
    data <- data.frame(
        sampleName = c("sample1", "sample2"),
        aggregate = "sample",
        counts = seq_len(8L)
    )
    object <- bcbio_geom_label_average(data, col = "counts")
    expect_is(object, "Layer")
})



test_that("bcbio_geom_label_repel", {
    object <- bcbio_geom_label_repel()
    expect_is(object, "Layer")

    # Single color mode.
    object <- bcbio_geom_label_repel(color = "orange")
    expect_identical(
        object = object[["aes_params"]][["colour"]],
        expected = "orange"
    )
})
