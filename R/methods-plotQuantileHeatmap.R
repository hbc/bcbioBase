#' Plot Heatmap with Quantile Breaks
#'
#' @name plotQuantileHeatmap
#' @family Plot Heatmap Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams plotHeatmap
#' @param n The number of breaks to create.
#'
#' @return Show heatmap and invisibly return a `list` of the components.
#'
#' @examples
#' # SummarizedExperiment ====
#' plotQuantileHeatmap(rse_small)
#'
#' # matrix ====
#' mat <- assay(rse_small)
#' plotQuantileHeatmap(mat)
NULL



# Constructors =================================================================
#' Create Breaks Based on Quantiles of the Data
#'
#' @keywords internal
#' @noRd
#'
#' @param x Numeric vector.
#' @param n The number of breaks to create.
#'
#' @return A vector of `n` quantile breaks.
.quantileBreaks <- function(object, n = 5L) {
    assert_is_matrix(object)
    assert_is_an_integer(n)
    assert_all_are_positive(n)
    q <- quantile(object, probs = seq(0L, 1L, length.out = n))
    q[!duplicated(q)]
}



.plotQuantileHeatmap.matrix <- function(  # nolint
    object,
    n = 5L,
    annotationCol = NULL,
    clusterCols = TRUE,
    clusterRows = TRUE,
    showColnames = TRUE,
    showRownames = FALSE,
    color = viridis,
    legendColor = viridis,
    borderColor = NULL,
    title = NULL,
    ...
) {
    assert_has_dims(object)
    assert_all_are_greater_than(nrow(object), 1L)
    assert_all_are_greater_than(ncol(object), 1L)
    object <- as.matrix(object)
    assertIsAnImplicitInteger(n)
    assertFormalAnnotationCol(object, annotationCol)
    assert_is_a_bool(clusterCols)
    assert_is_a_bool(clusterRows)
    assertIsHexColorFunctionOrNULL(color)
    assertIsHexColorFunctionOrNULL(legendColor)
    assertIsAStringOrNULL(borderColor)
    if (!is_a_string(borderColor)) {
        borderColor <- NA
    }
    assertIsAStringOrNULL(title)
    if (!is_a_string(title)) {
        title <- NA
    }

    # Calculate the quantile breaks
    breaks <- .quantileBreaks(object, n = n)

    annotationCol <- .pheatmapAnnotationCol(annotationCol)
    annotationColors <- .pheatmapAnnotationColors(
        annotationCol = annotationCol,
        legendColor = legendColor
    )
    color <- .pheatmapColor(color)

    # Return pretty heatmap with modified defaults
    args <- list(
        "mat" = object,
        "annotationCol" = annotationCol,
        "annotationColors" = annotationColors,
        "borderColor" = borderColor,
        "breaks" = breaks,
        "clusterCols" = clusterCols,
        "clusterRows" = clusterRows,
        "color" = color,
        "main" = title,
        "showColnames" = showColnames,
        "showRownames" = showRownames,
        ...
    )
    args <- .pheatmapArgs(args)
    do.call(pheatmap, args)
}



# Methods ======================================================================
#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("dgCMatrix"),
    .plotQuantileHeatmap.matrix
)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("dgTMatrix"),
    .plotQuantileHeatmap.matrix
)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("matrix"),
    .plotQuantileHeatmap.matrix
)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("SummarizedExperiment"),
    function(
        object,
        n = 5L,
        annotationCol,
        clusterCols = TRUE,
        clusterRows = TRUE,
        color = viridis,
        legendColor = viridis,
        borderColor = NULL,
        title = NULL,
        ...
    ) {
        if (missing(annotationCol)) {
            annotationCol <- sampleData(object)
        }
        plotQuantileHeatmap(
            object = assay(object),
            annotationCol = annotationCol,
            clusterCols = clusterCols,
            clusterRows = clusterRows,
            color = color,
            legendColor = legendColor,
            borderColor = borderColor,
            title = title,
            ...
        )
    }
)
