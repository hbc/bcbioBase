#' Plot Heatmap with Quantile Breaks
#'
#' @name plotQuantileHeatmap
#' @family Plot Heatmap Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams plotHeatmap
#' @param n The number of breaks to create.
#' @param legend Show the color legend.
#'
#' @return Show heatmap and invisibly return a `list` of the components.
#'
#' @examples
#' # SummarizedExperiment ====
#' plotQuantileHeatmap(rse_dds)
#'
#' # matrix ====
#' mat <- assay(rse_dds)
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
.quantileBreaks <- function(object, n = 10L) {
    assert_is_matrix(object)
    assert_is_an_integer(n)
    assert_all_are_positive(n)
    breaks <- quantile(object, probs = seq(0L, 1L, length.out = n))
    breaks[!duplicated(breaks)]
}



# Methods ======================================================================
#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("matrix"),
    function(
        object,
        n = 10L,
        annotationCol = NULL,
        clusterCols = TRUE,
        clusterRows = TRUE,
        showColnames = TRUE,
        showRownames = FALSE,
        legend = FALSE,
        color = viridis,
        legendColor = NULL,
        borderColor = NULL,
        title = NULL,
        ...
    ) {
        assert_has_dims(object)
        assert_all_are_greater_than(nrow(object), 1L)
        assert_all_are_greater_than(ncol(object), 1L)
        object <- as.matrix(object)
        assertIsAnImplicitInteger(n)
        n <- as.integer(n)
        assertFormalAnnotationCol(object, annotationCol)
        assert_is_a_bool(clusterCols)
        assert_is_a_bool(clusterRows)
        assertIsHexColorFunctionOrNULL(color)
        assert_is_a_bool(legend)
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
        assertFormalAnnotationCol(object, annotationCol)
        annotationColors <- .pheatmapAnnotationColors(
            annotationCol = annotationCol,
            legendColor = legendColor
        )
        color <- .pheatmapColor(color, n = length(breaks) - 1L)

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
            "legend" = legend,
            "legendBreaks" = breaks,
            "legendLabels" = round(breaks, digits = 2L),
            "main" = title,
            "scale" = "none",
            "showColnames" = showColnames,
            "showRownames" = showRownames,
            ...
        )
        args <- .pheatmapArgs(args)
        do.call(pheatmap, args)
    }
)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("dgCMatrix"),
    getMethod("plotQuantileHeatmap", "matrix")
)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("dgTMatrix"),
    getMethod("plotQuantileHeatmap", "matrix")
)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("SummarizedExperiment"),
    function(object, ...) {
        annotationCol <- sampleData(object, interestingGroups = NULL)
        plotQuantileHeatmap(
            object = assay(object),
            annotationCol = annotationCol,
            ...
        )
    }
)
