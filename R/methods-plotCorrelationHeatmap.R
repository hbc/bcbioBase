#' Plot Correlation Heatmap
#'
#' Construct a correlation heatmap comparing the columns of the matrix.
#'
#' @name plotCorrelationHeatmap
#' @family Plot Heatmap Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams plotHeatmap
#' @param method Correlation coefficient (or covariance) method to be computed.
#'   Defaults to "`pearson`" but "`spearman`" can also be used. Consult the
#'   [stats::cor()] documentation for more information.
#' @param clusteringMethod Clustering method. Accepts the same values as
#'   [stats::hclust()].
#'
#' @seealso
#' - [stats::cor()].
#' - [stats::hclust()].
#' - [pheatmap::pheatmap()].
#'
#' @return Show heatmap and invisibly return a `list` of the components.
#'
#' @examples
#' # SummarizedExperiment ====
#' plotCorrelationHeatmap(rse_small)
#'
#' # matrix ====
#' mat <- assay(rse_small)
#' plotCorrelationHeatmap(mat)
NULL



# Methods ======================================================================
#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("matrix"),
    function(
        object,
        method = c("pearson", "spearman"),
        clusteringMethod = "ward.D2",
        annotationCol = NULL,
        color = viridis,
        legendColor = NULL,
        borderColor = NULL,
        title = TRUE,
        ...
    ) {
        assert_has_dims(object)
        assert_all_are_greater_than(nrow(object), 1L)
        assert_all_are_greater_than(ncol(object), 1L)
        method <- match.arg(method)
        assertIsHexColorFunctionOrNULL(color)
        assertIsHexColorFunctionOrNULL(legendColor)
        assertIsAStringOrNULL(borderColor)
        if (!is_a_string(borderColor)) {
            borderColor <- NA
        }
        if (isTRUE(title)) {
            title <- paste(method, "correlation")
        } else if (!is_a_string(title)) {
            title <- NA
        }

        # Correlation matrix
        mat <- cor(object, method = method)

        annotationCol <- .pheatmapAnnotationCol(annotationCol)
        assertFormalAnnotationCol(object, annotationCol)
        annotationColors <- .pheatmapAnnotationColors(
            annotationCol = annotationCol,
            legendColor = legendColor
        )
        color <- .pheatmapColor(color)

        # Return pretty heatmap with modified defaults
        args <- list(
            "mat" = mat,
            "annotationCol" = annotationCol,
            "annotationColors" = annotationColors,
            "borderColor" = borderColor,
            "clusteringMethod" = clusteringMethod,
            "clusteringDistanceRows" = "correlation",
            "clusteringDistanceCols" = "correlation",
            "color" = color,
            "main" = title,
            "showColnames" = TRUE,
            "showRownames" = TRUE,
            ...
        )
        args <- .pheatmapArgs(args)
        do.call(pheatmap, args)
    }
)



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("dgCMatrix"),
    getMethod("plotCorrelationHeatmap", "matrix")
)



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("dgTMatrix"),
    getMethod("plotCorrelationHeatmap", "matrix")
)



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("SummarizedExperiment"),
    function(object, ...) {
        plotCorrelationHeatmap(
            object = assay(object),
            annotationCol = colData(object),
            ...
        )
    }
)
