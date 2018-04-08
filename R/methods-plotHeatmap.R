#' Plot Heatmap
#'
#' Construct a simple heatmap. By default, row- and column-wise hierarchical
#' clustering is performed using the Ward method, but this behavior can be
#' overrided by setting `clusterRows` or `clusterCols` to `FALSE`.
#'
#' @name plotHeatmap
#' @family Plot Heatmap Functions
#'
#' @inheritParams general
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the row direction or the column direction, or none. Corresponding
#'   values are "row", "column" and "none".
#' @param annotationCol *Optional.* `data.frame` that defines annotation
#'   mappings for the columns.
#' @param borderColor Border color.
#' @param clusterCols Logical determining if columns should be arranged with
#'   hierarchical clustering. Alternatively, can define an `hclust` object.
#' @param clusterRows Logical determining if rows should be arranged with
#'   hierarchical clustering. Alternatively, can define an `hclust` object.
#' @param showColnames Show column names.
#' @param showRownames Show row names.
#' @param color Colors to use for plot. Defaults to the [viridis()]
#'   palette.
#' @param legendColor Colors to use for legend labels. Defaults to the
#'   [viridis()] palette.
#' @param title *Optional.* Plot title.
#' @param ... Passthrough arguments to [pheatmap::pheatmap()]. The names of the
#'   arguments must be formatted in camel case, not snake case.
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Show heatmap and invisibly return a `list` of the components.
#'
#' @examples
#' # SummarizedExperiment ====
#' plotHeatmap(rse_small)
#'
#' # matrix ====
#' mat <- assay(rse_small)
#' plotHeatmap(mat)
NULL



# Methods ======================================================================
#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("matrix"),
    function(
        object,
        scale = c("row", "column", "none"),
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
        object <- as.matrix(object)
        assert_has_dims(object)
        assert_all_are_greater_than(nrow(object), 1L)
        assert_all_are_greater_than(ncol(object), 1L)
        scale <- match.arg(scale)
        assert_is_a_bool(clusterCols)
        assert_is_a_bool(clusterRows)
        assert_is_a_bool(showColnames)
        assert_is_a_bool(showRownames)
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

        if (scale == "row") {
            # Filter out any zero count rows
            object <- object %>%
                .[rowSums(.) > 0L, , drop = FALSE]
        }

        annotationCol <- .pheatmapAnnotationCol(annotationCol)
        assertFormalAnnotationCol(object, annotationCol)
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
            "clusterCols" = clusterCols,
            "clusterRows" = clusterRows,
            "color" = color,
            "main" = title,
            "scale" = scale,
            "showColnames" = showColnames,
            "showRownames" = showRownames,
            ...
        )
        args <- .pheatmapArgs(args)
        do.call(pheatmap, args)
    }
)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("dgCMatrix"),
    getMethod("plotHeatmap", "matrix")
)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("dgTMatrix"),
    getMethod("plotHeatmap", "matrix")
)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("SummarizedExperiment"),
    function(object, ...) {
        plotHeatmap(
            object = assay(object),
            annotationCol = colData(object),
            ...
        )
    }
)
