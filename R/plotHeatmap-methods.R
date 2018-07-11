#' Plot Heatmap
#'
#' Construct a simple heatmap. By default, row- and column-wise hierarchical
#' clustering is performed using the Ward method, but this behavior can be
#' overrided by setting `clusterRows` or `clusterCols` to `FALSE`.
#'
#' @name plotHeatmap
#' @family Plot Functions
#'
#' @inheritParams general
#' @param scale `character` indicating whether the values should be centered and
#'   scaled in either the row or column direction ("`row`", "`column`"), or
#'   remain unscaled ("`none`").
#' @param annotationCol *Optional.* `data.frame` that defines annotation
#'   mappings for the columns.
#' @param clusterRows,clusterCols `logical` determining if rows or columns
#'   should be arranged with hierarchical clustering. Alternatively, can define
#'   an `hclust` object.
#' @param showRownames,showColnames Show row or column names.
#' @param treeheightRow,treeheightCol Size of the row and column dendrograms.
#'   Use "`0`" to disable.
#' @param color Hexadecimal color function to use for plot. We recommend any of
#'   these from the viridis package:
#'   - [viridis::viridis()] (*default*).
#'   - [viridis::inferno()].
#'   - [viridis::magma()].
#'   - [viridis::plasma()].
#' @param legendColor Hexadecimal color function to use for legend labels.
#' @param borderColor *Optional.* Border color. Disabled by default for
#'   improved aesthetics.
#' @param title *Optional.* Plot title.
#' @param ... Passthrough arguments to [pheatmap::pheatmap()]. The names of the
#'   arguments should be formatted in camel case, not snake case.
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Show heatmap and invisibly return a `list` of the components.
#'
#' @examples
#' # SummarizedExperiment ====
#' plotHeatmap(rse_dds, interestingGroups = "condition")
#'
#' # Disable column clustering
#' plotHeatmap(rse_dds, clusterCols = FALSE)
NULL



# Methods ======================================================================
#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("SummarizedExperiment"),
    function(
        object,
        interestingGroups,
        scale = c("row", "column", "none"),
        clusterRows = TRUE,
        clusterCols = TRUE,
        showRownames = FALSE,
        showColnames = TRUE,
        treeheightRow = 0L,
        treeheightCol = 50L,
        color = viridis,
        legendColor = NULL,
        borderColor = NULL,
        title = NULL,
        ...
    ) {
        assert_all_are_greater_than(nrow(object), 1L)
        assert_all_are_greater_than(ncol(object), 1L)
        scale <- match.arg(scale)
        assert_is_a_bool(clusterCols)
        assert_is_a_bool(clusterRows)
        assert_is_a_bool(showColnames)
        assert_is_a_bool(showRownames)
        assert_is_a_number(treeheightRow)
        assert_is_a_number(treeheightCol)
        assert_all_are_non_negative(treeheightRow, treeheightCol)
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

        object <- suppressWarnings(convertGenesToSymbols(object))

        mat <- as.matrix(assay(object))
        if (scale == "row") {
            # Filter out any zero count rows
            mat <- mat[rowSums(mat) > 0L, , drop = FALSE]
        }

        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        if (length(interestingGroups)) {
            annotationCol <- colData(object)[, interestingGroups, drop = FALSE]
        } else {
            annotationCol <- NULL
        }

        # Use `sampleName`, if defined
        sampleName <- colData(object)[["sampleName"]]
        if (length(sampleName)) {
            colnames(mat) <- sampleName
            if (length(annotationCol)) {
                rownames(annotationCol) <- sampleName
            }
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
            mat = mat,
            annotationCol = annotationCol,
            annotationColors = annotationColors,
            borderColor = borderColor,
            clusterCols = clusterCols,
            clusterRows = clusterRows,
            color = color,
            main = title,
            scale = scale,
            showColnames = showColnames,
            showRownames = showRownames,
            treeheightCol = treeheightCol,
            treeheightRow = treeheightRow,
            ...
        )
        args <- .pheatmapArgs(args)
        do.call(pheatmap, args)
    }
)
