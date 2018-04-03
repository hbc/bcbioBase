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
#' @param color Colors to use for plot. Defaults to the [viridis()]
#'   palette.
#' @param legendColor Colors to use for legend labels. Defaults to the
#'   [viridis()] palette.
#' @param title *Optional.* Plot title.
#' @param ... Passthrough arguments to [pheatmap::pheatmap()].
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Show heatmap and return `list`, containing `gtable`.
#'
#' @examples
#' # SummarizedExperiment ====
#' plotHeatmap(rse_small)
#'
#' # matrix ====
#' mat <- assay(rse_small)
#' plotHeatmap(mat)
NULL



# Constructors =================================================================
.prepareAnnotationCol <- function(object) {
    # pheatmap requires `NA` argument if empty
    if (!has_dims(object)) {
        return(NA)
    }
    assertHasRownames(object)
    object %>%
        as.data.frame() %>%
        # Remove sample name columns
        .[, setdiff(colnames(.), metadataPriorityCols), drop = FALSE] %>%
        rownames_to_column() %>%
        # Ensure all columns are factor
        mutate_all(factor) %>%
        column_to_rownames()
}



.plotHeatmap.matrix <- function(  # nolint
    object,
    scale = c("row", "column", "none"),
    annotationCol = NULL,
    clusterCols = TRUE,
    clusterRows = TRUE,
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
    scale <- match.arg(scale)
    assertFormalAnnotationCol(object, annotationCol)
    annotationCol <- .prepareAnnotationCol(annotationCol)
    assert_is_a_bool(clusterCols)
    assert_is_a_bool(clusterRows)
    assertIsHexColorFunctionOrNULL(color)
    assertIsHexColorFunctionOrNULL(legendColor)
    assertIsAStringOrNULL(borderColor)
    assertIsAStringOrNULL(title)

    if (scale == "row") {
        # Filter out any zero count rows
        object <- object %>%
            .[rowSums(.) > 0L, , drop = FALSE]
    }

    # Define colors for each annotation column, if desired
    if (is.data.frame(annotationCol) && is.function(legendColor)) {
        annotationColors <- lapply(
            X = seq_along(colnames(annotationCol)),
            FUN = function(a) {
                col <- annotationCol[[a]] %>%
                    levels()
                colors <- annotationCol[[a]] %>%
                    levels() %>%
                    length() %>%
                    legendColor
                names(colors) <- col
                colors
            }
        ) %>%
            set_names(colnames(annotationCol))
    } else {
        annotationColors <- NULL
    }

    # If `color = NULL`, use the pheatmap default
    nColor <- 256L
    if (!is.function(color)) {
        color <- colorRampPalette(rev(
            brewer.pal(n = 7L, name = "RdYlBu")
        ))(nColor)
    } else {
        color <- color(nColor)
    }

    if (is.null(borderColor)) {
        borderColor <- NA
    }

    # Dynamic column and row labeling
    if (ncol(object) <= 50L) {
        showColnames <- TRUE
    } else{
        showColnames <- FALSE
    }
    if (nrow(object) <= 50L) {
        showRownames <- TRUE
    } else {
        showRownames <- FALSE
    }

    # pheatmap will error if `NULL` title is passed as `main`
    if (is.null(title)) {
        title <- ""
    }

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
    # Sanitize all argument names into snake case
    names(args) <- snake(names(args))
    assert_is_subset(names(args), formalArgs(pheatmap))

    do.call(pheatmap, args)
}



# Methods ======================================================================
#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("dgCMatrix"),
    .plotHeatmap.matrix
)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("dgTMatrix"),
    .plotHeatmap.matrix
)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("matrix"),
    .plotHeatmap.matrix
)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("SummarizedExperiment"),
    function(
        object,
        scale = c("row", "column", "none"),
        annotationCol,
        clusterCols = TRUE,
        clusterRows = TRUE,
        color = viridis,
        legendColor = viridis,
        borderColor = NULL,
        title = NULL,
        ...
    ) {
        scale <- match.arg(scale)
        if (missing(annotationCol)) {
            annotationCol <- colData(object)
        }
        plotHeatmap(
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
