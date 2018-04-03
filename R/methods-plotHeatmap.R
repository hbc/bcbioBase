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



# Constructors =================================================================
# Sanitize formals into snake case and abort on duplicates.
# Duplicates may arise if user is mixing and matching camel/snake case.
.pheatmapArgs <- function(args) {
    assert_is_list(args)
    assert_has_names(args)
    # Abort on snake case formatted formalArgs
    invalidNames <- grep("[._]", names(args), value = TRUE)
    if (length(invalidNames)) {
        abort(paste(
            "Define formalArgs in camel case:",
            toString(invalidNames)
        ))
    }
    # Abort on duplicate arguments
    names(args) <- snake(names(args))
    if (any(duplicated(names(args)))) {
        abort(paste(
            "Duplicate formalArgs detected:",
            toString(camel(names(args)[duplicated(names(args))]))
        ))
    }
    assert_is_subset(names(args), formalArgs(pheatmap))
    args
}



.plotHeatmap.matrix <- function(  # nolint
    object,
    scale = c("row", "column", "none"),
    annotationCol = NULL,
    clusterCols = TRUE,
    clusterRows = TRUE,
    showColnames = TRUE,
    showRownames = TRUE,
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
    assert_is_a_bool(showColnames)
    assert_is_a_bool(showRownames)
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
    args <- .pheatmapArgs(args)
    do.call(pheatmap, args)
}



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
