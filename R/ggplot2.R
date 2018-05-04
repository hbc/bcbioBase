# nolint start

#' ggplot2 Functions
#'
#' Convenience functions with modified defaults for
#' [ggplot2](http://ggplot2.org).
#'
#' @name ggplot2
#' @author Michael Steinbaugh
#'
#' @inheritParams ggplot2::geom_label
#'
#' @seealso
#' - [ggplot2::geom_label()].
#' - [ggrepel::geom_label_repel()].
#'
#' @return `ggproto`.
NULL



#' @rdname ggplot2
#' @section bcbio_geom_abline:
#'
#' Horizontal or vertical cutoff line.
#'
#' @export
#'
#' @examples
#' # x-axis line
#' geom <- bcbio_geom_abline(xintercept = 1L)
#' geom
#'
#' # y-axis line
#' geom <- bcbio_geom_abline(yintercept = 1L)
#' geom
bcbio_geom_abline <- function(
    xintercept = NULL,
    yintercept = NULL
) {
    alpha <- 0.75
    color <- "black"
    linetype <- "dashed"
    size <- 1L
    if (is.null(xintercept) && is.null(yintercept)) {
        stop("`xintercept` and `yintercept` are both NULL")
    } else if (is.numeric(xintercept) && is.numeric(yintercept)) {
        stop("Specifcly only `xintercept` or `yintercept` as numeric")
    } else if (is.numeric(xintercept)) {
        geom_vline(
            xintercept = xintercept,
            alpha = alpha,
            color = color,
            linetype = linetype,
            size = size
        )
    } else if (is.numeric(yintercept)) {
        geom_hline(
            yintercept = yintercept,
            alpha = alpha,
            color = color,
            linetype = linetype,
            size = size
        )
    }
}



#' @rdname ggplot2
#' @section bcbio_geom_label:
#'
#' Modified version of [ggplot2::geom_label()].
#'
#' @export
#'
#' @examples
#' geom <- bcbio_geom_label()
#' geom
bcbio_geom_label <- function(
    data = NULL,
    mapping = NULL,
    ...
) {
    geom_label(
        data = data,
        mapping = mapping,
        alpha = 0.75,
        color = "white",
        fill = "black",
        fontface = "bold",
        label.padding = unit(0.2, "lines"),
        label.size = NA,
        show.legend = FALSE,
        ...
    )
}



#' @rdname ggplot2
#' @section bcbio_geom_label_average:
#'
#' Add average labels to a plot. For example, `col` can be `nGene`. Median or
#' mean values are always calculated per sample (`sampleName`).
#'
#' @param data `data.frame`.
#' @param col Column.
#' @param fun Function to use for average calculation. Currently supports
#'   "`mean`" or "`median`".
#' @param digits Number of significant digits to use. Defaults to rounded.
#'
#' @export
#'
#' @examples
#' geom <- bcbio_geom_label_average(
#'     data = colData(rse_bcb),
#'     col = "exonicRate",
#'     fun = "median"
#' )
#' geom
bcbio_geom_label_average <- function(
    data,
    col,
    fun = c("median", "median"),
    digits = 0L,
    ...
) {
    data <- as.data.frame(data)
    assert_is_a_string(col)
    assert_is_subset(col, colnames(data))
    assert_is_an_integer(digits)
    fun <- match.arg(fun)
    fun <- get(fun)
    assert_is_function(fun)

    data <- aggregate(
        formula = as.formula(paste(col, "sampleName", sep = " ~ ")),
        data = data,
        FUN = fun
    )
    data[["roundedAverage"]] <- round(data[[col]], digits = digits)

    # Add `aggregate` column for facet wrapping, if necessary
    if ("aggregate" %in% colnames(data)) {
        sampleFacet <- data[, c("sampleName", "aggregate")] %>%
            unique()
        data <- merge(
            x = data,
            y = sampleFacet,
            by = "sampleName",
            all.x = TRUE
        )
    }

    bcbio_geom_label(
        data = data,
        mapping = aes_string(label = "roundedAverage"),
        ...
    )
}



#' @rdname ggplot2
#' @section bcbio_geom_label_repel:
#'
#' Repulsive textual annotations. Modified bcbio version of
#' [ggrepel::geom_label_repel()]. If advanced customization of the text labels
#' is required, simply use the ggrepel version instead.
#'
#' @export
#'
#' @examples
#' geom <- bcbio_geom_label_repel()
#' geom
bcbio_geom_label_repel <- function(
    data = NULL,
    mapping = NULL,
    color = NULL,
    size = 4L,
    ...
) {
    geom <- geom_label_repel(
        data = data,
        mapping = mapping,
        arrow = arrow(length = unit(0.01, "npc")),
        box.padding = unit(0.5, "lines"),
        fill = "white",
        fontface = "bold",
        force = 1L,
        point.padding = unit(0.75, "lines"),
        segment.size = 0.5,
        show.legend = FALSE,
        size = size
    )
    if (is.character(color)) {
        geom[["aes_params"]][["colour"]] <- color
    }
    geom
}



# nolint end
