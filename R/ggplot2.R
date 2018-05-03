# nolint start



#' Repulsive Textual Annotations
#'
#' Modified bcbio version of [ggrepel::geom_label_repel()]. If advanced
#' customization of the text labels is required, simply use the ggrepel version
#' instead. This is a convenience function.
#'
#' @author Michael Steinbaugh
#'
#' @seealso [ggrepel::geom_label_repel()].
#'
#' @return `ggproto`.
#' @export
#'
#' @examples
#' x <- geom_label_bcbio()
#' class(x)
geom_label_bcbio <- function(
    data = NULL,
    mapping = NULL,
    color = NULL,
    size = 4L
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



geom_cutoff_line_bcbio <- function(
    xintercept = NULL,
    yintercept = NULL,
    alpha = 0.75,
    color = "black",
    linetype = "dashed",
    size = 1L
) {
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



# nolint end
