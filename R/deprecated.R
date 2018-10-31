# nocov start
# nolint start



#' Defunct Functions
#'
#' @name defunct
#' @keywords internal
#'
#' @return No value.
NULL



#' Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @return No value.
NULL



# v0.99.0 ======================================================================
#' @importFrom basejump.experiment minimalSampleData
#' @export
basejump.experiment::minimalSampleData

#' @rdname deprecated
#' @importFrom basejump.plots basejump_geom_abline
#' @export
bcbio_geom_abline <- function(...) {
    .Deprecated("basejump.plots::basejump_geom_abline")
    basejump.plots::basejump_geom_abline(...)
}

#' @rdname deprecated
#' @importFrom basejump.plots basejump_geom_label
#' @export
bcbio_geom_label <- function(...) {
    .Deprecated("basejump.plots::basejump_geom_label")
    basejump.plots::basejump_geom_label(...)
}

#' @rdname deprecated
#' @importFrom basejump.plots basejump_geom_label_average
#' @export
bcbio_geom_label_average <- function(...) {
    .Deprecated("basejump.plots::basejump_geom_label_average")
    basejump.plots::basejump_geom_label_average(...)
}

#' @rdname deprecated
#' @importFrom basejump.plots basejump_geom_label_repel
#' @export
bcbio_geom_label_repel <- function(...) {
    .Deprecated("basejump.plots::basejump_geom_label_repel")
    basejump.plots::basejump_geom_label_repel(...)
}

#' @rdname deprecated
#' @importFrom basejump.io import
#' @export
readLog <- function(file) {
    .Deprecated("basejump.io::import")
    basejump.io::import(file)
}

#' @rdname deprecated
#' @export
readTx2gene <- function(...) {
    .Deprecated("readTx2Gene")
    readTx2Gene(...)
}

#' @rdname defunct
#' @export
readYAMLSampleData <- function(...) {
    .Defunct("getSampleDataFromYAML")
}

#' @rdname defunct
#' @export
readYAMLSampleMetrics <- function(...) {
    .Defunct("getMetricsFromYAML")
}



# nolint end
# nocov end
