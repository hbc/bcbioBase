# nocov start
# nolint start



#' @name defunct
#' @inherit basejump::defunct
#' @keywords internal
NULL

#' @name deprecated
#' @inherit basejump::deprecated
#' @keywords internal
NULL



# v0.4.1 =======================================================================
#' @importFrom basejump prepareTemplate
#' @export
basejump::prepareTemplate



# v0.5.0 =======================================================================
#' @importFrom basejump lanePattern
#' @export
basejump::lanePattern

#' @importFrom basejump minimalSampleData
#' @export
basejump::minimalSampleData

#' @rdname deprecated
#' @importFrom basejump basejump_geom_abline
#' @export
bcbio_geom_abline <- function(...) {
    .Deprecated("basejump::basejump_geom_abline")
    basejump::basejump_geom_abline(...)
}

#' @rdname deprecated
#' @importFrom basejump basejump_geom_label
#' @export
bcbio_geom_label <- function(...) {
    .Deprecated("basejump::basejump_geom_label")
    basejump::basejump_geom_label(...)
}

#' @rdname deprecated
#' @importFrom basejump basejump_geom_label_average
#' @export
bcbio_geom_label_average <- function(...) {
    .Deprecated("basejump::basejump_geom_label_average")
    basejump::basejump_geom_label_average(...)
}

#' @rdname deprecated
#' @importFrom basejump basejump_geom_label_repel
#' @export
bcbio_geom_label_repel <- function(...) {
    .Deprecated("basejump::basejump_geom_label_repel")
    basejump::basejump_geom_label_repel(...)
}

#' @rdname deprecated
#' @importFrom basejump import
#' @export
readLog <- function(file) {
    .Deprecated("basejump::import")
    basejump::import(file)
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



# 0.5.14 =======================================================================
#' @importFrom basejump metadataBlacklist
#' @export
basejump::metadataBlacklist



# nolint end
# nocov end
