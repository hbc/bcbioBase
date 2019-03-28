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

#' @rdname defunct
#' @export
bcbio_geom_abline <- function(...) {
    .Defunct("minimalism::acid_geom_abline")
}

#' @rdname defunct
#' @export
bcbio_geom_label <- function(...) {
    .Defunct("minimalism::acid_geom_label")
}

#' @rdname defunct
#' @export
bcbio_geom_label_average <- function(...) {
    .Defunct("minimalism::acid_geom_label_average")
}

#' @rdname defunct
#' @export
bcbio_geom_label_repel <- function(...) {
    .Defunct("minimalism::acid_geom_label_repel")
}

#' @rdname defunct
#' @export
readLog <- function(file) {
    .Defunct("basejump::import")
}

#' @rdname defunct
#' @export
readTx2gene <- function(...) {
    .Defunct("readTx2Gene")
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
