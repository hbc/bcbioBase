# nocov start
# nolint start



#' @name defunct
#' @inherit basejump::defunct
#' @inheritParams params
#' @keywords internal
NULL

#' @name deprecated
#' @inherit basejump::deprecated
#' @inheritParams params
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
    # .Defunct("acidplots::acid_geom_abline")  # nolint
    requireNamespace("acidplots", quietly = TRUE)
    acidplots::acid_geom_abline(...)
}

#' @rdname defunct
#' @export
bcbio_geom_label <- function(...) {
    # .Defunct("acidplots::acid_geom_label")  # nolint
    requireNamespace("acidplots", quietly = TRUE)
    acidplots::acid_geom_label(...)
}

#' @rdname defunct
#' @export
bcbio_geom_label_average <- function(...) {
    # .Defunct("acidplots::acid_geom_label_average")  # nolint
    requireNamespace("acidplots", quietly = TRUE)
    acidplots::acid_geom_label_average(...)
}

#' @rdname defunct
#' @export
bcbio_geom_label_repel <- function(...) {
    # .Defunct("acidplots::acid_geom_label_repel")  # nolint
    requireNamespace("acidplots", quietly = TRUE)
    acidplots::acid_geom_label_repel(...)
}

#' @rdname deprecated
#' @export
readLog <- function(file) {
    # .Defunct("basejump::import")  # nolint
    requireNamespace("basejump", quietly = TRUE)
    basejump::import(file)
}

#' @rdname defunct
#' @export
readTx2gene <- function(...) {
    .Defunct("readTx2Gene")
}

#' @rdname defunct
#' @export
readYAMLSampleData <- function(file) {
    # .Defunct("getSampleDataFromYAML")  # nolint
    requireNamespace("basejump", quietly = TRUE)
    yaml <- basejump::import(file)
    getSampleDataFromYAML(yaml)
}

#' @rdname defunct
#' @export
readYAMLSampleMetrics <- function(file) {
    # .Defunct("getMetricsFromYAML")  # nolint
    requireNamespace("basejump", quietly = TRUE)
    yaml <- basejump::import(file)
    getMetricsFromYAML(yaml)
}



# 0.5.14 =======================================================================
#' @importFrom basejump metadataBlacklist
#' @export
basejump::metadataBlacklist



# nolint end
# nocov end
