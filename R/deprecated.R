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



# v0.2.0 =======================================================================
#' @rdname defunct
#' @export
checkInterestingGroups <- function(...) {
    .Defunct("assertFormalInterestingGroups")
}

#' @rdname defunct
#' @export
prepareSampleMetadata <- function(...) {
    .Defunct()
}

#' @rdname deprecated
#' @export
readLogFile <- function(...) {
    .Deprecated("readLog")
    readLog(...)
}

#' @rdname deprecated
#' @export
readSampleMetadataFile <- function(...) {
    .Deprecated("readSampleData")
    readSampleData(...)
}



# v0.2.4 =======================================================================
#' @rdname deprecated
#' @export
sampleMetadata <- function(object, ...) {
    .Deprecated("sampleData")
    sampleData(object, ...)
}

#' @rdname deprecated
#' @export
`sampleMetadata<-` <- function(object, value) {
    .Deprecated("sampleData<-")
    sampleData(object) <- value
    object
}

#' @rdname defunct
#' @export
sampleYAML <- function(...) {
    .Defunct("readYAMLSampleData or readYAMLSampleMetrics")
}

#' @rdname defunct
#' @export
sampleYAMLMetadata <- function(...) {
    .Defunct("readYAMLSampleData")
}

#' @rdname defunct
#' @export
sampleYAMLMetrics <- function(...) {
    .Defunct("readSampleYAMLMetrics")
}



# v0.2.5 =======================================================================
#' @rdname defunct
#' @export
prepareSampleData <- function(...) {
    .Defunct()
}



# v0.3.0 =======================================================================
#' @rdname defunct
#' @export
assertFormalAnnotationCol <- function(...) {
    .Defunct()
}

#' @importFrom basejump assertFormalInterestingGroups
#' @export
basejump::assertFormalInterestingGroups

#' @importFrom basejump convertGenesToSymbols
#' @export
basejump::convertGenesToSymbols

#' @importFrom basejump gene2symbol
#' @export
basejump::gene2symbol

#' @importFrom basejump interestingGroups
#' @export
basejump::interestingGroups

#' @importFrom basejump interestingGroups<-
#' @export
basejump::`interestingGroups<-`

#' @importFrom basejump sampleData
#' @export
basejump::sampleData

#' @importFrom basejump sampleData<-
#' @export
basejump::`sampleData<-`

#' @importFrom basejump sampleNames
#' @export
basejump::sampleNames

#' @importFrom basejump sanitizeSampleData
#' @export
basejump::sanitizeSampleData

#' @importFrom basejump selectSamples
#' @export
basejump::selectSamples

#' @importFrom basejump uniteInterestingGroups
#' @export
basejump::uniteInterestingGroups



# v0.3.2 =======================================================================
#' @importFrom basejump separatorBar
#' @export
basejump::separatorBar

#' @importFrom basejump updateMessage
#' @export
basejump::updateMessage

#' @rdname defunct
#' @export
setGeneric(
    "bcbio",
    function(object, ...) {
        .Deprecated()
        standardGeneric("bcbio")
    }
)

#' @rdname defunct
#' @export
setGeneric(
    "bcbio<-",
    function(object, ..., value) {
        .Deprecated()
        standardGeneric("bcbio<-")
    }
)



# v0.4.0 =======================================================================
#' @rdname deprecated
#' @export
setGeneric(
    name = "flatFiles",
    def = getGeneric("flatFiles", package = "basejump")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "metrics",
    def = getGeneric("metrics", package = "basejump")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotCorrelationHeatmap",
    def = getGeneric("plotCorrelationHeatmap", package = "basejump")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotGene",
    def = getGeneric("plotGene", package = "basejump")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotHeatmap",
    def = getGeneric("plotHeatmap", package = "basejump")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotQC",
    def = getGeneric("plotQC", package = "basejump")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotQuantileHeatmap",
    def = getGeneric("plotQuantileHeatmap", package = "basejump")
)

#' @rdname deprecated
#' @export
prepareSummarizedExperiment <- function(...) {
    basejump::makeSummarizedExperiment(...)
}



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



# nolint end
# nocov end
