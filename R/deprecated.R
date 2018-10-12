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



# v0.2.0 =======================================================================
#' @rdname defunct
#' @export
checkInterestingGroups <- function(...) {
    .Defunct("basejump::assertFormalInterestingGroups")
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
    .Deprecated("basejump::sampleData")
    basejump::sampleData(object, ...)
}

#' @rdname deprecated
#' @export
`sampleMetadata<-` <- function(object, value) {
    .Deprecated("basejump::sampleData<-")
    basejump::sampleData(object) <- value
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



# v0.4.0 =======================================================================
#' @rdname deprecated
#' @importFrom basejump makeSummarizedExperiment
#' @export
prepareSummarizedExperiment <- function(...) {
    .Deprecated("basejump::makeSummarizedExperiment")
    basejump::makeSummarizedExperiment(...)
}



# v0.4.1 =======================================================================
#' @importFrom basejump prepareTemplate
#' @export
basejump::prepareTemplate



# v0.4.2 =======================================================================
#' @importFrom basejump lanePattern
#' @export
basejump::lanePattern



# v0.99.0 ======================================================================
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
