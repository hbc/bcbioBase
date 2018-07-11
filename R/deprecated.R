# nocov start



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



# nocov end
