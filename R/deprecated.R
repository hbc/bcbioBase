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

#' @rdname defunct
#' @export
readLogFile <- function(...) {
    .Defunct("readLog")
}

#' @rdname defunct
#' @export
readSampleMetadataFile <- function(...) {
    .Defunct("readSampleData")
}



# v0.2.4 =======================================================================
#' @rdname defunct
#' @export
sampleMetadata <- function(object, ...) {
    .Defunct("sampleData")
}

#' @rdname defunct
#' @export
`sampleMetadata<-` <- function(object, value) {
    .Defunct("sampleData<-")
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

#' @rdname deprecated
#' @name assertFormalInterestingGroups
#' @importFrom basejump assertFormalInterestingGroups
#' @usage NULL
#' @export
NULL

#' @rdname deprecated
#' @name convertGenesToSymbols
#' @importFrom basejump convertGenesToSymbols
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name gene2symbol
#' @importFrom basejump gene2symbol
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name interestingGroups
#' @importFrom basejump interestingGroups
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name interestingGroups<-
#' @importFrom basejump interestingGroups<-
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name sampleData
#' @importFrom basejump sampleData
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name sampleData<-
#' @importFrom basejump sampleData<-
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name sampleNames
#' @importFrom basejump sampleNames
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name sanitizeSampleData
#' @importFrom basejump sanitizeSampleData
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name selectSamples
#' @importFrom basejump selectSamples
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name uniteInterestingGroups
#' @importFrom basejump uniteInterestingGroups
#' @export
#' @usage NULL
NULL



# v0.3.2 =======================================================================
#' @rdname deprecated
#' @name separatorBar
#' @importFrom basejump separatorBar
#' @export
#' @usage NULL
NULL

#' @rdname deprecated
#' @name updateMessage
#' @importFrom basejump updateMessage
#' @export
#' @usage NULL
NULL

#' @rdname defunct
#' @export
setGeneric(
    "bcbio",
    function(object, ...) {
        .Defunct()
    }
)

#' @rdname defunct
#' @export
setGeneric(
    "bcbio<-",
    function(object, ..., value) {
        .Defunct()
    }
)



# v0.4.0 =======================================================================
#' @rdname deprecated
#' @export
setGeneric(
    name = "flatFiles",
    def = getGeneric("flatFiles")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "metrics",
    def = getGeneric("metrics")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotCorrelationHeatmap",
    def = getGeneric("plotCorrelationHeatmap")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotGene",
    def = getGeneric("plotGene")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotHeatmap",
    def = getGeneric("plotHeatmap")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotQC",
    def = getGeneric("plotQC")
)

#' @rdname deprecated
#' @export
setGeneric(
    name = "plotQuantileHeatmap",
    def = getGeneric("plotQuantileHeatmap")
)

#' @rdname deprecated
#' @export
prepareSummarizedExperiment <- function(...) {
    basejump::makeSummarizedExperiment(...)
}



# v0.4.1 =======================================================================
#' @rdname deprecated
#' @name prepareTemplate
#' @importFrom basejump prepareTemplate
#' @export
#' @usage NULL
NULL



# nocov end
