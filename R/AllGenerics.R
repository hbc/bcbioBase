#' S4 Generics
#'
#' @name AllGenerics
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return Varies, depending upon the method.
NULL



#' @rdname AllGenerics
#' @export
setGeneric(
    "bcbio",
    function(object, ...) {
        standardGeneric("bcbio")
    }
)



#' @rdname AllGenerics
#' @export
setGeneric(
    "bcbio<-",
    function(object, ..., value) {
        standardGeneric("bcbio<-")
    }
)



#' @rdname flatFiles
#' @export
setGeneric(
    "flatFiles",
    function(object, ...) {
        standardGeneric("flatFiles")
    }
)



#' @rdname gene2symbol
#' @export
setGeneric(
    "gene2symbol",
    function(object) {
        standardGeneric("gene2symbol")
    }
)



#' @rdname interestingGroups
#' @export
setGeneric(
    "interestingGroups",
    function(object, ...) {
        standardGeneric("interestingGroups")
    }
)



#' @rdname interestingGroups
#' @export
setGeneric(
    "interestingGroups<-",
    function(object, ..., value) {
        standardGeneric("interestingGroups<-")
    }
)



#' @rdname AllGenerics
#' @export
setGeneric(
    "metrics",
    function(object, ...) {
        standardGeneric("metrics")
    }
)



#' @rdname plotCorrelationHeatmap
#' @export
setGeneric(
    "plotCorrelationHeatmap",
    function(object, ...) {
        standardGeneric("plotCorrelationHeatmap")
    }
)



#' @rdname AllGenerics
#' @export
setGeneric(
    "plotGene",
    function(object, ...) {
        standardGeneric("plotGene")
    }
)



#' @rdname plotHeatmap
#' @export
setGeneric(
    "plotHeatmap",
    function(object, ...) {
        standardGeneric("plotHeatmap")
    }
)



#' @rdname plotQuantileHeatmap
#' @export
setGeneric(
    "plotQuantileHeatmap",
    function(object, ...) {
        standardGeneric("plotQuantileHeatmap")
    }
)



#' @rdname AllGenerics
#' @export
setGeneric(
    "plotQC",
    function(object, ...) {
        standardGeneric("plotQC")
    }
)



#' @rdname sampleData
#' @export
setGeneric(
    "sampleData",
    function(object, ...) {
        standardGeneric("sampleData")
    }
)



#' @rdname sampleData
#' @export
setGeneric(
    "sampleData<-",
    function(object, ..., value) {
        standardGeneric("sampleData<-")
    }
)



#' @rdname selectSamples
#' @export
setGeneric(
    "selectSamples",
    function(object, ...) {
        standardGeneric("selectSamples")
    }
)



#' @rdname uniteInterestingGroups
#' @export
setGeneric(
    "uniteInterestingGroups",
    function(object, ...) {
        standardGeneric("uniteInterestingGroups")
    }
)
