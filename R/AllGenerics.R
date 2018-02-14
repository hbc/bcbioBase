#' S4 Generics
#'
#' @rdname AllGenerics
#' @name AllGenerics
#'
#' @inheritParams general
#'
#' @return Varies, depending upon the method.
NULL



#' @rdname AllGenerics
#' @export
setGeneric("bcbio", function(object, ...) {
    standardGeneric("bcbio")
})



#' @rdname AllGenerics
#' @export
setGeneric("bcbio<-", function(object, ..., value) {
    standardGeneric("bcbio<-")
})



#' @rdname AllGenerics
#' @export
setGeneric("flatFiles", function(object, ...) {
    standardGeneric("flatFiles")
})



#' @rdname AllGenerics
#' @export
setGeneric("interestingGroups", function(object, ...) {
    standardGeneric("interestingGroups")
})



#' @rdname AllGenerics
#' @export
setGeneric("interestingGroups<-", function(object, ..., value) {
    standardGeneric("interestingGroups<-")
})



#' @rdname AllGenerics
#' @export
setGeneric("metrics", function(object, ...) {
    standardGeneric("metrics")
})



#' @rdname AllGenerics
#' @export
setGeneric("plotDot", function(object, ...) {
    standardGeneric("plotDot")
})



#' @rdname AllGenerics
#' @export
setGeneric("plotGene", function(object, ...) {
    standardGeneric("plotGene")
})



#' @rdname AllGenerics
#' @export
setGeneric("plotQC", function(object, ...) {
    standardGeneric("plotQC")
})



#' @rdname AllGenerics
#' @export
setGeneric("plotViolin", function(object, ...) {
    standardGeneric("plotViolin")
})



#' @rdname prepareSummarizedExperiment
#' @export
setGeneric("prepareSummarizedExperiment", function(assays, ...) {
    standardGeneric("prepareSummarizedExperiment")
})



#' @rdname prepareTemplate
#' @export
setGeneric("prepareTemplate", function(object, ...) {
    standardGeneric("prepareTemplate")
})



#' @rdname AllGenerics
#' @export
setGeneric("sampleMetadata", function(object, ...) {
    standardGeneric("sampleMetadata")
})



#' @rdname sampleYAML
#' @export
setGeneric("sampleYAML", function(yaml, keys, ...) {
    standardGeneric("sampleYAML")
})



#' @rdname sampleYAMLMetadata
#' @export
setGeneric("sampleYAMLMetadata", function(yaml, ...) {
    standardGeneric("sampleYAMLMetadata")
})



#' @rdname sampleYAMLMetrics
#' @export
setGeneric("sampleYAMLMetrics", function(yaml, ...) {
    standardGeneric("sampleYAMLMetrics")
})



#' @rdname AllGenerics
#' @export
setGeneric("selectSamples", function(object, ...) {
    standardGeneric("selectSamples")
})



#' @rdname AllGenerics
#' @export
setGeneric("tpm", function(object) {
    standardGeneric("tpm")
})
