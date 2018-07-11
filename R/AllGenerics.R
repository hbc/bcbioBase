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
    "plotDot",
    function(object, ...) {
        standardGeneric("plotDot")
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



#' @rdname AllGenerics
#' @export
setGeneric(
    "plotQC",
    function(object, ...) {
        standardGeneric("plotQC")
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
    "plotViolin",
    function(object, ...) {
        standardGeneric("plotViolin")
    }
)
