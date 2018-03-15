#' @rdname deprecated
#' @export
setGeneric("sampleMetadata", function(object, ...) {
    .Deprecated("sampleData")
    standardGeneric("sampleMetadata")
})

#' @rdname deprecated
#' @export
setGeneric("sampleMetadata<-", function(object, ..., value) {
    .Deprecated("sampleData<-")
    standardGeneric("sampleMetadata<-")
})

#' @rdname deprecated
#' @export
setMethod(
    "sampleMetadata",
    signature("SummarizedExperiment"),
    function(object, ...) {
        sampleData(object, ...)
    }
)

#' @rdname deprecated
#' @export
setMethod(
    "sampleMetadata<-",
    signature(
        object = "SummarizedExperiment",
        value = "ANY"
    ),
    function(object, ..., value) {
        value <- as(value, "DataFrame")
        sampleData(object) <- value
        object
    }
)
