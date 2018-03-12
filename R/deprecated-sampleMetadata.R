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
