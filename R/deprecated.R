# nocov start

#' Defunct or Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @return No value.
NULL

# v0.2.0 =======================================================================
#' @rdname deprecated
#' @export
sampleMetadata <- function(object, ...) {
    .Deprecated("sampleData")
    sampleData(object, ...)
}

#' @rdname deprecated
#' @export
`sampleMetadata<-` <- function(object, ..., value) {
    .Deprecated("sampleData<-")
    sampleData(object, ..., value)
}

#' @rdname deprecated
#' @export
flatFiles <- function(object) {
    .Deprecated("as(object, \"list\")")
    return(NULL)
    stopifnot(is(object, "SummarizedExperiment"))
    as(object, "list")
}



# future =======================================================================
# checkInterestingGroups

# nocov end
