# nocov start

#' Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @return No value.
NULL

# v0.1.5 =======================================================================
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



# future =======================================================================
# checkInterestingGroups

# nocov end
