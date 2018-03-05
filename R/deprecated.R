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
sampleMetadata <- function(...) {
    .Deprecated("sampleData")
    sampleData(...)
}

#' @rdname deprecated
#' @export
`sampleMetadata<-` <- function(...) {
    .Deprecated("sampleData<-")
    sampleData(...)
}

# nocov end
