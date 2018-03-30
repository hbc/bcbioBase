# nocov start

#' Defunct or Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @return No value.
NULL



# v0.2.0 =======================================================================
# `annotable()` made defunct in basejump v0.4.0

#' @rdname deprecated
#' @export
prepareSampleMetadata <- function(...) {
    .Deprecated("prepareSampleData")
    prepareSampleData(...)
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

# nocov end
