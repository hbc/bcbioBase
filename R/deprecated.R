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



# v0.4.0 =======================================================================
#' @rdname deprecated
#' @export
prepareSummarizedExperiment <- function(...) {
    basejump::makeSummarizedExperiment(...)
}



# v0.4.1 =======================================================================
#' @importFrom basejump prepareTemplate
#' @export
basejump::prepareTemplate



# nocov end
