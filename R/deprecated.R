## nocov start
## nolint start



#' @name defunct
#' @inherit acidroxygen::defunct description examples return seealso title
#' @inheritParams acidroxygen::params
#' @keywords internal
NULL



#' @name deprecated
#' @inherit acidroxygen::deprecated description examples return seealso title
#' @inheritParams acidroxygen::params
#' @keywords internal
NULL



## v0.5.14 =====================================================================
#' @importFrom basejump metadataBlacklist
#' @export
basejump::metadataBlacklist



# v0.6.9 =======================================================================
#' @importFrom basejump readSampleData
#' @export
basejump::readSampleData

#' @importFrom basejump readTx2Gene
#' @export
basejump::readTx2Gene



# v0.6.10 ======================================================================
#' @rdname deprecated
#' @export
readDataVersions <- function(...) {
    .Deprecated("importDataVersions")
    importDataVersions(...)
}

#' @rdname deprecated
#' @export
readProgramVersions <- function(...) {
    .Deprecated("importProgramVersions")
    importProgramVersions(...)
}



## nolint end
## nocov end
