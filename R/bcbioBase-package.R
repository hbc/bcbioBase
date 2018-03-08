#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @import methods
#' @importFrom rlang !!! !! .data abort inform sym syms warn
#' @importFrom utils globalVariables
"_PACKAGE"



globalVariables(".")
metadataPriorityCols <- c("sampleID", "sampleName", "description")



#' Project Directory Grep Pattern
#' @keywords internal
#' @export
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
