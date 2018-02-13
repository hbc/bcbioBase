#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @import methods
#' @importFrom rlang !!! !! .data abort inform sym syms warn
"_PACKAGE"

#' @importFrom utils globalVariables
globalVariables(".")

metadataPriorityCols <- c("sampleID", "sampleName", "description")
