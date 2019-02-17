globalVariables(".")



yamlFlatCols <- c(
    "description",
    "genome_build",
    "sam_ref"
)

metricsBlacklist <- c(
    camel(yamlFlatCols),
    "name"
)

#' Metadata Blacklist
#' @keywords internal
#' @export
#' @examples
#' metadataBlacklist
metadataBlacklist <- sort(c(
    metricsBlacklist,
    "aggregate",
    "fileName",
    "index",
    "qualityFormat",
    "revcomp",
    "sampleID",
    "sequence",
    "sequenceLength"
))



#' Project Directory Grep Pattern
#' @keywords internal
#' @export
#' @examples
#' projectDirPattern
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
