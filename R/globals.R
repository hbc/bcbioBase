globalVariables(".")

packageVersion <- packageVersion("bcbioBase")

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioBaseCacheURL
bcbioBaseCacheURL <- paste0(
    "http://bcbiobase.seq.cloud/",
    "v", packageVersion$major, ".", packageVersion$minor  # nolint
)

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
projectDirPattern <- "^([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})_([^/]+)$"
