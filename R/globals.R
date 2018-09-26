globalVariables(".")

version <- packageVersion("bcbioBase")
versionDir <- paste0("v", version$major, ".", version$minor)

# TODO Migrate S3 to a bucket that supports https.

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioBaseCacheURL
bcbioBaseCacheURL <- paste0("http://bcbiobase.seq.cloud/", versionDir)

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
