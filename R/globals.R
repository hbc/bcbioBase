globalVariables(".")

packageVersion <- packageVersion("bcbioBase")

#' Cache URL
#' @export
#' @examples
#' bcbioBaseCacheURL
bcbioBaseCacheURL <- paste0(
    "http://bcbiobase.seq.cloud/",
    "v", packageVersion$major, ".", packageVersion$minor  # nolint
)

#' Sample Metadata Blacklist
#' @export
#' @examples
#' metadataBlacklist
metadataBlacklist <- c(
    # Too vague.
    "ID", "Id", "id",
    # Generated automatically.
    "interestingGroups",
    # Use "sampleName" instead.
    "name",
    # Generated automatically from "sequence" column.
    "revcomp",
    # Used internally by dplyr.
    "rowname",
    # Use "sampleName" instead.
    "sample",
    # "sampleID" is set automatically, for multiplexed/cell-level data.
    "sampleID", "sampleId", "sampleid",
    # Use "sampleName" instead.
    "samplename"
)

#' Project Directory Grep Pattern
#' @export
#' @examples
#' projectDirPattern
projectDirPattern <- "^([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})_([^/]+)$"
