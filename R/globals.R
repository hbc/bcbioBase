globalVariables(".")

metadataBlacklist <- c(
    "description",
    "fileName",
    "genomeBuild",
    "name",  # from metrics YAML
    "samRef",
    "sampleID",
    "sampleName"
)

#' Lane Grep Pattern
#' @keywords internal
#' @export
#' @examples
#' lanePattern
lanePattern <- "_L(\\d{3})"

#' Project Directory Grep Pattern
#' @keywords internal
#' @export
#' @examples
#' projectDirPattern
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"

#' Separator Bar
#' @keywords internal
#' @export
#' @examples
#' separatorBar
separatorBar <- "============================================================"

#' Update Message
#' @keywords internal
#' @export
#' @examples
#' updateMessage
updateMessage <- "Run `updateObject()` to update your object"
