globalVariables(".")

#' Lane Grep Pattern
#' @keywords internal
#' @export
#' @examples
#' lanePattern
lanePattern <- "_L(\\d{3})"

#' Metadata Priority Columns
#' @keywords internal
#' @export
#' @examples
#' metadataPriorityCols
metadataPriorityCols <- c("sampleID", "sampleName", "description")

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
