#' Sample Metadata from YAML
#'
#' @rdname sampleYAMLMetadata
#' @name sampleYAMLMetadata
#' @family YAML Utilities
#'
#' @inherit sampleYAML
#'
#' @examples
#' url <- paste(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "project-summary.yaml",
#'     sep = "/")
#' yaml <- basejump::readYAML(url)
#' sampleYAMLMetadata(yaml) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_all
.sampleYAMLMetadata <- function(yaml) {
    sampleYAML(yaml = yaml, keys = "metadata") %>%
        mutate_all(as.factor) %>%
        .prepareSampleMetadata()
}



# Methods ======================================================================
#' @rdname sampleYAMLMetadata
#' @export
setMethod(
    "sampleYAMLMetadata",
    signature("list"),
    .sampleYAMLMetadata)
