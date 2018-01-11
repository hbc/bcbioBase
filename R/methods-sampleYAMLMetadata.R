#' Sample Metadata from YAML
#'
#' @rdname sampleYAMLMetadata
#' @name sampleYAMLMetadata
#' @family YAML Utilities
#'
#' @inherit sampleYAML
#' 
#' @inheritParams AllGenerics
#'
#' @examples
#' url <- file.path(
#'     "http://basejump.seq.cloud",
#'     "bcbio",
#'     "project-summary.yaml")
#' yaml <- basejump::readYAML(url)
#' sampleYAMLMetadata(yaml)
NULL



# Methods ======================================================================
#' @rdname sampleYAMLMetadata
#' @importFrom dplyr mutate_all
#' @export
setMethod(
    "sampleYAMLMetadata",
    signature("list"),
    function(yaml) {
        sampleYAML(
            yaml = yaml,
            keys = "metadata"
        ) %>%
            mutate_all(as.factor) %>%
            .prepareSampleMetadata()
    })
