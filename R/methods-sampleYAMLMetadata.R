#' Sample Metadata from YAML
#'
#' @rdname sampleYAMLMetadata
#' @name sampleYAMLMetadata
#' @family YAML Utilities
#'
#' @inherit sampleYAML
#'
#' @examples
#' url <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "bcbio",
#'     "project-summary.yaml")
#' yaml <- basejump::readYAML(url)
#' sampleYAMLMetadata(yaml)
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_all
.sampleYAMLMetadata <- function(yaml) {
    data <- sampleYAML(yaml = yaml, keys = "metadata")
    # TODO Abort here?
    if (is.null(data)) {
        warn("Failed to obtain sample metadata")
        return(NULL)
    }
    data %>%
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
