#' Sample Metadata from YAML
#'
#' @name sampleYAMLMetadata
#' @family YAML Functions
#' @author Michael Steinbaugh
#'
#' @inherit sampleYAML
#'
#' @examples
#' yaml <- basejump::readYAML(
#'     "http://bcbiobase.seq.cloud/project-summary.yaml"
#' )
#' sampleYAMLMetadata(yaml) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_all
.sampleYAMLMetadata <- function(yaml) {
    sampleYAML(yaml = yaml, keys = "metadata") %>%
        prepareSampleData()
}



# Methods ======================================================================
#' @rdname sampleYAMLMetadata
#' @export
setMethod(
    "sampleYAMLMetadata",
    signature("list"),
    .sampleYAMLMetadata
)
