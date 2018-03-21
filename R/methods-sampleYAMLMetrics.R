#' Sample Metrics from YAML File
#'
#' @note bcbio doesn't calculate sample metrics when run in fast mode.
#'
#' @name sampleYAMLMetrics
#' @family YAML Functions
#' @author Michael Steinbaugh
#'
#' @inherit sampleYAML
#'
#' @examples
#' yaml <- basejump::readYAML(
#'     "http://bcbiobase.seq.cloud/project-summary.yaml"
#' )
#' sampleYAMLMetrics(yaml) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_if
.sampleYAMLMetrics <- function(yaml) {
    # Early return on NULL metrics (fast mode)
    fastMode <- "Fast mode detected: No sample metrics were calculated"

    if (is.null(yaml[["samples"]][[1L]][["summary"]][["metrics"]])) {
        warn(fastMode)
        return(NULL)
    }

    data <- sampleYAML(
        yaml = yaml,
        keys = c("summary", "metrics")
    )
    assert_is_tbl(data)

    if (identical(colnames(data), "description")) {
        warn(fastMode)
        return(NULL)
    }

    # Fix numerics set as characters
    numericAsCharacter <- function(x) {
        any(grepl(x = x, pattern = "^[0-9\\.]+$"))
    }

    data %>%
        mutate_if(is.factor, as.character) %>%
        mutate_if(numericAsCharacter, as.numeric) %>%
        mutate_if(is.character, as.factor) %>%
        as.data.frame %>%
        set_rownames(make.names(.[["description"]], unique = TRUE)) %>%
        # Drop any metadata columns
        .[,
          sort(setdiff(
              x = colnames(.),
              y = c(metadataPriorityCols, "name")
          )),
          drop = FALSE]
}



# Methods ======================================================================
#' @rdname sampleYAMLMetrics
#' @export
setMethod(
    "sampleYAMLMetrics",
    signature("list"),
    .sampleYAMLMetrics
)
