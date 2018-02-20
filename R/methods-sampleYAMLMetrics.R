#' Sample Metrics from YAML File
#'
#' @rdname sampleYAMLMetrics
#' @name sampleYAMLMetrics
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
#' sampleYAMLMetrics(yaml) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_if
.sampleYAMLMetrics <- function(yaml) {
    # The fast mode RNA-seq pipeline doesn't report metrics generated from
    # STAR featureCounts output with MultiQC
    fastMode <- "Fast mode detected: No sample metrics were calculated"

    # Early return on NULL metrics (fast mode)
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
        .prepareSampleMetadata() %>%
        # Drop any sample metadata ID columns
        .[, setdiff(
            x = colnames(.),
            y = c(metadataPriorityCols, "name")
        ), drop = FALSE]
}



# Methods ======================================================================
#' @rdname sampleYAMLMetrics
#' @export
setMethod(
    "sampleYAMLMetrics",
    signature("list"),
    .sampleYAMLMetrics)
