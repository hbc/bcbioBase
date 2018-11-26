#' Get GTF File from YAML
#'
#' @author Michael Steinbaugh
#' @inheritParams basejump::params
#' @inheritParams params
#' @export
#'
#' @return `string` or `NULL`. File path if the file exists. `NULL` if the file
#'   does not exist.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "summary.yaml")
#' yaml <- import(file)
#' x <- getGTFFileFromYAML(yaml)
#' print(x)
getGTFFileFromYAML <- function(yaml) {
    .assertIsSummaryYAML(yaml)
    # Assume all samples are using the same GTF file.
    file <- yaml %>%
        .[["samples"]] %>%
        .[[1L]] %>%
        .[["genome_resources"]] %>%
        .[["rnaseq"]] %>%
        .[["transcripts"]]
    assert_is_a_string(file)
    message(paste("bcbio GTF file:", file))
    if (!file.exists(file)) {
        message("GTF file does not exist. Skipping.")
        NULL
    } else {
        file
    }
}
