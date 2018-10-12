# FIXME Need to add working example.



#' Get GTF File from YAML
#'
#' @param yaml `list`. Project summary YAML.
#'
#' @return `string`. GTF file path.
#' @export
#'
#' @examples
#' \dontrun{
#' getGTFFileFromYAML(yaml)
#' }
getGTFFileFromYAML <- function(yaml) {
    assert_is_list(yaml)
    # Assume all samples are using the same GTF file.
    file <- yaml %>%
        .[["samples"]] %>%
        .[[1L]] %>%
        .[["genome_resources"]] %>%
        .[["rnaseq"]] %>%
        .[["transcripts"]]
    assert_is_a_string(file)
    assert_all_are_existing_files(file)
    assert_are_identical(basename(file), "ref-transcripts.gtf")
    file
}
