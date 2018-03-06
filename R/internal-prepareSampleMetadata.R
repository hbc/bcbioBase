#' @importFrom dplyr arrange everything mutate mutate_if select
#' @importFrom magrittr set_rownames
#' @importFrom tibble as_tibble
.prepareSampleMetadata <- function(object) {
    assert_has_dimnames(object)
    object <- as_tibble(object)
    assert_is_subset("description", colnames(object))

    # Set `sampleName`, if necessary
    if (!"sampleName" %in% colnames(object)) {
        object[["sampleName"]] <- object[["description"]]
    }

    # Set `sampleID`, if necessary
    if (!"sampleID" %in% colnames(object)) {
        object[["sampleID"]] <- object[["sampleName"]]
    }

    # Ensure `sampleID` has valid names. This allows for input of samples
    # beginning with numbers or containing hyphens for example, which aren't
    # valid names in R.
    object[["sampleID"]] <- gsub(
        x = make.names(object[["sampleID"]], unique = TRUE),
        pattern = "\\.",
        replacement = "_"
    )

    object %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        select(metadataPriorityCols, everything()) %>%
        arrange(!!!syms(metadataPriorityCols)) %>%
        as.data.frame() %>%
        set_rownames(.[["sampleID"]])
}
