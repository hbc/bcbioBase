#' Prepare Sample Data
#'
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `data.frame`.
#' @export
prepareSampleData <- function(object) {
    assert_has_dimnames(object)
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
    object[["sampleID"]] <- makeNames(object[["sampleID"]], unique = TRUE)

    object %>%
        as_tibble() %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        .[, unique(c(metadataPriorityCols, colnames(.)))] %>%
        arrange(!!!syms(metadataPriorityCols)) %>%
        as.data.frame() %>%
        set_rownames(.[["sampleID"]])
}
