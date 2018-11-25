.makeSampleData <- function(object) {
    assert_is_subset("description", colnames(object))
    assert_are_disjoint_sets(colnames(object), metadataBlacklist)
    # Set sampleName from description, if necessary.
    if (!"sampleName" %in% colnames(object)) {
        object[["sampleName"]] <- object[["description"]]
    }
    # Ensure `sampleID` has valid names here. This allows for input of samples
    # beginning with numbers or containing hyphens for example, which aren't
    # valid names in R.
    rownames(object) <- makeNames(object[["sampleName"]], unique = TRUE)
    makeSampleData(object)
}
