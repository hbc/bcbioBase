.assertIsSampleData <- function(object) {
    # Stop on detection of blacklisted columns.
    intersect <- intersect(colnames(object), metadataBlacklist)
    if (length(intersect) > 0L) {
        stop(paste0(
            paste("Invalid columns:", toString(intersect)), "\n",
            "Recommended columns:\n",
            "  - fileName: FASTQ file name (optional, but recommended).\n",
            "  - description: Sample description per file (required).\n",
            "  - sampleName: Unique sample name",
            " (multiplexed samples only).\n",
            "The `description` column is sanitized into the sample ID ",
            "for demultiplexed samples."
        ))
    }
    # Check for required columns.
    assert_is_subset(x = "description", y = colnames(object))
    TRUE
}



.makeSampleData <- function(object) {
    object <- as(object, "DataFrame")
    # Note that we want to call `.assertIsSampleData()` earlier, so we can
    # detect if the user is attempting to input automatic columns, such as
    # "revcomp". At this point, automatic columsns are allowed, so we don't want
    # to check for them again here.
    assert_is_subset("description", colnames(object))
    # Set sampleName from description, if necessary.
    if (!"sampleName" %in% colnames(object)) {
        object[["sampleName"]] <- object[["description"]]
    }
    # Set the sample IDs as rownames, using the "description" column. Here we
    # are ensuring that the names are valid in R. This allows for input of
    # samples beginning with numbers or containing hyphens for example, which
    # aren't valid names in R. Note that periods and underscores are valid.
    rownames(object) <- makeNames(object[["description"]], unique = TRUE)
    makeSampleData(object)
}
