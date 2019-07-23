## Sample metadata assert check for goalie engine.
.isSampleData <- function(object) {
    ok <- isAny(object, c("data.frame", "DataFrame"))
    if (!isTRUE(ok)) return(ok)

    ## Check for blacklisted columns.
    intersect <- intersect(colnames(object), metadataBlacklist)
    ok <- !hasLength(intersect)
    if (!isTRUE(ok)) {
        return(false(paste0(
            "Blacklist detection: ", toString(intersect), "\n\n",
            "Recommended columns:\n",
            "  - fileName: FASTQ file name (optional, but recommended).\n",
            "  - description: Sample description per file (required).\n",
            "  - sampleName: Unique sample name",
            " (multiplexed samples only).\n\n",
            "Refer to bcbioBase::readSampleData() for formatting requirements."
        )))
    }

    ## Check for required columns (e.g. description).
    required <- "description"
    ok <- isSubset(required, colnames(object))
    if (!isTRUE(ok)) {
        setdiff <- setdiff(required, colnames(object))
        return(false(paste0(
            "Required columns missing: ", setdiff, "\n\n",
            "Refer to bcbioBase::readSampleData() for formatting requirements."
        )))
    }

    TRUE
}



## Wrap `makeSampleData()` call with bcbio-specific additions.
.makeSampleData <- function(object) {
    object <- as(object, "DataFrame")

    ## Note that we want to call `.assertIsSampleData()` earlier, so we can
    ## detect if the user is attempting to input automatic columns, such as
    ## "revcomp". At this point, automatic columsns are allowed, so we don't
    ## want to check for them again here.
    assert(isSubset("description", colnames(object)))

    ## Set `sampleName` from `description`, if necessary.
    if (!"sampleName" %in% colnames(object)) {
        object[["sampleName"]] <- object[["description"]]
    }

    ## Set the sample IDs as rownames, using the "description" column. Here we
    ## are ensuring that the names are valid in R. This allows for input of
    ## samples beginning with numbers or containing hyphens for example, which
    ## aren't valid names in R. Note that periods and underscores are valid.
    rownames(object) <- makeNames(object[["description"]], unique = TRUE)

    makeSampleData(object)
}
