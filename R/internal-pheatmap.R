# Sanitize formals into snake case and abort on duplicates.
# Duplicates may arise if user is mixing and matching camel/snake case.
.pheatmapArgs <- function(args) {
    assert_is_list(args)
    assert_has_names(args)
    # Abort on snake case formatted formalArgs
    invalidNames <- grep("[._]", names(args), value = TRUE)
    if (length(invalidNames)) {
        stop(paste(
            "Define formalArgs in camel case:",
            toString(invalidNames)
        ))
    }
    names(args) <- snake(names(args))
    assert_is_subset(names(args), formalArgs(pheatmap))
    args
}



# Automatically handle the annotation columns.
# Factors with a single level are automatically dropped.
.pheatmapAnnotationCol <- function(data) {
    # pheatmap requires `NA` argument if empty
    if (!has_dims(data)) {
        return(NA)
    }
    assertHasRownames(data)
    blacklist <- c("sampleName", metadataBlacklist)

    data <- data %>%
        as.data.frame() %>%
        # Remove sample name columns
        .[, setdiff(colnames(.), blacklist), drop = FALSE] %>%
        rownames_to_column() %>%
        # Ensure all strings are factor
        mutate_if(is.character, as.factor) %>%
        # Ensure unwanted columns like `sizeFactor` are dropped
        select_if(is.factor)

    # Drop any remaining factor columns that contain a single level
    hasLevels <- vapply(
        data,
        FUN = function(x) {
            length(levels(x)) > 1L
        },
        FUN.VALUE = logical(1L)
    )
    data <- data[, hasLevels, drop = FALSE]

    data <- column_to_rownames(data)

    if (ncol(data)) {
        data
    } else {
        warning("No valid annotation columns matched")
        NA
    }
}



# Define colors for each annotation column
.pheatmapAnnotationColors <- function(
    annotationCol,
    legendColor
) {
    if (is.data.frame(annotationCol) && is.function(legendColor)) {
        colors <- lapply(
            X = annotationCol,
            FUN = function(col) {
                assert_is_factor(col)
                levels <- levels(col)
                legendColor(length(levels))
            })
        names(colors) <- colnames(annotationCol)
        colors
    } else {
        NA
    }
}



# If `color = NULL`, use the pheatmap default palette
.pheatmapColor <- function(color = NULL, n = 256L) {
    if (!is.function(color)) {
        colorRampPalette(rev(
            brewer.pal(n = 7L, name = "RdYlBu")
        ))(n)
    } else {
        color(n)
    }
}
