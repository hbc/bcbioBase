# Sanitize formals into snake case and abort on duplicates.
# Duplicates may arise if user is mixing and matching camel/snake case.
.pheatmapArgs <- function(args) {
    assert_is_list(args)
    assert_has_names(args)
    # Abort on snake case formatted formalArgs
    invalidNames <- grep("[._]", names(args), value = TRUE)
    if (length(invalidNames)) {
        abort(paste(
            "Define formalArgs in camel case:",
            toString(invalidNames)
        ))
    }
    names(args) <- snake(names(args))
    assert_is_subset(names(args), formalArgs(pheatmap))
    args
}



# Automatically handle the annotation columns
.pheatmapAnnotationCol <- function(object) {
    # pheatmap requires `NA` argument if empty
    if (!has_dims(object)) {
        return(NA)
    }
    assertHasRownames(object)
    object %>%
        as.data.frame() %>%
        # Remove sample name columns
        .[, setdiff(colnames(.), metadataPriorityCols), drop = FALSE] %>%
        rownames_to_column() %>%
        # Ensure all strings are factor
        mutate_if(is.character, as.factor) %>%
        # Ensure unwanted columns like `sizeFactor` are dropped
        select_if(is.factor) %>%
        column_to_rownames()
}



# Define colors for each annotation column
.pheatmapAnnotationColors <- function(
    annotationCol,
    legendColor
) {
    if (is.data.frame(annotationCol) && is.function(legendColor)) {
        # FIXME Switch to mapply
        lapply(
            seq_along(colnames(annotationCol)), function(a) {
                col <- levels(annotationCol[[a]])
                colors <- annotationCol[[a]] %>%
                    levels() %>%
                    length() %>%
                    legendColor()
                names(colors) <- col
                colors
            }) %>%
            set_names(colnames(annotationCol))
    } else {
        NA
    }
}



# If `color = NULL`, use the pheatmap default palette
.pheatmapColor <- function(color) {
    nColor <- 256L
    if (!is.function(color)) {
        color <- colorRampPalette(rev(
            brewer.pal(n = 7L, name = "RdYlBu")
        ))(nColor)
    } else {
        color <- color(nColor)
    }
}
