genes <- paste0("gene", seq_len(4L))
samples <- paste0("sample", seq_len(4L))

mat <- matrix(
    data = seq(1L:16L),
    nrow = 4L,
    ncol = 4L,
    byrow = FALSE,
    dimnames = list(genes, samples)
)
dgc <- as(mat, "dgCMatrix")

heatmapList <- c("tree_row", "tree_col", "kmeans", "gtable")

rr <- GRanges(
    seqnames = replicate(n = 4L, expr = "1"),
    ranges = IRanges(
        start = c(1L, 101L, 201L, 301L),
        end = c(100L, 200L, 300L, 400L)
    )
)
names(rr) <- genes

cd <- data.frame(
    "genotype" = c(
        "wildtype",
        "wildtype",
        "knockout",
        "knockout"
    ),
    "age" = c(3L, 6L, 3L, 6L),
    row.names = samples
)
