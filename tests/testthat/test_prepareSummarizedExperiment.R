context("prepareSummarizedExperiment")

# Create the matrix with invalid names. We'll sanitize these
# into snake_case later.
mat <- matrix(
    seq(1L:16L),
    nrow = 4L,
    ncol = 4L,
    dimnames = list(
        c(
            "ENSMUSG00000000001",
            "ENSMUSG00000000003",
            "ENSMUSG00000000028",
            "ENSMUSG00000000031"),
        c(
            "sample_1",
            "sample_2",
            "sample_3",
            "sample_4")
        )
)
# Check handling of rowData (annotable) mismatch
rowData <- data.frame(
    ensgene = c(
        "ENSMUSG00000000001",
        "ENSMUSG00000000003",
        "ENSMUSG00000000028",
        "ENSMUSG00000000031"),
    biotype = c(
        "coding",
        "coding",
        "coding",
        "coding"),
    row.names = c(
        "ENSMUSG00000000001",
        "ENSMUSG00000000003",
        "ENSMUSG00000000028",
        "ENSMUSG00000000031")
)
colData <- data.frame(
    genotype = c(
        "wildtype",
        "wildtype",
        "knockout",
        "knockout"),
    age = c(3L, 6L, 3L, 6L),
    row.names = colnames(mat)
)


test_that("SummarizedExperiment", {
    se <- prepareSummarizedExperiment(
        assays = list(assay = mat),
        rowData = rowData,
        colData = colData)
    expect_s4_class(se, "SummarizedExperiment")
    expect_identical(dim(se), c(4L, 4L))
    expect_identical(
        lapply(metadata(se), class),
        list(
            "date" = "Date",
            "wd" = c("fs_path", "character"),
            "utilsSessionInfo" = "sessionInfo",
            "devtoolsSessionInfo" = "session_info"
        )
    )
})

test_that("RangedSummarizedExperiment", {
    rowData <- genes("Mus musculus")
    expect_s4_class(rowData, "GRanges")
    se <- prepareSummarizedExperiment(
        assays = list(assay = mat),
        rowData = rowData,
        colData = colData)
    expect_s4_class(se, "RangedSummarizedExperiment")
})

test_that("Empty row and/or column data", {
    noAnno <- suppressWarnings(
        prepareSummarizedExperiment(assays = list(assay = mat))
    )
    expect_warning(
        prepareSummarizedExperiment(
        assays = list(assay = mat),
        colData = colData),
        "Summarizing experiment without row data"
    )
    expect_warning(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowData),
        "Summarizing experiment without column data"
    )
    expect_identical(length(slot(noAnno, "elementMetadata")), 0L)
    expect_identical(length(slot(noAnno, "colData")), 0L)
})

# This checks to see if there are any dashes in the names
test_that("Enforce strict names", {
    matBadRows <- mat
    rownames(matBadRows) <- paste0(rownames(matBadRows), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matBadRows),
            rowData = rowData,
            colData = colData),
        paste(
            "are_identical :",
            "make.names\\(rownames\\(assay\\), unique = TRUE, allow_ = TRUE\\)",
            "and rownames\\(assay\\) are not identical."
        )
    )
    matBadCols <- mat
    colnames(matBadCols) <- paste0(colnames(matBadCols), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matBadCols),
            rowData = rowData,
            colData = colData),
        paste(
            "are_identical :",
            "make.names\\(colnames\\(assay\\), unique = TRUE, allow_ = TRUE\\)",
            "and colnames\\(assay\\) are not identical."
        )
    )
})

test_that("Duplicate names", {
    matDupeRows <- mat
    rownames(matDupeRows) <- c(
        "ENSMUSG00000000001",
        "ENSMUSG00000000001",
        "ENSMUSG00000000003",
        "ENSMUSG00000000003")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matDupeRows),
            rowData = rowData,
            colData = colData),
        paste(
            "has_no_duplicates :",
            "rownames\\(assay\\) has duplicates at positions 2, 4."
        )
    )
    matDupeCols <- mat
    colnames(matDupeCols) <- c(
        "sample_1",
        "sample_1",
        "sample_2",
        "sample_2")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matDupeCols),
            rowData = rowData,
            colData = colData),
        paste(
            "has_no_duplicates :",
            "colnames\\(assay\\) has duplicates at positions 2, 4."
        )
    )
})

test_that("Column data pass-in failure in assays", {
    # Bad pass-in of objects not supporting `dimnames()`
    expect_error(
        prepareSummarizedExperiment(
            assays = list(c(xxx = "yyy")),
            rowData = rowData,
            colData = colData),
        "has_dimnames : The dimension names of assay are NULL."
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowData,
            colData = c(xxx = "yyy")
        ),
        "is2 : colData"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = c(xxx = "yyy"),
            colData = colData
        ),
        "is2 : rowData"
    )
})

test_that("Invalid metadata", {
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowData,
            colData = colData,
            metadata = Sys.Date()
        ),
        "is2 : metadata")
})

test_that("Dimension mismatch", {
    matExtraCol <- cbind(mat, "sample_5" = seq(17L, 20L))
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matExtraCol),
            rowData = rowData,
            colData = colData),
        paste(
            "are_identical :",
            "colnames\\(assay\\) and rownames\\(colData\\) are not identical."
        )
    )
    matExtraRow <- rbind(mat, "ENSMUSG00000000037" = seq(17L, 20L))
    expect_warning(
        prepareSummarizedExperiment(
            assays = list(assay = matExtraRow),
            rowData = rowData,
            colData = colData),
        "1 unannotated rows detected")
    seExtraRow <- suppressWarnings(
        prepareSummarizedExperiment(
            assays = list(assay = matExtraRow),
            rowData = rowData,
            colData = colData)
    )
    expect_identical(
        metadata(seExtraRow)[["unannotatedRows"]],
        "ENSMUSG00000000037"
    )
})

test_that("Missing rownames", {
    # Missing rownames
    matnorownames <- mat
    rownames(matnorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matnorownames),
            rowData = rowData,
            colData = colData),
        "has_rownames : The row names of assay are NULL."
    )
    matnocolnames <- mat
    colnames(matnocolnames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matnocolnames),
            rowData = rowData,
            colData = colData),
        "has_colnames : The column names of assay are NULL."
    )
    rowdatanorownames <- rowData
    rownames(rowdatanorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowdatanorownames,
            colData = colData),
        paste(
            "are_intersecting_sets :",
            "rownames\\(assay\\) and rownames\\(rowData\\)",
            "have no common elements."
        )
    )
    coldatanorownames <- colData
    rownames(coldatanorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowData,
            colData = coldatanorownames),
        paste(
            "are_identical :",
            "colnames\\(assay\\) and rownames\\(colData\\) are not identical."
        )
    )
})
