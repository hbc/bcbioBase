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
            "sample_4")))
# Check handling of rowData (annotable) mismatch
rowdata <- data.frame(
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
        "ENSMUSG00000000031"))
coldata <- data.frame(
    genotype = c(
        "wildtype",
        "wildtype",
        "knockout",
        "knockout"),
    age = c(3L, 6L, 3L, 6L),
    row.names = colnames(mat))
se <- prepareSummarizedExperiment(
    assays = list(assay = mat),
    rowData = rowdata,
    colData = coldata)

test_that("Valid SummarizedExperiment", {
    expect_identical(
        dim(se),
        c(4L, 4L))
    expect_identical(
        names(metadata(se)),
        c(
            "date",
            "wd",
            "utilsSessionInfo",
            "devtoolsSessionInfo"
        )
    )
})

test_that("Empty row and/or column data", {
    noanno <- suppressWarnings(
        prepareSummarizedExperiment(assays = list(assay = mat))
    )
    expect_warning(
        prepareSummarizedExperiment(
        assays = list(assay = mat),
        colData = coldata),
        paste(
            "Summarizing experiment without row data",
            "\\(e.g. gene annotations\\)"
        )
    )
    expect_warning(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowdata),
        paste(
            "Summarizing experiment without column data",
            "\\(e.g. sample metadata\\)"
        )
    )
    expect_identical(length(slot(noanno, "elementMetadata")), 0L)
    expect_identical(length(slot(noanno, "colData")), 0L)
})

test_that("Ensure `assays()` requires a list", {
    expect_error(
        prepareSummarizedExperiment(
            assays = mat,
            rowData = rowdata,
            colData = coldata),
        paste(
            "unable to find an inherited method for function",
            "'prepareSummarizedExperiment' for signature '\"matrix\"'"
        )
    )
})

# This checks to see if there are any dashes in the names
test_that("Enforce strict names", {
    matbadrows <- mat
    rownames(matbadrows) <- paste0(rownames(matbadrows), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matbadrows),
            rowData = rowdata,
            colData = coldata),
        paste(
            "are_identical :",
            "make.names\\(rownames\\(assay\\), unique = TRUE, allow_ = TRUE\\)",
            "and rownames\\(assay\\) are not identical."
        )
    )
    matbadcols <- mat
    colnames(matbadcols) <- paste0(colnames(matbadcols), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matbadcols),
            rowData = rowdata,
            colData = coldata),
        paste(
            "are_identical :",
            "make.names\\(colnames\\(assay\\), unique = TRUE, allow_ = TRUE\\)",
            "and colnames\\(assay\\) are not identical."
        )
    )
})

test_that("Duplicate names", {
    matduperows <- mat
    rownames(matduperows) <- c(
        "ENSMUSG00000000001",
        "ENSMUSG00000000001",
        "ENSMUSG00000000003",
        "ENSMUSG00000000003")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matduperows),
            rowData = rowdata,
            colData = coldata),
        paste(
            "has_no_duplicates :",
            "rownames\\(assay\\) has duplicates at positions 2, 4."
        )
    )
    matdupecols <- mat
    colnames(matdupecols) <- c(
        "sample_1",
        "sample_1",
        "sample_2",
        "sample_2")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matdupecols),
            rowData = rowdata,
            colData = coldata),
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
            rowData = rowdata,
            colData = coldata),
        "has_dimnames : The dimension names of assay are NULL."
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowdata,
            colData = c(xxx = "yyy")
        ),
        "is2 : colData"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = c(xxx = "yyy"),
            colData = coldata
        ),
        "is2 : rowData"
    )
})

test_that("Invalid metadata", {
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowdata,
            colData = coldata,
            metadata = Sys.Date()
        ),
        "is2 : metadata")
})

test_that("Dimension mismatch", {
    matextracol <- cbind(mat, "sample_5" = seq(17L, 20L))
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matextracol),
            rowData = rowdata,
            colData = coldata),
        paste(
            "are_identical :",
            "colnames\\(assay\\) and rownames\\(colData\\) are not identical."
        )
    )
    matextrarow <- rbind(mat, "ENSMUSG00000000037" = seq(17L, 20L))
    expect_warning(
        prepareSummarizedExperiment(
            assays = list(assay = matextrarow),
            rowData = rowdata,
            colData = coldata),
        "Unannotated genes detected in assay \\(20%\\)")
    seextrarow <- suppressWarnings(
        prepareSummarizedExperiment(
            assays = list(assay = matextrarow),
            rowData = rowdata,
            colData = coldata)
    )
    expect_identical(
        metadata(seextrarow)[["unannotatedGenes"]],
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
            rowData = rowdata,
            colData = coldata),
        "has_rownames : The row names of assay are NULL."
    )
    matnocolnames <- mat
    colnames(matnocolnames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matnocolnames),
            rowData = rowdata,
            colData = coldata),
        "has_colnames : The column names of assay are NULL."
    )
    rowdatanorownames <- rowdata
    rownames(rowdatanorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowdatanorownames,
            colData = coldata),
        paste(
            "are_intersecting_sets :",
            "rownames\\(assay\\) and rownames\\(rowData\\)",
            "have no common elements."
        )
    )
    coldatanorownames <- coldata
    rownames(coldatanorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowdata,
            colData = coldatanorownames),
        paste(
            "are_identical :",
            "colnames\\(assay\\) and rownames\\(colData\\) are not identical."
        )
    )
})
