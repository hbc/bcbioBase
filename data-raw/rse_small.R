library(devtools)
library(pryr)
library(DESeq2)

# bcbioRNASeq coercion
# rse_small <- as(
#     object = bcbioRNASeq::bcb_small,
#     Class = "RangedSummarizedExperiment"
# )
# assays(rse_small) <- assays(rse_small)[1L]

# DESeqDataSet coercion
dds_small <- makeExampleDESeqDataSet()
rse_small <- as(dds_small, "RangedSummarizedExperiment")

stopifnot(identical(
    names(assays(rse_small)),
    "counts"
))
object_size(rse_small)
use_data(rse_small, compress = "xz", overwrite = TRUE)
