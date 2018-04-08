library(devtools)
library(pryr)
library(bcbioRNASeq)
library(DESeq2)

# DESeqDataSet coercion
dds <- makeExampleDESeqDataSet()
rse_dds <- as(dds, "RangedSummarizedExperiment")
object_size(rse_dds)

# bcbioRNASeq coercion
bcb <- bcbioRNASeq::bcb_small
rse_bcb <- as(bcb, "RangedSummarizedExperiment")
assays(rse_bcb) <- assays(rse_bcb)[1L]
object_size(rse_bcb)

use_data(rse_dds, rse_bcb, compress = "xz", overwrite = TRUE)
