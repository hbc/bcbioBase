library(devtools)
library(pryr)
library(bcbioRNASeq)

bcb <- bcbioRNASeq::bcb_small
rse_bcb <- as(bcb, "RangedSummarizedExperiment")
assays(rse_bcb) <- assays(rse_bcb)[1L]
object_size(rse_bcb)

use_data(rse_bcb, compress = "xz", overwrite = TRUE)
