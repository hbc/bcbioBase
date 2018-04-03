library(devtools)
library(pryr)
library(bcbioRNASeq)
rse_small <- as(
    object = bcbioRNASeq::bcb_small,
    Class = "RangedSummarizedExperiment"
)
assays(rse_small) <- assays(rse_small)[1L]
object_size(rse_small)
use_data(rse_small, compress = "xz", overwrite = TRUE)
