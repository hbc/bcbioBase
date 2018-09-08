cacheURL <- "http://bcbiobase.seq.cloud"
files <- c(
    "bcbio_legacy_samplename.csv",
    "bcbio_nextgen.log",
    "cellranger_metadata.csv",
    "data_versions.csv",
    "demultiplexed.csv",
    "demultiplexed_duplicated_description.csv",
    "demultiplexed_missing_columns.csv",
    "multiplexed.csv",
    "multiplexed_duplicated_sampleName.csv",
    "multiplexed_missing_columns.csv",
    "multiplexed_sampleID_column_defined.csv",
    "project_summary.yaml",
    "project_summary_metrics_mismatch.yaml",
    "project_summary_nested_metadata.yaml",
    "programs.txt",
    "tx2gene.csv"
)
mapply(
    FUN = function(cacheURL, file, envir) {
        if (!file.exists(file)) {
            utils::download.file(
                url = paste(cacheURL, file, sep = "/"),
                destfile = file)
        }
    },
    file = files,
    MoreArgs = list(cacheURL = cacheURL, envir = environment())
)
