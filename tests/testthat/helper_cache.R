invisible(lapply(
    X = c(
        "bcbio_legacy_samplename.csv",
        "bcbio_nextgen.log",
        "data_versions.csv",
        "demultiplexed.csv",
        "demultiplexed_duplicated_description.csv",
        "demultiplexed_missing_columns.csv",
        "multiplexed.csv",
        "multiplexed_cellranger.csv",
        "multiplexed_duplicated_sampleName.csv",
        "multiplexed_missing_columns.csv",
        "multiplexed_sampleID_column_defined.csv",
        "project_summary.yaml",
        "project_summary_metrics_mismatch.yaml",
        "project_summary_nested_metadata.yaml",
        "programs.txt",
        "tx2gene.csv"
    ),
    FUN = function(file, url) {
        if (!file.exists(file)) {
            utils::download.file(
                url = paste(url, file, sep = "/"),
                destfile = file
            )
        }
    },
    url = bcbioBaseCacheURL
))
