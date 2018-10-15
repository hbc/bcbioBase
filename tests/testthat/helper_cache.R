invisible(lapply(
    X = c(
        "bcbio_nextgen.log",
        "data_versions.csv",
        "demultiplexed_invalid_duplicated.csv",
        "demultiplexed_invalid_legacy_samplename.csv",
        "demultiplexed_invalid_missing.csv",
        "demultiplexed_invalid_sample_id.csv",
        "demultiplexed.csv",
        "multiplexed_cellranger.csv",
        "multiplexed_indrops.csv",
        "multiplexed_invalid_duplicated.csv",
        "multiplexed_invalid_missing.csv",
        "programs.txt",
        "summary_invalid_metrics_mismatch.yaml",
        "summary_nested_metadata.yaml",
        "summary.yaml",
        "surecell_commands.log",
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
