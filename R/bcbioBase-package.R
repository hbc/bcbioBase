#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @import S4Vectors methods
#'
#' @importFrom Biostrings reverseComplement
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom SummarizedExperiment colData colData<- rowData
#'   SummarizedExperiment
#' @importFrom basejump camel localOrRemoteFile makeNames readFileByExtension
#'   removeNA sanitizeSampleData
#' @importFrom dplyr arrange bind_rows group_by left_join mutate mutate_all
#'   mutate_if ungroup
#' @importFrom knitr kable
#' @importFrom magrittr set_rownames
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#' @importFrom readr read_csv read_lines
#' @importFrom rlang !!! !! abort inform sym syms warn
#' @importFrom scales percent
#' @importFrom sessioninfo session_info
#' @importFrom stringr str_pad str_trunc
#' @importFrom tibble as_tibble has_rownames is_tibble tibble
#' @importFrom tidyr unite
#' @importFrom utils globalVariables sessionInfo
"_PACKAGE"



globalVariables(".")
metadataPriorityCols <- c("sampleID", "sampleName", "description")



#' Project Directory Grep Pattern
#' @keywords internal
#' @export
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
