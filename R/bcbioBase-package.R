#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @keywords internal
#'
#' @importFrom Biostrings reverseComplement
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors as.data.frame cor mcols mcols<- metadata<-
#' @importFrom SummarizedExperiment assay colData colData<- rowData
#'   SummarizedExperiment
#' @importFrom basejump camel convertGenesToSymbols emptyRanges
#'   localOrRemoteFile makeNames readFileByExtension removeNA sanitizeSampleData
#'   snake
#' @importFrom dendsort dendsort
#' @importFrom dplyr arrange bind_rows group_by left_join mutate mutate_all
#'   mutate_if select_if ungroup
#' @importFrom grDevices colorRampPalette
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom methods .hasSlot as formalArgs getMethod is slotNames validObject
#' @importFrom pheatmap pheatmap
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#' @importFrom readr read_csv read_lines
#' @importFrom rlang !!! !! sym syms
#' @importFrom scales percent
#' @importFrom sessioninfo session_info
#' @importFrom stats hclust quantile
#' @importFrom stringr str_pad str_trunc
#' @importFrom tibble as_tibble column_to_rownames has_rownames is_tibble
#'   rownames_to_column tibble
#' @importFrom tidyr unite
#' @importFrom utils globalVariables sessionInfo
#' @importFrom viridis viridis
"_PACKAGE"
