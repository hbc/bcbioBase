#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importFrom assertive.base assert_are_identical assert_is_identical_to_na
#' @importFrom assertive.data assert_all_are_hex_colors
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_greater_than
#'   assert_all_are_in_range assert_all_are_non_negative assert_all_are_positive
#' @importFrom assertive.properties assert_has_colnames assert_has_dimnames
#'   assert_has_dims assert_has_names assert_has_no_duplicates
#'   assert_has_rownames assert_is_atomic assert_is_non_empty has_dims
#' @importFrom assertive.sets assert_are_disjoint_sets
#'   assert_are_intersecting_sets assert_is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#'   assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive.types assert_is_a_bool assert_is_a_number
#'   assert_is_a_string assert_is_all_of assert_is_an_integer assert_is_any_of
#'   assert_is_character assert_is_factor assert_is_function assert_is_integer
#'   assert_is_list assert_is_matrix assert_is_tbl is_a_string
#' @importFrom basejump assertFormalInterestingGroups assertHasRownames
#'   assertIsAnImplicitInteger assertIsAStringOrNULL
#'   assertIsHexColorFunctionOrNULL assertIsTx2gene camel emptyRanges fixNA
#'   localOrRemoteFile makeNames printString readFileByExtension readYAML
#'   removeNA snake
#' @importFrom Biostrings reverseComplement
#' @importFrom dendsort dendsort
#' @importFrom dplyr arrange everything funs group_by left_join mutate
#'   mutate_all mutate_at mutate_if select select_if ungroup
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom ggplot2 aes geom_hline geom_label geom_vline
#' @importFrom ggrepel geom_label_repel
#' @importFrom grDevices colorRampPalette
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom methods .hasSlot as formalArgs getMethod is setAs slotNames
#'   validObject
#' @importFrom pheatmap pheatmap
#' @importFrom plyr ldply
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#' @importFrom readr read_csv read_lines
#' @importFrom rlang !!! !! sym syms
#' @importFrom S4Vectors aggregate as.data.frame cor mcols mcols<- merge
#'   metadata<-
#' @importFrom scales percent
#' @importFrom sessioninfo session_info
#' @importFrom stats as.formula hclust quantile
#' @importFrom stringr str_pad str_trunc
#' @importFrom SummarizedExperiment assay assays colData colData<- rowData
#'   SummarizedExperiment
#' @importFrom tibble as_tibble column_to_rownames has_rownames is_tibble
#'   rownames_to_column tibble
#' @importFrom tidyr expand
#' @importFrom utils globalVariables sessionInfo
#' @importFrom viridis viridis
"_PACKAGE"
