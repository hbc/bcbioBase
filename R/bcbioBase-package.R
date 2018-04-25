#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @importFrom Biostrings reverseComplement
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors as.data.frame cor mcols mcols<- metadata<-
#' @importFrom SummarizedExperiment assay assays colData colData<- rowData
#'   SummarizedExperiment
#' @importFrom basejump camel convertGenesToSymbols emptyRanges fixNA
#'   localOrRemoteFile makeNames readFileByExtension readYAML removeNA
#'   sanitizeSampleData snake
#' @importFrom dendsort dendsort
#' @importFrom dplyr arrange everything funs group_by mutate mutate_all
#'   mutate_at mutate_if select select_if ungroup
#' @importFrom grDevices colorRampPalette
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom methods .hasSlot as formalArgs getMethod is setAs slotNames
#'   validObject
#' @importFrom pheatmap pheatmap
#' @importFrom plyr ldply
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
#' @importFrom tidyr expand unite
#' @importFrom utils globalVariables sessionInfo
#' @importFrom viridis viridis
#'
#' @importFrom assertive assert_all_are_dirs
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom assertive assert_all_are_greater_than
#' @importFrom assertive assert_all_are_in_range
#' @importFrom assertive assert_all_are_matching_regex
#' @importFrom assertive assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive assert_all_are_non_negative
#' @importFrom assertive assert_all_are_positive
#' @importFrom assertive assert_are_disjoint_sets
#' @importFrom assertive assert_are_identical
#' @importFrom assertive assert_are_intersecting_sets
#' @importFrom assertive assert_has_colnames
#' @importFrom assertive assert_has_dimnames
#' @importFrom assertive assert_has_dims
#' @importFrom assertive assert_has_names
#' @importFrom assertive assert_has_rownames
#' @importFrom assertive assert_has_no_duplicates
#' @importFrom assertive assert_is_a_bool
#' @importFrom assertive assert_is_a_number
#' @importFrom assertive assert_is_a_string
#' @importFrom assertive assert_is_all_of
#' @importFrom assertive assert_is_an_integer
#' @importFrom assertive assert_is_any_of
#' @importFrom assertive assert_is_atomic
#' @importFrom assertive assert_is_character
#' @importFrom assertive assert_is_factor
#' @importFrom assertive assert_is_identical_to_na
#' @importFrom assertive assert_is_integer
#' @importFrom assertive assert_is_list
#' @importFrom assertive assert_is_matrix
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_subset
#' @importFrom assertive assert_is_tbl
#' @importFrom assertive has_dims
#' @importFrom assertive is_a_string
#'
#' @importFrom basejump assertHasRownames
#' @importFrom basejump assertIsAnImplicitInteger
#' @importFrom basejump assertIsAStringOrNULL
#' @importFrom basejump assertIsCharacterOrNULL
#' @importFrom basejump assertIsGene2symbol
#' @importFrom basejump assertIsHexColorFunctionOrNULL
#' @importFrom basejump assertIsTx2gene
"_PACKAGE"
