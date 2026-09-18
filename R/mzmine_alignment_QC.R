#' mzmine_to_se
#'
#' Gets the alignment scores from mzmine feature table
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% mutate group_by
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @param path_to_file path to the mzmine feature table
#' @param rowData_cols columns to be used as rowData
#'
#' @returns a data.frame wit alignemtn_scores
#' @export
#'

mzmine_alignment_qc <- function(
    path_to_file,
    rowData_cols = c("id","rt","mz")
) {

  # load feature data

  print("Loading mzMine feature table")

  features <- read.csv(
    path_to_file
  ) %>%
    mutate(
      id = as.character(id) # convert id from int to character
    )

  # extract columns with height or area

  print("Processing feature table")

  alignment_table <- features[
    ,c(
      rowData_cols,
      names(features)[grepl(paste0("alignment_scores[.]"),names(features))]
    )
  ]
  names(alignment_table) <- gsub(
    "alignment_scores[.]", "", names(alignment_table)
  )

  return(alignment_table)

}

mzmine_alignment_QC <- function(...) mzmine_alignment_qc(...)


