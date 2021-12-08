#' Create a CyTOF panel to be fill
#'
#'
#'
#' @param path Path to search for cytometry files.
#' @param glob A wildcard aka globbing pattern (e.g. *.fcs) passed on
#'  to [base::grep()] to filter paths.
#' @param flowcore_params Extra arguments passed to [flowCore::read.flowSet()].
#' @param cytofkit2_params Extra arguments passed to [cytofkit2::cytof_exprsMerge()].
#'
#' @return A tibble with 3 columns:
#'  1. **fcs_colname**: Column filled it with channels found in data.
#'  2. **antigen**: Filled with NAs.
#'  3. **marker_class**: Filled with NAs.
#' @export
create_panel <- function(path,
                         glob = "*.fcs",
                         flowcore_params = list(
                           truncate_max_range = FALSE
                         ),
                         cytofkit2_params = list(
                           mergeMethod = "all"
                         )) {

  ## Identify files to read ----
  fcs_files <- path |>
    fs::dir_ls(
      glob = glob,
      type = "file"
    )

  ## Extract channels names ----
  flowset_columns <- purrr::exec(
    .fn = flowCore::read.flowSet,
    files = fcs_files,
    !!!flowcore_params
  ) |>
    flowCore::colnames()

  cytofkit2_columns <- purrr::exec(
    .fn = cytofkit2::cytof_exprsMerge(
      fcsFiles = fcs_files,
      !!!cytofkit2_params
    )
  )

  extra_columns <- cytofkit2_columns |>
    stringr::str_remove("<.+>") |>
    setdiff(x = flowset_columns)

  ## Create panel table ----
  tibble::tibble(
    fcs_colname = c(cytofkit2_columns, extra_columns),
    antigen = NA,
    marker_class = NA
  )
}

#' Write CyTOF panel
#'
#'
#'
#' @param panel A data frame or data frame extension (e.g. a tibble).
#' @param file Path to the output file.
#' @inheritDotParams xlsx::write.xlsx -x -file -showNA -row.names -col.names -append
#'
#' @return NULL
#' @export
write_panel <- function(panel, file, ...) {
  xlsx::write.xlsx(
    x = as.data.frame(panel),
    file = file,
    showNA = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    append = FALSE,
    ...
  )
}
