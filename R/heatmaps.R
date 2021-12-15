#' Create three avera expression heatmaps using a SCE object
#'
#' This function is a wrapper around [CATALYST::plotExprHeatmap()] and
#' [ComplexHeatmap::Heatmap()] that allows more styling flexibility.
#'
#' @param sce SingleCellExperiment object.
#' @param by Data to use in addition to markers.
#' @param by_keep Elements to subset before plotting (no features).
#' @param by_label Title to accompany annotation.
#' @param by_colors Color to use in annotation.
#' @param catalyst_args Extra arguments for [CATALYST::plotExprHeatmap()].
#'   Except x, by, and scale.
#' @param heatmap_col A vector of colors to be interpolated.
#' @param rescale_limits Output range (numeric vector of length two).
#' @param special_text Logical. If TRUE [ComplexHeatmap::Heatmap()] uses
#'   [ComplexHeatmap::gt_render()] for rownames and colnames of heatmap.
#' @inheritDotParams ComplexHeatmap::Heatmap -name -matrix -col
#'
#' @return A named list with three [ComplexHeatmap::Heatmap()] objects:
#'   1. normalized_heatmap: Use resulting matrix from
#'      [CATALYST::plotExprHeatmap()] without modification.
#'   2. scaled_heatmap: Matrix scaled without center and then rescaled using
#'      [scales::rescale()].
#'   3. scaled_markerwise_heatmap: Matrix scaled without center and then
#'      rescaled by each marker using [scales::rescale()].
#'
#' @examples
#' @md
catalyst_expression_heatmaps <- function(sce,
                                         by = "sample_id",
                                         by_keep = NULL,
                                         by_label = by,
                                         by_colors = NULL,
                                         catalyst_args = list(),
                                         heatmap_col = NULL,
                                         rescale_limits = c(0, 1),
                                         special_text = FALSE,
                                         ...) {
  ## Calculate matrix to plot ----
  mat <- purrr::exec(
    .fn = CATALYST::plotExprHeatmap,
    x = sce,
    by = by,
    scale = "never",
    !!!catalyst_args
  ) |>
    purrr::pluck("matrix")

  ## Get extra parameters to handle.
  extra_heatmap_params <- list(...)

  ## Find m or k to by_keep ----
  if (is.null(by_keep)) {
    by_keep <- rownames(mat)
  }

  ## Set plot defaults ----
  ### Column annotation pallete ----
  preexisting_top_annotation <- "top_annotation" %in% names(extra_heatmap_params)

  if (!preexisting_top_annotation) {
    if (is.null(by_colors)) {
      by_colors <- CATALYST:::.cluster_cols[1:nrow(mat)]
      names(by_colors) <- rownames(mat)
      by_colors <- by_colors[by_keep]
    }

    top_annotation <- ComplexHeatmap::columnAnnotation(
      col_annotation = ComplexHeatmap::anno_simple(
        x = by_keep,
        col = by_colors
      ),
      annotation_label = by_label
    )

    extra_heatmap_params[["top_annotation"]] <- top_annotation
  }

  ### Heatmap continious palette ----
  custom_color_palette <- function(colors, matrix) {
    circlize::colorRamp2(
      breaks = seq(
        from = min(matrix, na.rm = TRUE),
        to = max(matrix, na.rm = TRUE),
        length.out = length(colors)
      ),
      colors = colors,
      space = "LAB"
    )
  }

  if (is.null(heatmap_col)) {
    heatmap_col <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  }

  ### Function annotation ----
  fun <- purrr::pluck(
    .x = catalyst_args,
    "fun",
    .default = "median"
  )

  ### Labels annotation ----
  if (special_text) {
    if (!"column_labels" %in% names(extra_heatmap_params)) {
      extra_heatmap_params[["column_labels"]] <- mat[by_keep, ] |>
        t() |>
        colnames() |>
        ComplexHeatmap::gt_render()
    }

    if (!"row_labels" %in% names(extra_heatmap_params)) {
      extra_heatmap_params[["row_labels"]] <- mat[by_keep, ] |>
        t() |>
        rownames() |>
        ComplexHeatmap::gt_render()
    }
  }

  ## Create shared function for plotting ----
  custom_heatmap <- function(matrix,
                             name) {
    purrr::exec(
      .fn = ComplexHeatmap::Heatmap,
      name = name,
      matrix = matrix,
      col = custom_color_palette(heatmap_col, matrix),
      !!!extra_heatmap_params
    )
  }

  ## Plot matrices ----
  normalized_heatmap <- mat[by_keep, ] |>
    t() |>
    (\(x) {
      custom_heatmap(
        name = glue::glue("{stringr::str_to_title(fun)}\nexpression"),
        matrix = x
      )
    })()

  scaled_heatmap <- mat |>
    scale(center = FALSE) |>
    scales::rescale(to = rescale_limits) |>
    (\(x) {
      x[by_keep, ]
    })() |>
    t() |>
    (\(x) {
      custom_heatmap(
        name = glue::glue("Scaled\n{fun}\nexpression"),
        matrix = x
      )
    })()

  scaled_marker_wise <- mat |>
    scale(center = FALSE) |>
    asplit(MARGIN = 2) |>
    lapply(function(x) scales::rescale(x = x, to = rescale_limits)) |>
    as.data.frame() |>
    as.matrix() |>
    (\(x) {
      x[by_keep, ]
    })() |>
    t() |>
    (\(x){
      custom_heatmap(
        name = glue::glue("Scaled\n{fun}\nexpression\nby marker"),
        matrix = x
      )
    })()

  list(
    normalized_heatmap = normalized_heatmap,
    scaled_heatmap = scaled_heatmap,
    scaled_marker_wise = scaled_marker_wise
  )
}
