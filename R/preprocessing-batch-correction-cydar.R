#' Split matrix by batch and condition
#'
#' @param mat
#' @param metadata_sep
#' @param batch_position
#' @param condition_position
#' @param control_condition_string
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
split_by_batch_condition <- function(mat,
                                     metadata_sep = "_",
                                     batch_position = 1,
                                     condition_position = 2,
                                     control_condition_string = "Sp",
                                     ...) {
  samples_rownames <- rownames(mat)

  ## Prepare data
  samples_rownames |>
    stringr::str_remove(pattern = "_\\d+$") |>
    unique() |>
    (\(x) {
      tibble::tibble(
        sample = x
      )
    })() |>
    dplyr::mutate(
      metadata_info = stringr::str_split(
        string = sample,
        pattern = metadata_sep,
        simplify = TRUE
      ),
      batch = metadata_info[, batch_position],
      batch = stringr::str_c(batch, metadata_sep),
      condition = metadata_info[, condition_position],
    ) |>
    dplyr::mutate(
      condition = dplyr::case_when(
        stringr::str_detect(
          string = condition,
          pattern = control_condition_string
        ) ~ "Control",
        TRUE ~ "Sample"
      ),
      condition = factor(
        x = condition,
        levels = c("Control", "Sample")
      ),
      sample_data = purrr::map(
        .x = sample,
        .f = ~ mat[grep(.x, samples_rownames), ]
      )
    ) |>
    dplyr::arrange(batch, dplyr::desc(condition)) |>
    dplyr::select(
      -metadata_info
    ) |>
    identity()
}

#' Prepare data to cydar batch correction format
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
prepare_batch_correction_format <- function(data) {
  data |>
    dplyr::arrange(dplyr::desc(condition)) |>
    dplyr::group_by(batch) |>
    dplyr::summarize(
      batch_x = list(sample_data),
      batch_x = purrr::map(
        .x = batch_x,
        .f = ~ {
          purrr::set_names(x = .x, nm = sample) |>
            list() |>
            purrr::set_names(nm = dplyr::cur_group())
        }
      ),
      batch_comp = list(condition),
      batch_comp = purrr::map(
        .x = batch_comp,
        .f = ~ {
          purrr::set_names(x = .x, nm = condition) |>
            list() |>
            purrr::set_names(nm = dplyr::cur_group())
        }
      ),
      .groups = "drop"
    )
}


#' Run cydar batch correction
#'
#' @param batch_correction_df
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
run_batch_correction <- function(batch_correction_df,
                                 ...) {
  batch_x <- batch_correction_df |>
    dplyr::pull(batch_x) |>
    unlist(recursive = FALSE)

  batch_comp <- batch_correction_df |>
    dplyr::pull(batch_comp) |>
    unlist(recursive = FALSE)

  batch_corrected_results <- cydar::normalizeBatch(
    batch.x = batch_x,
    batch.comp = batch_comp,
    ...
  )

  batch_correction_df |>
    dplyr::transmute(
      batch,
      sample = purrr::map_chr(
        .x = batch_corrected_results,
        .f = ~ names(.x)[1]
      ),
      data_x = purrr::map(
        .x = batch_x,
        .f = ~ .x[[1]][[1]]
      ),
      data_x = purrr::set_names(
        x = data_x,
        nm = sample
      ),
      data_corrected = purrr::map(
        .x = batch_corrected_results,
        .f = dplyr::first
      ),
      data_corrected = purrr::set_names(
        x = data_corrected,
        nm = sample
      )
    )
}


#' Title
#'
#' @param batch_corrected_df
#'
#' @return
#' @export
#'
#' @examples
batch_corrected_drop_na <- function(batch_corrected_df) {
  batch_corrected_df |>
    dplyr::rowwise(batch, sample) |>
    dplyr::summarize(
      data_x = list(data_x) |>
        purrr::set_names(
          nm = dplyr::cur_group()[["sample"]]
        ),
      data_corrected = data_corrected |>
        (\(x){
          na_channels = apply(X = x, MARGIN = 2, FUN = anyNA)
          data_corrected[, !na_channels]
        })() |>
        list() |>
        purrr::set_names(
          nm = dplyr::cur_group()[["sample"]]
        ),
      .groups = "drop"
    )
}


#' Save cydar batch corrected results
#'
#' @param batch_corrected_df
#' @param output_dir
#' @param prefix
#' @param suffix
#'
#' @return
#' @export
#'
#' @examples
save_batch_corrected_results <- function(batch_corrected_df,
                                         output_dir,
                                         prefix = "",
                                         suffix = "") {
  fs::dir_create(output_dir)

  batch_corrected_df |>
    dplyr::mutate(
      output_file = fs::path(
        output_dir,
        glue::glue("{prefix}{sample}{suffix}.fcs")
      )
    ) |>
    dplyr::rowwise(output_file) |>
    dplyr::summarize(
      new(
        Class = "flowFrame",
        exprs = as.matrix(data_corrected)
      ) |>
        flowCore::write.FCS(
          filename = output_file
        ),
      .groups = "drop"
    )

  invisible(batch_corrected_df)
}

#' Plot before and after cydar batch correction
#'
#' @param batch_corrected_df
#' @param output_file
#'
#' @return
#' @export
#'
#' @examples
plot_batch_corrected_results <- function(batch_corrected_df) {

}
