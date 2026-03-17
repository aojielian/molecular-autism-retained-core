# scoring_utils.R
# Helper functions for donor-level aggregation and module-score handling

aggregate_donor_scores <- function(meta_df, donor_col, group_col, score_cols) {
  stopifnot(all(c(donor_col, group_col, score_cols) %in% colnames(meta_df)))
  dplyr::as_tibble(meta_df) %>%
    dplyr::mutate(Donor = .data[[donor_col]], Group = .data[[group_col]]) %>%
    dplyr::group_by(Donor, Group) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(score_cols), ~mean(.x, na.rm = TRUE)),
      n_cells = dplyr::n(),
      .groups = "drop"
    )
}
