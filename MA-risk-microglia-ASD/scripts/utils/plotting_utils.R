# plotting_utils.R
# Common plotting helpers for Molecular Autism revision project

theme_ma <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
}
