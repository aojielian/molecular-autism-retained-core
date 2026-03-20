#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(magick)
})

# ------------------------------------------------------------------------------
# Supplementary_Figure_S7.R
# Build Supplementary Figure S7 from pre-rendered panel images:
#   A: all-cluster UMAP
#   B: canonical marker dot plot
#   C: donor-level all-cluster abundance forest
#   D: resolution benchmark overview
#
# Output:
#   Supplementary_Figure_S7.png
#   Supplementary_Figure_S7.pdf
#
# Example:
#   Rscript Supplementary_Figure_S7.R \
#     --root /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision \
#     --outdir /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Supplementary_Figure_S7 \
#     --width 5200 \
#     --height 3600
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) if (!is.null(x)) x else y

parse_args_simple <- function() {
  x <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

args <- parse_args_simple()

root   <- args$root   %||% stop("Missing required argument: --root", call. = FALSE)
outdir <- args$outdir %||% stop("Missing required argument: --outdir", call. = FALSE)
width  <- as.integer(args$width  %||% "5200")
height <- as.integer(args$height %||% "3600")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pkg1  <- file.path(root, "Package1_Figure1_DiscoveryAudit_reviewer_plus")
pkg1c <- file.path(root, "Package1c_resolution_sensitivity_benchmark")

panel_paths <- list(
  A = file.path(pkg1,  "01_microglia_allcluster_umap.png"),
  B = file.path(pkg1,  "20_canonical_allcluster_marker_dotplot.png"),
  C = file.path(pkg1,  "24_allcluster_abundance_reviewer_forest.png"),
  D = file.path(pkg1c, "06_resolution_benchmark_overview.png")
)

missing <- names(panel_paths)[!file.exists(unlist(panel_paths))]
if (length(missing) > 0) {
  stop(
    "Missing input panel files:\n",
    paste(sprintf("  %s: %s", missing, unlist(panel_paths)[missing]), collapse = "\n"),
    call. = FALSE
  )
}

message("[INFO] root   = ", root)
message("[INFO] outdir = ", outdir)
for (nm in names(panel_paths)) {
  message("[INFO] panel ", nm, " = ", panel_paths[[nm]])
}

trim_and_prepare <- function(path, target_w, target_h, fuzz = 12) {
  img <- image_read(path)
  img <- image_background(img, "white", flatten = TRUE)
  img <- image_trim(img, fuzz = fuzz)
  info <- image_info(img)
  scale <- min(target_w / info$width, target_h / info$height)
  new_w <- max(1L, floor(info$width  * scale))
  new_h <- max(1L, floor(info$height * scale))
  img <- image_resize(img, sprintf("%dx%d!", new_w, new_h))
  img
}

make_panel <- function(img, label, panel_w, panel_h,
                       label_band_w = 52,
                       label_band_h = 84,
                       pad_r = 24,
                       pad_b = 24,
                       label_size = 72,
                       img_x_extra = 0,
                       img_y_extra = 0) {
  panel <- image_blank(width = panel_w, height = panel_h, color = "white")
  avail_w <- panel_w - label_band_w - pad_r - img_x_extra
  avail_h <- panel_h - label_band_h - pad_b - img_y_extra
  img <- trim_and_prepare(img, target_w = avail_w, target_h = avail_h)

  x_off <- label_band_w + img_x_extra
  y_off <- label_band_h + img_y_extra
  panel <- image_composite(panel, img, offset = sprintf("+%d+%d", x_off, y_off))

  # Put label in a dedicated gutter so it never overlaps the plot body.
  panel <- image_annotate(
    panel, text = label, size = label_size, weight = 700,
    gravity = "northwest", location = "+8+0", color = "black"
  )
  panel
}

# 2 x 2 layout
left   <- 60L
top    <- 50L
right  <- 60L
bottom <- 50L
gutter_x <- 60L
gutter_y <- 60L

panel_w <- floor((width  - left - right  - gutter_x) / 2)
panel_h <- floor((height - top  - bottom - gutter_y) / 2)

# Per-panel nudges can be tuned if needed without touching the overall layout.
panels <- list(
  A = make_panel(panel_paths$A, "A", panel_w, panel_h, img_x_extra = 4,  img_y_extra = 2),
  B = make_panel(panel_paths$B, "B", panel_w, panel_h, img_x_extra = 10, img_y_extra = 2),
  C = make_panel(panel_paths$C, "C", panel_w, panel_h, img_x_extra = 4,  img_y_extra = 2),
  D = make_panel(panel_paths$D, "D", panel_w, panel_h, img_x_extra = 4,  img_y_extra = 2)
)

canvas <- image_blank(width = width, height = height, color = "white")

positions <- list(
  A = c(left, top),
  B = c(left + panel_w + gutter_x, top),
  C = c(left, top + panel_h + gutter_y),
  D = c(left + panel_w + gutter_x, top + panel_h + gutter_y)
)

for (nm in names(panels)) {
  pos <- positions[[nm]]
  canvas <- image_composite(canvas, panels[[nm]], offset = sprintf("+%d+%d", pos[[1]], pos[[2]]))
}

png_out <- file.path(outdir, "Supplementary_Figure_S7.png")
pdf_out <- file.path(outdir, "Supplementary_Figure_S7.pdf")

image_write(canvas, path = png_out, format = "png", density = "300x300")
image_write(canvas, path = pdf_out, format = "pdf", density = "300x300")

message("[INFO] wrote: ", png_out)
message("[INFO] wrote: ", pdf_out)
