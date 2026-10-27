## Figure 1 (chessboard toy example) re-drawn from data recovered from the
## original vector PDFs (see fig1_extract_from_pdf.py).  No data are simulated.
## Usage: Rscript fig1_plot.R [<dir with fig1_*_from_pdf.csv.gz> [<outdir>]]
##   or, after installing ReMPCA:
##   source(system.file("reproducibility", "figure1", "fig1_plot.R", package = "ReMPCA"))
## Figure 1 is a conceptual illustration. This folder provides the plotted values
## (recovered from the vector graphics of the figure in the earlier version of the
## manuscript by fig1_extract_from_pdf.py) and the script that draws the figure.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- system.file("reproducibility", "figure1", package = "ReMPCA")
indir  <- if (length(args) >= 1) args[1] else if (nzchar(pkg_dir)) pkg_dir else "."
outdir <- if (length(args) >= 2) args[2] else Sys.getenv("REMPCA_OUTPUT_DIR", "ReMPCA_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
SHOW_KEY <- TRUE          # small key "True vector / Estimated vector" under the grid
FONT     <- if (capabilities("cairo")) "Arial" else "sans"   # Arial is embedded by cairo_pdf
FS_TITLE <- 8.5           # panel titles, row labels
FS_SMALL <- 7.5           # tick labels, key

H <- read.csv(file.path(indir, "fig1_heatmaps_from_pdf.csv.gz"), stringsAsFactors = FALSE)
L <- read.csv(file.path(indir, "fig1_vectors_from_pdf.csv.gz"),  stringsAsFactors = FALSE)

img_lev   <- c("True", "Noisy", "SVD", "Smooth", "Group lasso", "Smooth + group lasso")
img_lab   <- c("True", "Noisy", "SVD", "Smooth", "Group lasso", "Smooth +\ngroup lasso")
meth_lev  <- c("SVD", "Smooth", "Group lasso", "Smooth + group lasso")
meth_col  <- c("SVD" = "#E69F00", "Smooth" = "#009E73",
               "Group lasso" = "#CC79A7", "Smooth + group lasso" = "#0072B2")
truth_col <- "grey70"
lwd2ggp   <- function(lwd) lwd / .pt      # R lwd units -> ggplot2 linewidth
LW_TRUTH  <- lwd2ggp(1.6); LW_EST <- lwd2ggp(0.8)

base_theme <- theme_void(base_family = FONT, base_size = FS_TITLE) +
  theme(strip.text = element_text(size = FS_TITLE, colour = "black",
                                  margin = margin(0, 0, 2, 0)),
        plot.margin = margin(1, 2, 1, 2))

## ---- top row: six heatmaps (row 1 at top, column 1 at left; gray = recovered fill) ----
H$panel <- factor(img_lab[match(H$image, img_lev)], levels = img_lab)
H$fill  <- grDevices::rgb(H$gray, H$gray, H$gray)
nr <- max(H$row); nc <- max(H$col)
p_img <- ggplot(H, aes(col, row, fill = fill)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  facet_wrap(~panel, nrow = 1) +
  base_theme +
  theme(aspect.ratio = 1,
        strip.text = element_text(size = FS_TITLE, vjust = 0, lineheight = 0.9,
                                  margin = margin(0, 0, 2.5, 0)),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.spacing.x = unit(0.075, "in"),
        plot.margin = margin(-4, 2, 7, 2))

## ---- bottom: 2 x 4 grid of singular-vector plots ----
## value_std: per-strip affine rescaling of the original page y such that the true
## vector's two levels sit at 0 and 1 (the original strips had no y axis).
L$method <- factor(L$method, levels = meth_lev)
L$role   <- factor(L$role, levels = c("truth", "estimate"))
row_lab  <- c(u = 'bolditalic(u)~"(rows)"', v = 'bolditalic(v)~"(columns)"')

vec_plot <- function(vec, show_titles, breaks) {
  d  <- L[L$vector == vec, ]
  d$rowlab <- row_lab[vec]
  n  <- max(d$index)
  ggplot(d, aes(index, value_std)) +
    geom_line(data = d[d$role == "truth", ], colour = truth_col, linewidth = LW_TRUTH,
              lineend = "round", linejoin = "round") +
    geom_line(data = d[d$role == "estimate", ], aes(colour = method), linewidth = LW_EST,
              lineend = "round", linejoin = "round") +
    scale_colour_manual(values = meth_col, guide = "none") +
    scale_x_continuous(breaks = breaks, limits = c(1, n), expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = 0.06)) +
    facet_grid(rowlab ~ method, switch = "y",
               labeller = labeller(rowlab = label_parsed, method = label_value)) +
    base_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
          panel.spacing.x = unit(0.13, "in"),
          strip.placement = "outside",
          strip.text.x = if (show_titles) element_text(size = FS_TITLE, margin = margin(0, 0, 2.5, 0))
                         else element_blank(),
          strip.text.y.left = element_text(size = FS_TITLE, angle = 90, margin = margin(0, 3, 0, 0)),
          axis.text.x = element_text(size = FS_SMALL, colour = "grey20", margin = margin(1.2, 0, 0, 0)),
          axis.ticks.x = element_line(colour = "black", linewidth = 0.25),
          axis.ticks.length.x = unit(1.6, "pt"))
}
p_u <- vec_plot("u", TRUE,  c(1, 20, 40, 60, 80))
p_v <- vec_plot("v", FALSE, c(1, 20, 40, 60))

## ---- key ----
key_df <- data.frame(x = 0, y = 0)
p_key <- ggplot() +
  annotate("segment", x = 0.00, xend = 0.045, y = 0, yend = 0, colour = truth_col,
           linewidth = LW_TRUTH, lineend = "round") +
  annotate("text", x = 0.055, y = 0, label = "True vector", hjust = 0, size = FS_SMALL / .pt,
           family = FONT) +
  annotate("segment", x = 0.24 + 0.03 * (0:3), xend = 0.24 + 0.03 * (0:3) + 0.022, y = 0, yend = 0,
           colour = unname(meth_col), linewidth = LW_EST * 1.6, lineend = "round") +
  annotate("text", x = 0.365, y = 0, label = "Estimated vector (colour = method)", hjust = 0,
           size = FS_SMALL / .pt, family = FONT) +
  coord_cartesian(xlim = c(-0.02, 1), ylim = c(-1, 1), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(3, 2, 0, 2))

plots <- list(p_img, p_u, p_v)
## heights are panel-area heights: heatmap row fixed (square panels filling the width), line rows equal
hts <- c(65.6, 1, 1); hun <- c("pt", "null", "null")
if (SHOW_KEY) { plots <- c(plots, list(p_key)); hts <- c(hts, 0.10); hun <- c(hun, "null") }
fig <- wrap_plots(plots, ncol = 1, heights = unit(hts, hun))

W <- 6.3; Ht <- 3.75
out_pdf <- file.path(outdir, "fig1_chessboard.pdf")
if (capabilities("cairo")) {
  cairo_pdf(out_pdf, width = W, height = Ht, family = FONT)
} else {
  pdf(out_pdf, width = W, height = Ht)   # falls back to the default sans font
}
print(fig)
invisible(dev.off())
cat("wrote", out_pdf, "\n")
