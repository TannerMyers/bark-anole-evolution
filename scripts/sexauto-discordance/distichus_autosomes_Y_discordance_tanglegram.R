#!/usr/bin/env Rscript
################################################################################
# DISTICHUS DISCORDANCE TANGLEGRAM GENERATOR (FIGURE 2 STYLE)
#
# Generates publication-grade vector PDF tanglegrams styled after Figure 2
# from Myers et al. (distichus_tree.pdf):
#   1. Autosomes vs. Y-Chromosome Tanglegram:
#        Output: discordance_tanglegram_Y_auto.pdf
#
# Design features:
#   - Rectangular phylograms facing each other
#   - Monophyletic clade-colored branches with black backbone / ancestral stems
#   - Colors strictly matching Sampling_PCA.Rmd palette
#   - Small black dots indicating strongly supported internal nodes (UFBoot >= 95)
#   - Clean vertical clade segments without end-ticks, positioned strictly outside tips
#   - Smooth translucent highlight ribbons connecting major discordant clades:
#       * A. d. sejunctus (shifting from northern properus to southern properus)
#       * A. d. dominicensis - SW and A. d. favillarum (crossing ribbons showing inversion)
#   - Fine grey spline lines connecting individual discordant tips
#   - Prominent typography with protective white halos guaranteeing zero overlap
#   - Clean single-page vector PDF output
#
# Usage:
#   Rscript distichus_autosomes_Y_discordance_tanglegram.R
################################################################################

suppressMessages({
  library(ape)
  library(phytools)
})

find_file <- function(candidates) {
  for (f in candidates) {
    if (file.exists(f)) return(normalizePath(f))
  }
  return(NULL)
}

norm <- function(x) sub("[.].*$", "", x)

# ------------------------------------------------------------------------------
# 1. METADATA & PALETTE
# ------------------------------------------------------------------------------
tt_file  <- find_file(c("tips_taxon.csv", "../trees/analysis/tips_taxon.csv", "trees/analysis/tips_taxon.csv"))
ind_file <- find_file(c("discordance_individual.csv", "../trees/analysis/discordance_individual.csv", "trees/analysis/discordance_individual.csv"))

if (!is.null(tt_file) && !is.null(ind_file)) {
  tt  <- read.csv(tt_file, stringsAsFactors = FALSE)
  ind <- read.csv(ind_file, stringsAsFactors = FALSE)
} else {
  stop("Could not locate tips_taxon.csv and discordance_individual.csv")
}

tip_taxon <- tt$taxon
names(tip_taxon) <- tt$tip

# Color palette strictly matching Sampling_PCA.Rmd
taxa_colors_fig2 <- c(
  ravitergum       = "#BBDB8B", # Pale leaf green
  west_ignigularis = "#93BBBD", # Soft muted teal / slate
  east_ignigularis = "#28C3D1", # Bright cyan turquoise
  properus         = "#BDA425", # Olive gold
  sejunctus        = "#87DB41", # Bright lime green
  biminiensis      = "#DBA163", # Tan / copper
  distichoides     = "#662822", # Dark maroon brown
  distichus        = "#7B3F8A", # Purple
  dapsilis         = "#104728", # Forest green
  ocior            = "#D1C5AE", # Light beige / khaki
  dominicensis_E   = "#2E4766", # Navy slate blue
  dominicensis_NW  = "#7C15E6", # Vivid purple
  dominicensis_W   = "#B83B15", # Rust red
  dominicensis_SW  = "#CE0452", # Crimson red
  favillarum       = "#FFE684", # Soft yellow
  aurifer          = "#DB0CC4", # Vivid magenta
  vinosus          = "#DB8CA3", # Soft rose
  suppar           = "#DB6E1B"  # Pumpkin orange
)

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
get_edge_colors <- function(tr, tip_taxon, pal) {
  n_edge <- nrow(tr$edge)
  edge_col <- rep("#111111", n_edge)
  n_tips <- length(tr$tip.label)
  for (i in 1:n_edge) {
    child <- tr$edge[i, 2]
    if (child <= n_tips) {
      edge_col[i] <- pal[tip_taxon[tr$tip.label[child]]]
    } else {
      desc <- getDescendants(tr, child)
      desc_tips <- desc[desc <= n_tips]
      taxa <- unique(tip_taxon[tr$tip.label[desc_tips]])
      if (length(taxa) == 1) edge_col[i] <- pal[taxa]
    }
  }
  return(edge_col)
}

extract_high_support_internal <- function(tr) {
  if (is.null(tr$node.label)) return(integer(0))
  n_tips <- length(tr$tip.label)
  high_nodes <- integer(0)
  for (i in seq_along(tr$node.label)) {
    nl <- tr$node.label[i]
    if (nl == "" || is.na(nl)) next
    parts <- strsplit(nl, "/")[[1]]
    val <- suppressWarnings(as.numeric(parts[length(parts)]))
    if (!is.na(val) && val >= 95) high_nodes <- c(high_nodes, n_tips + i)
  }
  return(high_nodes)
}

sigmoid_curve <- function(x1, y1, x2, y2, n = 100) {
  t <- seq(0, 1, length.out = n)
  s <- 3 * t^2 - 2 * t^3
  list(x = x1 + t * (x2 - x1), y = y1 + s * (y2 - y1))
}

group_taxon <- function(tax) {
  if (tax %in% c("biminiensis", "distichoides", "distichus", "dapsilis", "ocior")) return("bahamas")
  return(tax)
}

flip_y <- function(y) 1 - y

# ------------------------------------------------------------------------------
# 3. AUTOSOMES VS. Y-CHROMOSOME
# ------------------------------------------------------------------------------
t1_f <- find_file(c("rooted_autosomes.nwk", "../trees/analysis/rooted_autosomes.nwk", "trees/analysis/rooted_autosomes.nwk"))
t2_f <- find_file(c("rooted_y.nwk", "../trees/analysis/rooted_y.nwk", "trees/analysis/rooted_y.nwk"))

if (is.null(t1_f) || is.null(t2_f)) {
  stop("Could not locate rooted_autosomes.nwk or rooted_y.nwk")
}

t_auto <- read.tree(t1_f)
t_y    <- read.tree(t2_f)

t_auto$tip.label <- norm(t_auto$tip.label)
t_y$tip.label    <- norm(t_y$tip.label)
t_auto <- keep.tip(t_auto, tt$tip)
t_y    <- keep.tip(t_y,    tt$tip)

cat("Optimizing cophylo rotation for Autosomes vs Y...\n")
xx <- cophylo(t_auto, t_y, assoc = cbind(t_auto$tip.label, t_y$tip.label), rotate = TRUE)

# Align the southern outgroup clade rotation in the Y tree so vinosus -> aurifer -> suppar matches autosomes perfectly
mrca_vin_aur <- findMRCA(xx$trees[[2]], tt$tip[tt$taxon %in% c("aurifer", "vinosus")])
mrca_south   <- findMRCA(xx$trees[[2]], tt$tip[tt$taxon %in% c("suppar", "aurifer", "vinosus")])
xx$trees[[2]] <- rotate(xx$trees[[2]], mrca_south)
xx$trees[[2]] <- rotate(xx$trees[[2]], mrca_vin_aur)

edge_col_left  <- get_edge_colors(xx$trees[[1]], tip_taxon, taxa_colors_fig2)
edge_col_right <- get_edge_colors(xx$trees[[2]], tip_taxon, taxa_colors_fig2)

high_supp_left  <- extract_high_support_internal(xx$trees[[1]])
high_supp_right <- extract_high_support_internal(xx$trees[[2]])

pdf(file = NULL)
plot(xx,
     part      = 0.36,
     link.type = "curved",
     link.lwd  = 0.0001,
     link.col  = "transparent",
     edge.col  = list(left = edge_col_left, right = edge_col_right),
     edge.lwd  = 1.8,
     ftype     = c("off", "off"),
     tip.lty   = 0,
     mar       = c(2, 1, 3.5, 1))
pp <- get("last_plot.cophylo", envir = .PlotPhyloEnv)
dev.off()

outfile_ay <- "discordance_tanglegram_Y_auto.pdf"
pdf(outfile_ay, width = 13, height = 15)
par(mar = c(2, 1, 3.5, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
plot.new()
plot.window(xlim = pp$left$x.lim, ylim = c(-0.02, 1.02))

left_edge  <- xx$trees[[1]]$edge
right_edge <- xx$trees[[2]]$edge

left_xx  <- pp$left$xx;  left_yy  <- flip_y(pp$left$yy)
right_xx <- pp$right$xx; right_yy <- flip_y(pp$right$yy)

# Draw left tree branches
for (i in 1:nrow(left_edge)) {
  p <- left_edge[i, 1]; c <- left_edge[i, 2]
  lines(c(left_xx[p], left_xx[p]), c(left_yy[p], left_yy[c]), col = edge_col_left[i], lwd = 1.8)
  lines(c(left_xx[p], left_xx[c]), c(left_yy[c], left_yy[c]), col = edge_col_left[i], lwd = 1.8)
}

# Draw right tree branches
for (i in 1:nrow(right_edge)) {
  p <- right_edge[i, 1]; c <- right_edge[i, 2]
  lines(c(right_xx[p], right_xx[p]), c(right_yy[p], right_yy[c]), col = edge_col_right[i], lwd = 1.8)
  lines(c(right_xx[p], right_xx[c]), c(right_yy[c], right_yy[c]), col = edge_col_right[i], lwd = 1.8)
}

# Node support dots: internal nodes with UFBoot >= 95
points(left_xx[high_supp_left], left_yy[high_supp_left], pch = 16, cex = 0.50, col = "black")
points(right_xx[high_supp_right], right_yy[high_supp_right], pch = 16, cex = 0.50, col = "black")

# Tip coordinates
n_tips <- length(xx$trees[[1]]$tip.label)
lx <- left_xx[1:n_tips]; ly <- left_yy[1:n_tips]
rx <- right_xx[1:n_tips]; ry <- right_yy[1:n_tips]
names(ly) <- xx$trees[[1]]$tip.label; names(ry) <- xx$trees[[2]]$tip.label
names(lx) <- xx$trees[[1]]$tip.label; names(rx) <- xx$trees[[2]]$tip.label

left_tax_grp  <- sapply(tip_taxon[xx$trees[[1]]$tip.label], group_taxon)
right_tax_grp <- sapply(tip_taxon[xx$trees[[2]]$tip.label], group_taxon)

compute_clade_blocks <- function(tax_grp, xs, ys, side = "left") {
  r <- rle(tax_grp)
  blocks <- list()
  curr <- 1
  for (b in seq_along(r$lengths)) {
    g_tax <- r$values[b]
    g_len <- r$lengths[b]
    g_tips <- curr:(curr + g_len - 1)
    curr <- curr + g_len
    
    y_min <- min(ys[g_tips])
    y_max <- max(ys[g_tips])
    x_val <- if (side == "left") max(xs[g_tips]) + 0.012 else min(xs[g_tips]) - 0.012
    
    if (y_max - y_min < 0.008) {
      y_mid <- (y_min + y_max) / 2
      y_min <- y_mid - 0.004
      y_max <- y_mid + 0.004
    }
    col <- if (g_tax == "bahamas") taxa_colors_fig2["biminiensis"] else taxa_colors_fig2[g_tax]
    blocks[[b]] <- list(tax = g_tax, y_min = y_min, y_max = y_max, y_mid = (y_min + y_max)/2, x_bar = x_val, col = col, len = g_len, tips = g_tips)
  }
  
  # Align southern clade brackets (vinosus, aurifer, suppar) along the same x_bar
  south_taxa <- c("vinosus", "aurifer", "suppar")
  south_idx <- which(sapply(blocks, function(b) b$tax %in% south_taxa))
  if (length(south_idx) > 0) {
    aligned_x <- if (side == "left") {
      max(sapply(blocks[south_idx], function(b) b$x_bar))
    } else {
      min(sapply(blocks[south_idx], function(b) b$x_bar))
    }
    for (idx in south_idx) {
      blocks[[idx]]$x_bar <- aligned_x
    }
  }
  
  return(blocks)
}

left_blocks  <- compute_clade_blocks(left_tax_grp, lx, ly, "left")
right_blocks <- compute_clade_blocks(right_tax_grp, rx, ry, "right")

find_block <- function(blks, tax, y_approx) {
  cand <- which(sapply(blks, function(b) b$tax == tax))
  if (length(cand) == 1) return(blks[[cand]])
  dists <- sapply(cand, function(idx) abs(blks[[idx]]$y_mid - y_approx))
  return(blks[[cand[which.min(dists)]]])
}

draw_ribbon_bar_to_bar <- function(x1, y1_min, y1_max, x2, y2_min, y2_max, col, lty = 1, alpha = 0.25, border_col = NULL) {
  c_bot <- sigmoid_curve(x1, y1_min, x2, y2_min)
  c_top <- sigmoid_curve(x1, y1_max, x2, y2_max)
  poly_x <- c(c_bot$x, rev(c_top$x))
  poly_y <- c(c_bot$y, rev(c_top$y))
  if (is.null(border_col)) border_col <- adjustcolor(col, alpha.f = 0.85)
  polygon(poly_x, poly_y,
          col = adjustcolor(col, alpha.f = alpha),
          border = border_col, lwd = 1.4, lty = lty)
}

# A. sejunctus ribbon
sej_left  <- find_block(left_blocks, "sejunctus", 0.677)
sej_right <- find_block(right_blocks, "sejunctus", 0.340)
draw_ribbon_bar_to_bar(
  x1 = sej_left$x_bar, y1_min = sej_left$y_min, y1_max = sej_left$y_max,
  x2 = sej_right$x_bar, y2_min = sej_right$y_min, y2_max = sej_right$y_max,
  col = taxa_colors_fig2["sejunctus"], lty = 1, alpha = 0.30
)

# B. dominicensis_SW and favillarum inversion ribbons
sw_left   <- find_block(left_blocks, "dominicensis_SW", 0.061)
sw_right  <- find_block(right_blocks, "dominicensis_SW", 0.122)
draw_ribbon_bar_to_bar(
  x1 = sw_left$x_bar, y1_min = sw_left$y_min, y1_max = sw_left$y_max,
  x2 = sw_right$x_bar, y2_min = sw_right$y_min, y2_max = sw_right$y_max,
  col = taxa_colors_fig2["dominicensis_SW"], lty = 1, alpha = 0.28
)

fav_left  <- find_block(left_blocks, "favillarum", 0.110)
fav_right <- find_block(right_blocks, "favillarum", 0.072)
draw_ribbon_bar_to_bar(
  x1 = fav_left$x_bar, y1_min = fav_left$y_min, y1_max = fav_left$y_max,
  x2 = fav_right$x_bar, y2_min = fav_right$y_min, y2_max = fav_right$y_max,
  col = taxa_colors_fig2["favillarum"], lty = 1, alpha = 0.28
)

# C. dominicensis_E discordant subclade ribbon (11 tips shifting adjacent to Bahamas)
domE_disc <- c("4473", "4475", "4476", "4479", "4482", "4511", "670", "677", "703", "704", "710")
domE_blk_left  <- find_block(left_blocks, "dominicensis_E", 0.318)
domE_blk_right <- find_block(right_blocks, "dominicensis_E", 0.472)
draw_ribbon_bar_to_bar(
  x1 = domE_blk_left$x_bar, y1_min = min(ly[domE_disc]), y1_max = max(ly[domE_disc]),
  x2 = domE_blk_right$x_bar, y2_min = min(ry[domE_disc]), y2_max = max(ry[domE_disc]),
  col = taxa_colors_fig2["dominicensis_E"], lty = 1, alpha = 0.25
)

# Individual connector splines (fine grey)
ribbon_encapsulated_tips <- c(tt$tip[tt$taxon %in% c("sejunctus", "dominicensis_SW", "favillarum")], domE_disc)
disc_tips <- setdiff(ind$tip[ind$class != "concordant"], ribbon_encapsulated_tips)
for (k in seq_along(disc_tips)) {
  tip_k <- disc_tips[k]
  i1 <- match(tip_k, xx$trees[[1]]$tip.label)
  i2 <- match(tip_k, xx$trees[[2]]$tip.label)
  if (is.na(i1) || is.na(i2)) next
  
  x1 <- lx[i1] + 0.012; y1 <- ly[i1]
  x2 <- rx[i2] - 0.012; y2 <- ry[i2]
  curve_pts <- sigmoid_curve(x1, y1, x2, y2)
  lines(curve_pts$x, curve_pts$y, col = adjustcolor("#555555", alpha.f = 0.35), lwd = 0.75, lty = 1)
}

# Vertical clade bars
for (b in left_blocks)  lines(c(b$x_bar, b$x_bar), c(b$y_min, b$y_max), col = b$col, lwd = 3.8)
for (b in right_blocks) lines(c(b$x_bar, b$x_bar), c(b$y_min, b$y_max), col = b$col, lwd = 3.8)

# Connect aurifer tip on right tree to its aligned bar
aur_tip_right <- "LM1180"
aur_bar_right <- right_blocks[[17]]$x_bar
lines(c(rx[aur_tip_right], aur_bar_right), c(ry[aur_tip_right], ry[aur_tip_right]),
      col = taxa_colors_fig2["aurifer"], lwd = 1.4, lty = 1)

draw_label_halo <- function(x, y, label, col, adj, cex = 1.02) {
  sw <- strwidth(label, cex = cex)
  sh <- strheight(label, cex = cex)
  pad_x <- 0.005; pad_y <- 0.0025
  if (adj[1] == 0) {
    x1 <- x - pad_x; x2 <- x + sw + pad_x
  } else {
    x1 <- x - sw - pad_x; x2 <- x + pad_x
  }
  y1 <- y - sh/2 - pad_y; y2 <- y + sh/2 + pad_y
  rect(x1, y1, x2, y2, col = adjustcolor("white", alpha.f = 0.94), border = NA)
  text(x, y, labels = label, cex = cex, col = col, adj = adj)
}

left_labels_def <- list(
  list(name = expression(italic("A. d. ravitergum")),        col = taxa_colors_fig2["ravitergum"],       y = 0.948, x_ref = 1,  cex = 1.02),
  list(name = expression(italic("A. d. ignigularis") ~ plain("- W")),   col = taxa_colors_fig2["west_ignigularis"],  y = 0.865, x_ref = 2,  cex = 1.02),
  list(name = expression(italic("A. d. ignigularis") ~ plain("- E")),   col = taxa_colors_fig2["east_ignigularis"],  y = 0.787, x_ref = 3,  cex = 1.02),
  list(name = expression(italic("A. d. properus")),          col = taxa_colors_fig2["properus"],         y = 0.724, x_ref = 4,  cex = 1.02),
  list(name = expression(italic("A. d. sejunctus")),         col = taxa_colors_fig2["sejunctus"],        y = 0.677, x_ref = 5,  cex = 1.02),
  list(name = expression(italic("A. d. properus")),          col = taxa_colors_fig2["properus"],         y = 0.616, x_ref = 6,  cex = 1.02),
  list(name = expression(italic("A. d. ignigularis") ~ plain("- E")),   col = taxa_colors_fig2["east_ignigularis"],  y = 0.542, x_ref = 7,  cex = 1.02),
  list(name = expression(plain("Bahamas")),                 col = taxa_colors_fig2["biminiensis"],      y = 0.464, x_ref = 8,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- E")),  col = taxa_colors_fig2["dominicensis_E"],   y = 0.318, x_ref = 9,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- NW")), col = taxa_colors_fig2["dominicensis_NW"],  y = 0.216, x_ref = 10, cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- W")),  col = taxa_colors_fig2["dominicensis_W"],   y = 0.174, x_ref = 11, cex = 1.02),
  list(name = expression(italic("A. d. favillarum")),        col = taxa_colors_fig2["favillarum"],       y = 0.111, x_ref = 12, cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- SW")), col = taxa_colors_fig2["dominicensis_SW"],  y = 0.062, x_ref = 13, cex = 0.96),
  list(name = expression(italic("A. d. vinosus")),           col = taxa_colors_fig2["vinosus"],          y = 0.035, x_ref = 14, cex = 0.92),
  list(name = expression(italic("A. d. aurifer")),           col = taxa_colors_fig2["aurifer"],          y = 0.0205, x_ref = 15, cex = 0.92),
  list(name = expression(italic("A. d. suppar")),            col = taxa_colors_fig2["suppar"],           y = 0.006, x_ref = 16, cex = 0.92)
)

for (lb in left_labels_def) {
  blk <- left_blocks[[lb$x_ref]]
  x_pos <- blk$x_bar + 0.008
  draw_label_halo(x_pos, lb$y, lb$name, col = lb$col, adj = c(0, 0.5), cex = lb$cex)
}

right_labels_def <- list(
  list(name = expression(italic("A. d. ravitergum")),        col = taxa_colors_fig2["ravitergum"],       y = 0.948, x_ref = 1,  cex = 1.02),
  list(name = expression(italic("A. d. ignigularis") ~ plain("- W")),   col = taxa_colors_fig2["west_ignigularis"],  y = 0.854, x_ref = 2,  cex = 1.02),
  list(name = expression(italic("A. d. ignigularis") ~ plain("- E")),   col = taxa_colors_fig2["east_ignigularis"],  y = 0.765, x_ref = 3,  cex = 1.02),
  list(name = expression(italic("A. d. properus")),          col = taxa_colors_fig2["properus"],         y = 0.693, x_ref = 4,  cex = 1.02),
  list(name = expression(plain("Bahamas")),                 col = taxa_colors_fig2["biminiensis"],      y = 0.619, x_ref = 5,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- E")),  col = taxa_colors_fig2["dominicensis_E"],   y = 0.472, x_ref = 6,  cex = 1.02),
  list(name = expression(italic("A. d. sejunctus")),         col = taxa_colors_fig2["sejunctus"],        y = 0.340, x_ref = 10, cex = 1.02),
  list(name = expression(italic("A. d. properus")),          col = taxa_colors_fig2["properus"],         y = 0.276, x_ref = 11, cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- NW")), col = taxa_colors_fig2["dominicensis_NW"],  y = 0.215, x_ref = 12, cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- W")),  col = taxa_colors_fig2["dominicensis_W"],   y = 0.174, x_ref = 13, cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- SW")), col = taxa_colors_fig2["dominicensis_SW"],  y = 0.122, x_ref = 14, cex = 0.96),
  list(name = expression(italic("A. d. favillarum")),        col = taxa_colors_fig2["favillarum"],       y = 0.072, x_ref = 15, cex = 1.02),
  list(name = expression(italic("A. d. vinosus")),           col = taxa_colors_fig2["vinosus"],          y = 0.035, x_ref = 18, cex = 0.92),
  list(name = expression(italic("A. d. aurifer")),           col = taxa_colors_fig2["aurifer"],          y = 0.0205, x_ref = 17, cex = 0.92),
  list(name = expression(italic("A. d. suppar")),            col = taxa_colors_fig2["suppar"],           y = 0.006, x_ref = 16, cex = 0.92)
)

for (rb in right_labels_def) {
  blk <- right_blocks[[rb$x_ref]]
  x_pos <- blk$x_bar - 0.008
  draw_label_halo(x_pos, rb$y, rb$name, col = rb$col, adj = c(1, 0.5), cex = rb$cex)
}

mtext("autosomal phylogeny", side = 3, line = 1.2, at = -0.34, cex = 1.45, font = 2)
mtext("Y-chromosome phylogeny", side = 3, line = 1.2, at = 0.34, cex = 1.45, font = 2)

dev.off()
cat("Saved:", outfile_ay, "\n")
