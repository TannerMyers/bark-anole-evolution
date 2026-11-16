#!/usr/bin/env Rscript
################################################################################
# X1 VS. X2 SEX CHROMOSOME DISCORDANCE TANGLEGRAM (FIGURE 2 STYLE)
#
# Generates a publication-grade vector PDF tanglegram comparing the ancestral
# X (X1, 50% ddRAD) and neo-X (X2, 50% ddRAD) phylogenies (182 male tips),
# styled to match Figure 2 from Myers et al. (distichus_tree.pdf):
#   - Rectangular phylograms facing each other
#   - Colors strictly matching Sampling_PCA.Rmd palette
#   - Monophyletic clade-colored branches with black backbone / ancestral stems
#   - Small black dots indicating strongly supported internal nodes (UFBoot >= 95)
#   - Clean vertical clade segments without end-ticks, positioned strictly outside tips
#   - Clean open center (no ribbons for concordant clades)
#   - Fine grey spline lines showing only reticulating / discordant tips
#   - Prominent, publication-matched typography with protective white halos
#     (only scientific names italicized; cardinal directions and counts plain)
#   - Clean single-page vector PDF output
#
# Usage:
#   Rscript distichus_x1_x2_tanglegram.R
#   Rscript distichus_x1_x2_tanglegram.R [x1_tree] [x2_tree] [output_pdf]
################################################################################

suppressMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)

find_file <- function(candidates) {
  for (f in candidates) {
    if (file.exists(f)) return(normalizePath(f))
  }
  return(NULL)
}

# 1. Resolve Input Trees
if (length(args) >= 2) {
  x1_path <- args[1]
  x2_path <- args[2]
} else {
  x1_path <- find_file(c(
    "../trees/x1_distichus_50complete.treefile",
    "x1_distichus_50complete.treefile",
    "trees/x1_distichus_50complete.treefile",
    "../../trees/x1_distichus_50complete.treefile",
    "../x1_distichus_50complete.treefile"
  ))
  x2_path <- find_file(c(
    "../trees/x2_distichus_50complete.treefile",
    "x2_distichus_50complete.treefile",
    "trees/x2_distichus_50complete.treefile",
    "../../trees/x2_distichus_50complete.treefile",
    "../x2_distichus_50complete.treefile"
  ))
}

outfile <- if (length(args) >= 3) args[3] else "discordance_tanglegram_x1_x2.pdf"

if (is.null(x1_path) || is.null(x2_path)) {
  stop("Could not locate input tree files. Please provide them as arguments: Rscript distichus_x1_x2_tanglegram.R <x1_tree> <x2_tree> [out_pdf]")
}

cat("X1 tree:    ", x1_path, "\n")
cat("X2 tree:    ", x2_path, "\n")
cat("Output PDF: ", outfile, "\n")

norm <- function(x) sub("[.].*$", "", x)

# 2. Specimen Metadata
tt_path <- find_file(c(
  "tips_taxon.csv",
  "../trees/analysis/tips_taxon.csv",
  "trees/analysis/tips_taxon.csv",
  "../../trees/analysis/tips_taxon.csv"
))

if (!is.null(tt_path)) {
  tt <- read.csv(tt_path, stringsAsFactors = FALSE)
} else {
  stop("Could not locate tips_taxon.csv. Please ensure it is present in trees/analysis/.")
}

tip_taxon <- tt$taxon
names(tip_taxon) <- tt$tip

# Harmonized color palette strictly matching Sampling_PCA.Rmd
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

# 3. Load & Root Trees
tx1 <- read.tree(x1_path)
tx2 <- read.tree(x2_path)

tx1$tip.label <- norm(tx1$tip.label)
tx2$tip.label <- norm(tx2$tip.label)

auto_basal_tips <- c(
  "5594", "5593", "5586", "5683", "7459", "7458", "4329", "4250", "4242", "4257",
  "4041", "4276", "4277", "5100", "5099", "5101", "590", "5102", "5094", "606",
  "9840", "9839", "9837", "9833", "9832", "4549", "4548", "4428", "4384", "4382",
  "4030", "7697", "3705", "3703"
)

if (!is.rooted(tx1) || length(tx1$edge[tx1$edge[,1] == (length(tx1$tip.label)+1), 2]) > 2) {
  cat("Rooting X1 tree at 34-tip basal split...\n")
  tx1 <- root(tx1, outgroup = tx1$tip.label[tx1$tip.label %in% auto_basal_tips], resolve.root = TRUE)
}
if (!is.rooted(tx2) || length(tx2$edge[tx2$edge[,1] == (length(tx2$tip.label)+1), 2]) > 2) {
  cat("Rooting X2 tree at 34-tip basal split...\n")
  tx2 <- root(tx2, outgroup = tx2$tip.label[tx2$tip.label %in% auto_basal_tips], resolve.root = TRUE)
}

tx1 <- keep.tip(tx1, tt$tip)
tx2 <- keep.tip(tx2, tt$tip)

# Clade and backbone edge colors
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
      if (length(taxa) == 1) {
        edge_col[i] <- pal[taxa]
      }
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
    if (!is.na(val) && val >= 95) {
      high_nodes <- c(high_nodes, n_tips + i)
    }
  }
  return(high_nodes)
}

cat("Optimizing cophylo rotation for X1 vs X2...\n")
xx <- cophylo(tx1, tx2, assoc = cbind(tx1$tip.label, tx1$tip.label), rotate = TRUE)

edge_col_left  <- get_edge_colors(xx$trees[[1]], tip_taxon, taxa_colors_fig2)
edge_col_right <- get_edge_colors(xx$trees[[2]], tip_taxon, taxa_colors_fig2)

high_supp_left  <- extract_high_support_internal(xx$trees[[1]])
high_supp_right <- extract_high_support_internal(xx$trees[[2]])

# Compute cophylo coordinates off-screen
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

# In X1 vs X2, Tip 1 is aurifer (yy=0) and Tip 182 is west_ignigularis (yy=1)
flip_y <- function(y) y

# 4. Render Single-Page Vector PDF
pdf(outfile, width = 13, height = 15)
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

# Support dots on strongly supported internal nodes ONLY
points(left_xx[high_supp_left], left_yy[high_supp_left], pch = 16, cex = 0.50, col = "black")
points(right_xx[high_supp_right], right_yy[high_supp_right], pch = 16, cex = 0.50, col = "black")

n_tips <- length(xx$trees[[1]]$tip.label)
lx <- left_xx[1:n_tips]; ly <- left_yy[1:n_tips]
rx <- right_xx[1:n_tips]; ry <- right_yy[1:n_tips]
names(ly) <- xx$trees[[1]]$tip.label; names(ry) <- xx$trees[[2]]$tip.label
names(lx) <- xx$trees[[1]]$tip.label; names(rx) <- xx$trees[[2]]$tip.label

sigmoid_curve <- function(x1, y1, x2, y2, n = 100) {
  t <- seq(0, 1, length.out = n)
  s <- 3 * t^2 - 2 * t^3
  list(x = x1 + t * (x2 - x1), y = y1 + s * (y2 - y1))
}

group_taxon <- function(tax) {
  if (tax %in% c("biminiensis", "distichoides", "distichus", "dapsilis", "ocior")) return("bahamas")
  return(tax)
}

left_tax_grp  <- sapply(tip_taxon[xx$trees[[1]]$tip.label], group_taxon)
right_tax_grp <- sapply(tip_taxon[xx$trees[[2]]$tip.label], group_taxon)

# Compute clade blocks for vertical bars
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
  
  # Align southern clade brackets (vinosus, suppar, aurifer) along the same x_bar
  south_taxa <- c("vinosus", "suppar", "aurifer")
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

# Reticulating / discordant connector splines for tips that shift in the Eastern complex:
east_tips <- tt$tip[tt$taxon == "east_ignigularis"]
prop_tips <- tt$tip[tt$taxon == "properus"]

outside_east <- east_tips[ry[east_tips] < 0.768]
outside_prop <- prop_tips[ry[prop_tips] > 0.669]

reticulating_tips <- c(outside_east, outside_prop)

for (tip_id in reticulating_tips) {
  i1 <- match(tip_id, xx$trees[[1]]$tip.label)
  i2 <- match(tip_id, xx$trees[[2]]$tip.label)
  x1 <- lx[i1] + 0.012; y1 <- ly[i1]
  x2 <- rx[i2] - 0.012; y2 <- ry[i2]
  crv <- sigmoid_curve(x1, y1, x2, y2)
  lines(crv$x, crv$y, col = adjustcolor("#555555", alpha.f = 0.50), lwd = 0.95, lty = 1)
}

# Draw vertical clade bars (no ticks, strictly outside tips)
for (b in left_blocks)  lines(c(b$x_bar, b$x_bar), c(b$y_min, b$y_max), col = b$col, lwd = 3.8)
for (b in right_blocks) lines(c(b$x_bar, b$x_bar), c(b$y_min, b$y_max), col = b$col, lwd = 3.8)

# Halo label drawer
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

# Left labels with calibrated spacing (guaranteed zero overlap)
left_label_defs <- list(
  list(name = expression(italic("A. d. ignigularis") ~ plain("- W")), col = taxa_colors_fig2["west_ignigularis"], y = 0.964, x_ref = 19, cex = 1.02),
  list(name = expression(italic("A. d. ravitergum")),                  col = taxa_colors_fig2["ravitergum"],       y = 0.870, x_ref = 18, cex = 1.02),
  list(name = expression(italic("A. d. ignigularis") ~ plain("- E")), col = taxa_colors_fig2["east_ignigularis"], y = 0.770, x_ref = 17, cex = 1.02),
  list(name = expression(italic("A. d. properus")),                    col = taxa_colors_fig2["properus"],         y = 0.690, x_ref = 16, cex = 1.02),
  list(name = expression(italic("A. d. properus")),                    col = taxa_colors_fig2["properus"],         y = 0.620, x_ref = 14, cex = 1.02),
  list(name = expression(italic("A. d. sejunctus")),                   col = taxa_colors_fig2["sejunctus"],        y = 0.540, x_ref = 13, cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- E")), col = taxa_colors_fig2["dominicensis_E"],   y = 0.423, x_ref = 9,  cex = 1.02),
  list(name = expression(plain("Bahamas")),                           col = taxa_colors_fig2["biminiensis"],      y = 0.276, x_ref = 8,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- W")), col = taxa_colors_fig2["dominicensis_W"],   y = 0.191, x_ref = 7,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- NW")),col = taxa_colors_fig2["dominicensis_NW"],  y = 0.149, x_ref = 6,  cex = 1.02),
  list(name = expression(italic("A. d. favillarum")),                  col = taxa_colors_fig2["favillarum"],       y = 0.110, x_ref = 5,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- SW")),col = taxa_colors_fig2["dominicensis_SW"],  y = 0.062, x_ref = 4,  cex = 0.98),
  list(name = expression(italic("A. d. vinosus")),                     col = taxa_colors_fig2["vinosus"],          y = 0.042, x_ref = 3,  cex = 0.94),
  list(name = expression(italic("A. d. suppar")),                      col = taxa_colors_fig2["suppar"],           y = 0.022, x_ref = 2,  cex = 0.94),
  list(name = expression(italic("A. d. aurifer")),                     col = taxa_colors_fig2["aurifer"],          y = 0.003, x_ref = 1,  cex = 0.94)
)

for (lb in left_label_defs) {
  blk <- left_blocks[[lb$x_ref]]
  draw_label_halo(blk$x_bar + 0.008, lb$y, lb$name, col = lb$col, adj = c(0, 0.5), cex = lb$cex)
}

# Right labels with calibrated spacing (guaranteed zero overlap)
right_label_defs <- list(
  list(name = expression(italic("A. d. ignigularis") ~ plain("- W")), col = taxa_colors_fig2["west_ignigularis"], y = 0.964, x_ref = 21, cex = 1.02),
  list(name = expression(italic("A. d. ravitergum")),                  col = taxa_colors_fig2["ravitergum"],       y = 0.870, x_ref = 20, cex = 1.02),
  list(name = expression(italic("A. d. ignigularis") ~ plain("- E")), col = taxa_colors_fig2["east_ignigularis"], y = 0.770, x_ref = 19, cex = 1.02),
  list(name = expression(italic("A. d. properus")),                    col = taxa_colors_fig2["properus"],         y = 0.690, x_ref = 16, cex = 1.02),
  list(name = expression(italic("A. d. properus")),                    col = taxa_colors_fig2["properus"],         y = 0.620, x_ref = 14, cex = 1.02),
  list(name = expression(italic("A. d. sejunctus")),                   col = taxa_colors_fig2["sejunctus"],        y = 0.535, x_ref = 13, cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- E")), col = taxa_colors_fig2["dominicensis_E"],   y = 0.423, x_ref = 10, cex = 1.02),
  list(name = expression(plain("Bahamas")),                           col = taxa_colors_fig2["biminiensis"],      y = 0.276, x_ref = 9,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- W")), col = taxa_colors_fig2["dominicensis_W"],   y = 0.191, x_ref = 8,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- NW")),col = taxa_colors_fig2["dominicensis_NW"],  y = 0.149, x_ref = 7,  cex = 1.02),
  list(name = expression(italic("A. d. favillarum")),                  col = taxa_colors_fig2["favillarum"],       y = 0.110, x_ref = 6,  cex = 1.02),
  list(name = expression(italic("A. d. dominicensis") ~ plain("- SW")),col = taxa_colors_fig2["dominicensis_SW"],  y = 0.062, x_ref = 5,  cex = 0.98),
  list(name = expression(italic("A. d. vinosus")),                     col = taxa_colors_fig2["vinosus"],          y = 0.042, x_ref = 4,  cex = 0.94),
  list(name = expression(italic("A. d. suppar")),                      col = taxa_colors_fig2["suppar"],           y = 0.022, x_ref = 3,  cex = 0.94),
  list(name = expression(italic("A. d. aurifer")),                     col = taxa_colors_fig2["aurifer"],          y = 0.003, x_ref = 1,  cex = 0.94)
)

for (rb in right_label_defs) {
  blk <- right_blocks[[rb$x_ref]]
  draw_label_halo(blk$x_bar - 0.008, rb$y, rb$name, col = rb$col, adj = c(1, 0.5), cex = rb$cex)
}

# Publication Headers
mtext(expression(plain("ancestral X (") * italic("X1") * plain(") phylogeny")), side = 3, line = 1.2, at = -0.34, cex = 1.45, font = 2)
mtext(expression(plain("neo-X (") * italic("X2") * plain(") phylogeny")), side = 3, line = 1.2, at = 0.34, cex = 1.45, font = 2)

dev.off()
cat("Successfully generated single-page vector PDF:", outfile, "\n")
