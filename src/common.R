#---------------------------
# Load libraries and define
# plotting functions
# 25 Jan 2023
#---------------------------

library(ape)
library(corrplot)
library(diagram)
library(dplyr)
library(factoextra)
library(geiger)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(igraph)
library(memoise)
library(mvMORPH)
library(nlme)
library(phytools)
library(Rphylopars)
library(svglite)

# Load data
source("./load_data.R")

catn <- function(...) {
    cat(...)
    cat("\n")
}

save_pdf <- function(path, width=10, height=10) {
  # save the current plot to pdf
  pdf(path, width, height)
  dev.set(which = 2)
  dev.copy(which = 4)
  dev.off()
}

my_plot_loadings <- function(names, loadings) {
  names <- factor(names, levels = names)
  p <- ggplot(
    data.frame(
      x = names,
      y = loadings),
    aes(x = x, y = y)) +
    geom_segment(aes(x = x, xend = x, y = 0, yend = y), color = "grey") +
    geom_point(size = 3) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("") +
    ylab("Loading")
  return(p)
}

my_plot_pcor <- function(M, row_names, thresh_nodes, thresh_edges, main) {
  G <- round(cor2pcor(M), digits = 3)
  N <- graph_from_adjacency_matrix(G,
                                   weighted = T,
                                   mode = "undirected", diag = F)
  layout_g <- abs(G)
  layout_n <- graph_from_adjacency_matrix(layout_g,
                                          weighted = T,
                                          mode = "undirected", diag = F)

  my_palette <- colorRampPalette(c("red", "gray", "blue"))(n = 100)
  maxw <- max(E(N)$weight)
  E(layout_n)$color <- my_palette[1 + 99 * (E(N)$weight + maxw) / (2 * maxw)]
  E(layout_n)$edge.label <- sprintf("%.3g", E(N)$weight)
  E(layout_n)$edge.label[abs(E(N)$weight) < thresh_nodes] <- ""
  edge_width <- E(layout_n)$weight
  edge_width[edge_width < thresh_edges] <- 0
  plot(layout_n,
       rescale = T,
       edge.width = 10 * edge_width,
       edge.label = E(layout_n)$edge.label,
       edge.label.family = "Helvetica",
       edge.label.color = "black",
       edge.label.cex = 0.7,
       vertex.label = row_names,
       vertex.label.family = "Helvetica",
       vertex.label.font = 2,
       vertex.label.cex = 0.8,
       vertex.color = "gray",
       vertex.label.color = "black",
       vertex.frame.color = "gray",
       main = main)
}

my_plot_pca <- function(scores, values, labels, flip=c(1, 1)) {
  ggpubr::ggscatter(
    data.frame(
      x = scores[, 1] * flip[1],
      y = scores[, 2] * flip[2],
      taxon = taxa[cl[labels]],
      name = gsub("_", " ", labels)),
    x = "x", y = "y",
    xlab = sprintf("Dim 1 (%.1f%%)", 100 * values[1] / sum(values)),
    ylab = sprintf("Dim 2 (%.1f%%)", 100 * values[2] / sum(values)),
    color = "taxon", shape = "taxon", size = 5,
    label = "name", show.legend.text = F, repel = T,
    font.label = c(12, "plain", "black"),
    ggtheme = theme_minimal(),
    ellipse = T,
    ellipse.level = 0.5,
    ellipse.type = "norm",
    ellipse.alpha = 0.1,
    ellipse.border.remove = T) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "dashed") +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed") +
  scale_shape_manual(
    values = rep(15:18, len = 8)) +
  scale_color_grey(
    start = 0,
    end = 0.8)
}

my_plot_allometry <- function(x, y, labels, beta, xlab, ylab, main) {
  df <- data.frame(x = x, y = y)

  # linear regression
  reg <- lm(y ~ x, df)
  pval <- summary(reg)$coefficients[8]
  print(coef(reg))
  print(summary(reg))
  intercept <- coef(reg)[1] + min(x) * (coef(reg)[2] - beta)

  # orthogonal regression
  pc <- prcomp(~x + y, df)
  ro <- pc$rotation
  ce <- pc$center
  sd <- pc$sdev
  print(sd)

  dd <- data.frame(
    x = c(ce[1], ce[1] + ro[1, 1] * sd[1], NA, ce[1], ce[1] + ro[1, 2] * sd[2]),
    y = c(ce[2], ce[2] + ro[2, 1] * sd[1], NA, ce[2], ce[2] + ro[2, 2] * sd[2]))
  beta2 <- ro[2, 1] / ro[1, 1]
  intercept2 <- beta2 * (-ce[1]) + ce[2]

  xmax <- max(na.omit(x))
  xmin <- min(na.omit(x))
  ymax <- max(na.omit(y))
  ymin <- min(na.omit(y))
  xpos <- 0.7 * xmax + (1 - 0.7) * xmin
  ypos <- 0.1 * ymax + (1 - 0.1) * ymin

  strformat <- sprintf("LRslope = %.4g\nORslope = %.4g\n", coef(reg)[2], beta2)
  if (beta) {
    strformat <- paste(strformat, sprintf("Isometry = %.4g\n", beta), sep = "")
  }
  strformat <- paste(
    strformat,
    sprintf("Rsqr = %.4g\np-value = %.4g", summary(reg)$r.squared, pval),
    sep = "")

  dfcl <- taxa[cl[labels]]

  ggplot(data = df,
    aes(x, y, color = dfcl, label = labels)) +
    geom_smooth(method = "lm", color = "steelblue4") +
    scale_shape_manual(values = rep(15:18, len = length(unique(dfcl)))) +
    scale_color_grey(start = 0, end = 0.8) +
    geom_point(aes(shape = dfcl), size = 3, stroke = 1) +
    geom_text_repel(color = "black", max.overlaps = 15) +
    theme_minimal() +
    annotate("text", x = xpos, y = ypos, hjust = 0, label = strformat) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) + {
      if (beta) geom_abline(
      intercept = intercept,
      slope = beta,
      linetype = "dotted")
    } +
    geom_abline(intercept = intercept2, slope = beta2, linetype = "dashed")
}

my_plot_ancestral <- function(tree, tips, tipnames, states, mytitle) {
  par(xpd = NA)
  par(oma = c(4, 8, 4, 2))
  names(tips) <- tipnames
  obj <- contMap(tree, tips, method = "user", anc.states = states, plot = "F")
  obj <- setMap(obj, colors = c("blue", "white", "red"))
  plot.contMap(obj, ftype = "reg")
  axisPhylo()
  title(mytitle)
}
