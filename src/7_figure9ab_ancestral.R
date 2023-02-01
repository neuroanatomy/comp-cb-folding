#----------------------------
# Ancestral state estimation
# 25 January 2023
#----------------------------

source("./common.R")

# load mvMORPH imputation
load("../data/derived/fit_OU.RData")
dfOU <- imp_2$estimates
rownames(dfOU) <- df$English_name
anc <- anc_2$estimates

# fit OU model
Y <- as.matrix(dfOU[, c("Log10Area.cb", "Log10Length.cb",
    "Log10ThicknessMedian", "Log10PeriodMedian", "Log10WidthMedian",
    "Log10Area.ctx", "Log10Length.ctx",
    "Log10BodyWeight", "Log10BrainWeight")])

p_OU <- mvgls(Y~1,
    data = list(Y = Y),
    tree_en,
    model = "OU1",
    penalty = "RidgeArch", upper = Inf)

# phylogenetic pca

pca_OU <- mvgls.pca(p_OU, plot = FALSE)
tmp <- matrix(rep(p_OU$coefficients,
              each = nrow(Y)),
              nrow = nrow(Y))
imp_pc <- (Y - tmp) %*% pca_OU$vectors
tmp <- matrix(rep(p_OU$coefficients, each = nrow(anc)),
              nrow = nrow(anc))
anc_pc <- (anc - tmp) %*% pca_OU$vectors

# <<Figure 9a>>
my_plot_ancestral(tree_en, -imp_pc[, 1],
            rownames(Y), -anc_pc[, 1], "Ancestral PC1")
save_pdf("../data/derived/fig.9a1.pdf", width = 9, height = 12)
my_plot_ancestral(tree_en, -imp_pc[, 2],
                  rownames(Y), -anc_pc[, 2], "Ancestral PC2")
save_pdf("../data/derived/fig.9a2.pdf", width = 9, height = 12)


# Only brain phenotypes
#-----------------------

# mvMORPH imputation
load("../data/derived/fitb_OU.RData")
dfbOU <- impb_2$estimates
rownames(dfbOU) <- df$English_name
ancb <- ancb_2$estimates

# fit mvgls
Yb <- as.matrix(dfbOU[, c("Log10Area.cb", "Log10Length.cb",
          "Log10ThicknessMedian", "Log10PeriodMedian", "Log10WidthMedian",
          "Log10Area.ctx", "Log10Length.ctx")])

pb_OU <- mvgls(Y~1, data = list(Y = Yb), tree_en,
               model = "OU1", penalty = "RidgeArch")

pcab_OU <- mvgls.pca(pb_OU, plot = FALSE)
tmp <- matrix(rep(pb_OU$coefficients,
              each = nrow(Yb)),
              nrow = nrow(Yb))
impb_pc <- (Yb - tmp) %*% pcab_OU$vectors
tmp <- matrix(rep(pb_OU$coefficients, each = nrow(ancb)),
              nrow = nrow(ancb))
ancb_pc <- (ancb - tmp) %*% pcab_OU$vectors

# <<Figure 9b>>
# plot ancestral states of PC1 and PC2
my_plot_ancestral(tree_en, -impb_pc[, 1],
        rownames(Yb), -ancb_pc[, 1], "Ancestral PC1 (Brain)")
save_pdf("../data/derived/fig.9b1.pdf", width = 9, height = 12)
my_plot_ancestral(tree_en, impb_pc[, 2],
          rownames(Yb), ancb_pc[, 2], "Ancestral PC2 (Brain)")
save_pdf("../data/derived/fig.9b2.pdf", width = 9, height = 12)
