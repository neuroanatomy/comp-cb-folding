#------------------
# Phylogenetic PCA
# 25 Jan 2023
#------------------

source("./common.R")

# All phenotypes
#----------------

# load imputed data
load("../data/derived/fit_OU.RData")
dfOU <- imp_2$estimates
rownames(dfOU) <- df$Binomial_name_timetree

# fit mvgls
variables <- c(
  "Log10BodyWeight", "Log10BrainWeight",
  "Log10Area.cb", "Log10Area.ctx",
  "Log10Length.cb", "Log10Length.ctx",
  "Log10WidthMedian", "Log10PeriodMedian", "Log10ThicknessMedian")
Y <- as.matrix(dfOU[, variables])
p_OU <- mvgls(Y~1, data = list(Y = Y), tree,
    model = "OU1", penalty = "RidgeArch", upper = Inf)

# phylogenetic pca

pca_OU <- mvgls.pca(p_OU, plot = FALSE)
tmp <- matrix(rep(p_OU$coefficients,
              each = nrow(Y)),
              nrow = nrow(Y))
imp_pc <- (Y - tmp) %*% pca_OU$vectors

# save multivariate allometry (pc1)
pc1 <- data.frame(loadings = pca_OU$vectors[, 1], row.names = variables)
sink("../data/derived/pc1_loadings.txt")
print(pc1)
sink()

# <<Figure 8a>>
# pc1 vs pc2 plot
my_plot_pca(
  pca_OU$scores,
  pca_OU$values,
  df$English_name[match(rownames(pca_OU$scores), df$Binomial_name_timetree)],
  flip = c(1, -1)
  )
save_pdf("../data/derived/fig.8a.pdf", width = 14, height = 7)

# <<Figure 8b>>
# loadings
names <- c("Body weight", "Brain weight", "Cerebellar area",
           "Cerebral area", "Cerebellar length", "Cerebral length",
           "Folial width", "Folial perimeter", "Thickness")
p1 <- my_plot_loadings(names, pca_OU$vectors[, 1])
p2 <- my_plot_loadings(names, -pca_OU$vectors[, 2])
grid.arrange(p1, p2, nrow = 2)
save_pdf("../data/derived/fig.8b.pdf", width = 5, height = 5)

# Only brain phenotypes
#-----------------------

load("../data/derived/fitb_OU.RData")
dfbOU <- impb_2$estimates
rownames(dfbOU) <- df$Binomial_name_timetree

# fit mvgls
Yb <- as.matrix(dfOU[, c("Log10Area.cb", "Log10Area.ctx", "Log10Length.cb",
              "Log10Length.ctx", "Log10WidthMedian", "Log10PeriodMedian",
              "Log10ThicknessMedian")])
pb_OU <- mvgls(Y~1,
              data = list(Y = Yb),
              tree,
              model = "OU1",
              penalty = "RidgeArch", upper = Inf)

# phylogenetic pca
pcab_OU <- mvgls.pca(pb_OU, plot = FALSE)
tmp <- matrix(rep(pb_OU$coefficients,
              each = nrow(Yb)),
              nrow = nrow(Yb))
impb_pc <- (Yb - tmp) %*% pcab_OU$vectors

# <<Figure 8c>>
# pc1 vs pc2 plot
my_plot_pca(
  pcab_OU$scores,
  pcab_OU$values,
  df$English_name[match(rownames(pca_OU$scores), df$Binomial_name_timetree)],
  flip = c(-1, 1)
)
save_pdf("../data/derived/fig.8c.pdf", width = 14, height = 7)

# <<Figure 8d>>
# loadings
names <- c("Cerebellar area", "Cerebral area", "Cerebellar length",
           "Cerebral length", "Folial width", "Folial perimeter",
           "Thickness")
p1 <- my_plot_loadings(names, -pcab_OU$vectors[, 1])
p2 <- my_plot_loadings(names, pcab_OU$vectors[, 2])
grid.arrange(p1, p2, nrow = 2)
save_pdf("../data/derived/fig.8d.pdf", width = 5, height = 5)
