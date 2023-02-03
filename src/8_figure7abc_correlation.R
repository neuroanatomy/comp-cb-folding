#-----------------------
# Correlation structure
# 25 Jan 2023
#-----------------------

source("./common.R")

# load imputed data
load("../data/derived/fit_OU.RData")
dfOU <- imp_2$estimates
rownames(dfOU) <- df$Binomial_name_timetree

# fit OU
Y <- as.matrix(dfOU[, c("Log10Area.cb", "Log10Length.cb",
    "Log10ThicknessMedian", "Log10PeriodMedian", "Log10WidthMedian",
    "Log10Area.ctx", "Log10Length.ctx", "Log10BodyWeight",
    "Log10BrainWeight")])
fit_mvgls_OU <- mvgls(Y~1, tree = tree, model = "OU1",
                      penalty = "RidgeArch", upper = Inf)

# <<Figure 7c>>
# plot partial correlation diagram
M <- round(cor2pcor(fit_mvgls_OU$sigma$Pinv), digits = 3)
rownames(M) <- c("Cerebellar\nsection\narea", "Cerebellar\nsection\nlength",
    "Thickness", "Folial\nperiod", "Folial\nwidth",
    "Cerebral\nsection\narea", "Cerebral\nsection\nlength",
    "Body\nweight", "Brain\nweight")
diag(M) <- 0
arr.lwd <- (M * 6)^2
M[M < 0.01] <- 0
set.seed(3)
my_plot_pcor(fit_mvgls_OU$sigma$Pinv, rownames(M), 0.15, 0.1, "Partial correlations diagram")
save_pdf("../data/derived/fig.7c.pdf", width = 6, height = 6)

# <<Figure 7b>>
# plot partial correlation matrix
MM <- cor2pcor(fit_mvgls_OU$sigma$Pinv)
rownames(MM) <- rownames(M)
corrplot(MM, method = "circle", tl.col = "black", tl.cex = 0.8)
save_pdf("../data/derived/fig.7b.pdf", width = 6, height = 6)

# <<Figure 7a>>
# plot correlation matrix
MM <- cov2cor(fit_mvgls_OU$sigma$Pinv)
rownames(MM) <- rownames(M)
corrplot(MM, method = "circle", tl.col = "black", tl.cex = 0.8)
save_pdf("../data/derived/fig.7a.pdf", width = 6, height = 6)
