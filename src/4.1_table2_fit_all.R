#--------------------------------
# Fitting of evolutionary models
# All phenotypes
# 25 January 2023
#
# Warning: this script takes >1h
# to execute.
#--------------------------------

source("./common.R")

Y <- as.matrix(df[c("Log10Area.cb", "Log10Length.cb", "Log10ThicknessMedian",
 "Log10PeriodMedian", "Log10WidthMedian", "Log10Area.ctx", "Log10Length.ctx",
 "Log10BodyWeight", "Log10BrainWeight")])
rownames(Y) <- df$English_name
Y[is.nan(Y)] <- NA

# impute missing data

fit_1 <- mvBM(tree_en, Y, model = "BM1", control = list(maxit = 1000),
    param = list(constraint = "equaldiagonal"))
imp_1 <- estim(tree_en, Y, fit_1, asr = F)
anc_1 <- estim(tree_en, imp_1$estimates, fit_1, asr = T)
save(fit_1, imp_1, anc_1, file = "../data/derived/fit_BM.RData")

fit_2 <- mvOU(tree_en, Y, model = "OU1",
    control = list(maxit = 1000), param = list(decompSigma = "equaldiagonal"))
imp_2 <- estim(tree_en, Y, fit_2, asr = F)
anc_2 <- estim(tree_en, imp_2$estimates, fit_2, asr = T)
save(fit_2, imp_2, anc_2, file = "../data/derived/fit_OU.RData")

fit_3 <- mvEB(tree_en, Y, control = list(maxit = 1000))
imp_3 <- estim(tree_en, Y, fit_3, asr = F)
anc_3 <- estim(tree_en, imp_3$estimates, fit_3, asr = T)
save(fit_3, imp_3, anc_3, file = "../data/derived/fit_EB.RData")

# fit models

Y <- as.matrix(df[c("Log10Area.cb", "Log10Length.cb", "Log10ThicknessMedian",
 "Log10PeriodMedian", "Log10WidthMedian", "Log10Area.ctx", "Log10Length.ctx",
 "Log10BodyWeight", "Log10BrainWeight")])
rownames(Y) <- df$English_name
Y[is.nan(Y)] <- NA

data_BM <- list(Y = imp_1$estimates)
p_BM <- mvgls(Y~1, data = data_BM, tree_en, model = "BM", penalty = "RidgeArch")

data_OU <- list(Y = imp_2$estimates)
p_OU <- mvgls(Y~1, data = data_OU, tree_en, model = "OU",
    penalty = "RidgeArch", upper = Inf)

data_EB <- list(Y = imp_3$estimates)
p_EB <- mvgls(Y~1, data = data_EB, tree_en, model = "EB", penalty = "RidgeArch")

# <<Table 2>>
GIC(p_BM)
GIC(p_OU)
GIC(p_EB)
