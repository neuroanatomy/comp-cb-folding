#--------------------------------
# Fitting of evolutionary models
# Brain phenotypes only
# 25 January 2023
#
# Warning: this script takes >1h
# to execute.
#--------------------------------

source("./common.R")

Yb <- as.matrix(df[c("Log10Area.cb", "Log10Length.cb", "Log10ThicknessMedian",
                    "Log10PeriodMedian", "Log10WidthMedian", "Log10Area.ctx",
                    "Log10Length.ctx")])
rownames(Yb) <- df$English_name
Yb[is.nan(Yb)] <- NA

# impute missing data

fitb_1 <- mvBM(tree_en, Yb, model = "BM1", control = list(maxit = 1000),
    param = list(constraint = "equaldiagonal"))
impb_1 <- estim(tree_en, Yb, fitb_1, asr = F)
ancb_1 <- estim(tree_en, impb_1$estimates, fitb_1, asr = T)
save(fitb_1, impb_1, ancb_1, file = "../data/derived/fitb_BM.RData")

fitb_2 <- mvOU(tree_en, Yb, model = "OU1",
    control = list(maxit = 1000),
    param = list(decompSigma = "equaldiagonal"))
impb_2 <- estim(tree_en, Yb, fitb_2, asr = F)
ancb_2 <- estim(tree_en, impb_2$estimates, fitb_2, asr = T)
save(fitb_2, impb_2, ancb_2, file = "../data/derived/fitb_OU.RData")

fitb_3 <- mvEB(tree_en, Yb, control = list(maxit = 1000))
impb_3 <- estim(tree_en, Yb, fitb_3, asr = F)
ancb_3 <- estim(tree_en, impb_3$estimates, fitb_3, asr = T)
save(fitb_3, impb_3, ancb_3, file = "../data/derived/fitb_EB.RData")

# fit models

Yb <- as.matrix(df[c("Log10Area.cb", "Log10Length.cb", "Log10ThicknessMedian",
                    "Log10PeriodMedian", "Log10WidthMedian", "Log10Area.ctx",
                    "Log10Length.ctx")])
rownames(Yb) <- df$English_name
Yb[is.nan(Yb)] <- NA

datab_BM <- list(Y = impb_1$estimates)
pb_BM <- mvgls(Y~1, data = datab_BM, tree_en,
    model = "BM", penalty = "RidgeArch")

datab_OU <- list(Y = impb_2$estimates)
pb_OU <- mvgls(Y~1, data = datab_OU, tree_en,
    model = "OU", penalty = "RidgeArch", upper = Inf)

datab_EB <- list(Y = impb_3$estimates)
pb_EB <- mvgls(Y~1, data = datab_EB, tree_en,
    model = "EB", penalty = "RidgeArch")

GIC(pb_BM)
GIC(pb_OU)
GIC(pb_EB)
