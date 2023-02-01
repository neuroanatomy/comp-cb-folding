#--------------------------------------
# Compare length and area measurements
# to surface and volume measurements
# 25 Jan 2023
#--------------------------------------

source("./common.R")

# select variables

fw  <-  data.frame(
  `Binomial_name_timetree` = character(),
  `English_name` = character(),
  `Log10Brain` = double(),
  `Log10Cerebellum` = double(),
  `Log10MedialCb` = double(),
  `Log10Area.ctx` = double(),
  `Log10Length.ctx` = double(),
  `Log10Area.cb` = double(),
  `Log10Length.cb` = double(),
  `Log10WidthMedian.cb` = double(),
  `Log10PeriodMedian.cb` = double(),
  `Log10ThicknessMedian.cb` = double()
)
for (row in smaers$Binomial_name_timetree) {
  if (sum(df$Binomial_name_timetree == row)) {
    r <- c(
      row,
      as.matrix(df[df$Binomial_name_timetree == row, "English_name"]),
      log10(as.matrix(smaers[smaers$Binomial_name_timetree == row, "Brain"])),
      log10(as.matrix(smaers[smaers$Binomial_name_timetree == row,
        "Cerebellum"])),
      log10(as.matrix(smaers[smaers$Binomial_name_timetree == row, "Medial"])),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10Area.ctx"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10Length.ctx"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10Area.cb"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10Length.cb"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10WidthMedian"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10PeriodMedian"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10ThicknessMedian"])
    )
    fw[nrow(fw) + 1, ] <- r
  }
}
fw$Log10Brain <- as.numeric(fw$Log10Brain) + 3 # convert cm3 to mm3
fw$Log10Cerebellum <- as.numeric(fw$Log10Cerebellum) + 3 # convert cm3 to mm3
fw$Log10MedialCb <- as.numeric(fw$Log10MedialCb) + 3 # convert cm3 to mm3
fw$Log10Area.ctx <- as.numeric(fw$Log10Area.ctx)
fw$Log10Length.ctx <- as.numeric(fw$Log10Length.ctx)
fw$Log10Area.cb <- as.numeric(fw$Log10Area.cb)
fw$Log10Length.cb <- as.numeric(fw$Log10Length.cb)
fw$Log10WidthMedian.cb <- as.numeric(fw$Log10WidthMedian.cb)
fw$Log10PeriodMedian.cb <- as.numeric(fw$Log10PeriodMedian.cb)
fw$Log10ThicknessMedian.cb <- as.numeric(fw$Log10ThicknessMedian.cb)

# <<Figure 6a>>
my_plot_allometry(fw$Log10Brain, fw$Log10Area.ctx,
        fw$English_name, 0.67,
        xlab = expression("Smaers' cerebral volume (log"[10] * "mm"^3 * ")"),
        ylab = expression("Cerebral section area (log"[10]*"mm"^2*")"),
        main = "Cerebral section area vs. volume")
save_pdf("../data/derived/fig.6a.pdf")

# <<Figure 6b>>
my_plot_allometry(fw$Log10Cerebellum, fw$Log10Area.cb,
        fw$English_name, 0.67,
        xlab = expression("Smaers' cerebellar volume (log"[10] * "mm"^3 * ")"),
        ylab = expression("Cerebellar section area (log"[10]*"mm"^2*")"),
        main = "Cerebellar section area vs. volume")
save_pdf("../data/derived/fig.6b.pdf")
