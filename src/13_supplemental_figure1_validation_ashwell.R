#--------------------------------------
# Compare length and area measurements
# to surface and volume measurements from Ashwell
# 31 Jan 2023
#--------------------------------------

source("./common.R")

# select variables

fw  <-  data.frame(
  `Binomial_name_timetree` = character(),
  `English_name` = character(),
  `Log10Brain` = double(),
  `Log10Cerebellum` = double(),
  `Log10CbExtSurface`=double(),
  `Log10CbPialSurface`=double(),
  `Log10Area.ctx`=double(),
  `Log10Length.ctx`=double(),
  `Log10Area.cb`=double(),
  `Log10Length.cb`=double(),
  `Log10WidthMedian.cb`=double(),
  `Log10PeriodMedian.cb`=double(),
  `Log10ThicknessMedian.cb`=double()
)
for (row in ashwell$Binomial_name_timetree) {
  if (sum(df$Binomial_name_timetree == row)) {
    r <- c(
      row,
      as.matrix(df[df$Binomial_name_timetree == row, "English_name"]),
      log10(as.matrix(ashwell[ashwell$Binomial_name_timetree == row, "Brain.volume"])),
      log10(as.matrix(ashwell[ashwell$Binomial_name_timetree == row, "Total.Cb.volume"])),
      log10(as.matrix(ashwell[ashwell$Binomial_name_timetree == row, "Cb.ext.surface..ESA."])),
      log10(as.matrix(ashwell[ashwell$Binomial_name_timetree == row, "Cb.pial.surface..PSA."])),
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
fw$Log10Brain <- as.numeric(fw$Log10Brain)
fw$Log10Cerebellum <- as.numeric(fw$Log10Cerebellum)
fw$Log10CbExtSurface <- as.numeric(fw$Log10CbExtSurface)
fw$Log10CbPialSurface <- as.numeric(fw$Log10CbPialSurface)
fw$Log10Area.ctx <- as.numeric(fw$Log10Area.ctx)
fw$Log10Length.ctx <- as.numeric(fw$Log10Length.ctx)
fw$Log10Area.cb <- as.numeric(fw$Log10Area.cb)
fw$Log10Length.cb <- as.numeric(fw$Log10Length.cb)
fw$Log10WidthMedian.cb <- as.numeric(fw$Log10WidthMedian.cb)
fw$Log10PeriodMedian.cb <- as.numeric(fw$Log10PeriodMedian.cb)
fw$Log10ThicknessMedian.cb <- as.numeric(fw$Log10ThicknessMedian.cb)

# <<Supplemental Figure 1a validation Ashwell>>
my_plot_allometry(fw$Log10Cerebellum, fw$Log10Area.cb,
        fw$English_name, 0.67,
        xlab = expression("Ashwell cerebellar volume (log"[10] * "mm" ^3 * ")"),
        ylab = expression("Cerebellum section area (log"[10] * "mm" ^2 * ")"),
        main = "Cerebellum section area vs. volume")
save_pdf("../data/derived/supplemental_figure1a_validation_ashwell.pdf")

# <<Supplemental Figure 1b validation Ashwell>>
my_plot_allometry(fw$Log10CbPialSurface, fw$Log10Length.cb,
        fw$English_name, 0.5,
        xlab = expression("Ashwell cerebellar pial surface (log"[10] * "mm" ^2 * ")"),
        ylab = expression("Cerebellum section length (log"[10] * "mm" * ")"),
        main = "Cerebellum section length vs. cerebellar pial surface")
save_pdf("../data/derived/supplemental_figure1b_validation_ashwell.pdf")



