#------------------------------------
# Compare folial with measures
# to folial width measures from
# Sultan and Braitenberg (although
# they do not denote the same thing)
# 25 Jan 2023
#------------------------------------

source("./common.R")

# select variables

fw <- data.frame(
  `Binomial_name_timetree` = character(),
  `English_name` = character(),
  `Log10FolialLength` = double(),
  `Log10FolialWidth` = double(),
  `Log10CbSurface` = double(),
  `Log10CbWeight` = double(),
  `Log10BodyWeight` = double(),
  `Log10Area.cb` = double(),
  `Log10Length.cb` = double(),
  `Log10WidthMedian.cb` = double(),
  `Log10PeriodMedian.cb` = double(),
  `Log10BrainG` = double(),
  `Log10BodyG` = double()
)
for (row in sultan$Binomial_name_timetree) {
  if (sum(df$Binomial_name_timetree == row)) {
    r <- c(
      row,
      as.matrix(df[df$Binomial_name_timetree == row, "English_name"]),
      log10(as.matrix(sultan[sultan$Binomial_name_timetree == row,
        "Sum_of_all_folial_lengths_mm"])),
      log10(as.matrix(sultan[sultan$Binomial_name_timetree == row,
        "Average_width_of_the_folia_mm"])),
      log10(as.matrix(sultan[sultan$Binomial_name_timetree == row,
        "Cerebellar_surface_mm2"])),
      log10(as.matrix(sultan[sultan$Binomial_name_timetree == row,
        "Cerebellar_weight_g"])),
      log10(as.matrix(sultan[sultan$Binomial_name_timetree == row,
        "Body_weight_g"])),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10Area.cb"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10Length.cb"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10WidthMedian"]),
      as.matrix(df[df$Binomial_name_timetree == row, "Log10PeriodMedian"]),
      as.matrix(df[df$Binomial_name_timetree == row,
        "Log10BrainWeight"]) / log(10), # convert natural log to log10
      as.matrix(df[df$Binomial_name_timetree == row,
        "Log10BodyWeight"]) / log(10) # convert natural log to log10
    )
    fw[nrow(fw) + 1, ] <- r
  }
}
fw$Log10FolialLength <- as.numeric(fw$Log10FolialLength)
fw$Log10FolialWidth <- as.numeric(fw$Log10FolialWidth)
fw$Log10CbSurface <- as.numeric(fw$Log10CbSurface)
fw$Log10CbWeight <- as.numeric(fw$Log10CbWeight)
fw$Log10BodyWeight <- as.numeric(fw$Log10BodyWeight)
fw$Log10Area.cb <- as.numeric(fw$Log10Area.cb)
fw$Log10Length.cb <- as.numeric(fw$Log10Length.cb)
fw$Log10WidthMedian.cb <- as.numeric(fw$Log10WidthMedian.cb)
fw$Log10PeriodMedian.cb <- as.numeric(fw$Log10PeriodMedian.cb)
fw$Log10BrainG <- as.numeric(fw$Log10BrainG)
fw$Log10BodyG <- as.numeric(fw$Log10BodyG)

# <<Figure 6c>>
my_plot_allometry(fw$Log10FolialWidth, fw$Log10WidthMedian.cb,
               fw$English_name, 0,
               xlab = expression("Sultan's cerebellar folial width (log"[10]*"mm)"),
               ylab = expression("Cerebellar folial width (log"[10]*")"),
               main = "Folial width versus Sultan's folial width")
save_pdf("../data/derived/fig.6c.pdf")
