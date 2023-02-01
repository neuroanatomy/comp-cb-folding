#------------------
# Folial allometry
# 25 January 2023
#------------------

source("./common.R")

# bivariate scatterplots

# <<Figure 11a>>
# width versus length/area
my_plot_allometry(
  df$Log10Area.cb, df$Log10WidthMedian,
  df$English_name, 0,
  expression("Cerebellar section area (log"[10] * "mm"^2 * ")"),
  expression("Folial width (log"[10] * "mm)"),
  "Folial width allometry")
save_pdf("../data/derived/fig.11a.pdf")

# <<Figure 11b>>
# period versus length/area
my_plot_allometry(
  df$Log10Area.cb, df$Log10PeriodMedian,
  df$English_name, 0,
  expression("Cerebellar section area (log"[10] * "mm"^2 * ")"),
  expression("Folial period (log"[10] * "mm)"),
  "Folial period allometry")
save_pdf("../data/derived/fig.11b.pdf")

# <<Figure 11c>>
# thickness versus length/area
my_plot_allometry(
  df$Log10Area.cb, df$Log10ThicknessMedian,
  df$English_name, 0,
  expression("Cerebellar section area (log"[10] * "mm"^2 * ")"),
  expression("Molecular layer thickness (log"[10] * "mm)"),
  "Molecular layer thickness allometry")
save_pdf("../data/derived/fig.11c.pdf")
