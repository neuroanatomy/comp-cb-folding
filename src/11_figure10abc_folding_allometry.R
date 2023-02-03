#-------------------
# Folding allometry
# 25 January 2023
#-------------------

source("./common.R")

# bivariate scatterplots

# <<Figure 10a>>
# Cb area versus length
my_plot_allometry(
  df$Log10Area.cb, df$Log10Length.cb,
  df$English_name, 0.5,
  expression("Cerebellum section area (log"[10] * "mm"^2 * ")"),
  expression("Cerebellum section length (log"[10] * "mm)"),
  "Cerebellum section length vs. area allometry")
save_pdf("../data/derived/fig.10a.pdf")

# <<Figure 10b>>
# Ctx area versus length
my_plot_allometry(
  df$Log10Area.ctx, df$Log10Length.ctx,
  df$English_name, 0.5,
  expression("Cerebrum section area (log"[10] * "mm"^2 * ")"),
  expression("Cerebrum section length (log"[10] * "mm)"),
  "Cerebrum section length vs. area allometry")
save_pdf("../data/derived/fig.10b.pdf")

# <<Figure 10c>>
# Cb Area ~ Ctx Area
my_plot_allometry(
  df$Log10Area.ctx, df$Log10Area.cb,
  df$English_name, 1,
  expression("Cerebrum section area (log"[10] * "mm"^2 * ")"),
  expression("Cerebellum section area (log"[10] * "mm"^2 * ")"),
  "Cerebellum vs. cerebrum section area allometry")
save_pdf("../data/derived/fig.10c.pdf")
