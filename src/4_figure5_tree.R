#------------------------
# Plot phylogenetic tree
# 25 Jan 2023
#------------------------

source("./common.R")

# <<Figure 5>>
plot.phylo(tree_en, tip.color = cl, font = 1)
save_pdf("../data/derived/fig.5.pdf")
