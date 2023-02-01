#--------------------------
# Load data and configure
# globally used variablaes
# 25 January 2023
#--------------------------

# phenotypes
df <- read.csv("../data/derived/data.csv", row.names = 1)
# imputed <- read.csv("../data/derived/data_OU_imputed.csv", row.names = 1)

# phylogenetic tree
tree <- read.tree("../data/external/tree.nwk")

# external data used for validation
ashwell <- read.csv("../data/external/ashwell.csv", row.names = 1)
smaers <- read.csv("../data/external/smaers.csv", row.names = 1)
sultan <- read.csv("../data/external/sultan.csv", row.names = 1)

# tree taxa clusters
taxa <- c("Marsupials", "Eulipotyphla", "Scrotifera", "Glires", "Euarchonta",
        "Pilosa", "Paenungulata", "Afroinsectiphilia")
cl <- setNames(c(1,
        2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        7, 7,
        8, 8,
        6, 6, 6, 6),
        df$English_name[match(tree$tip.label, df$Binomial_name_timetree)])
scl <- taxa[cl]

# smaller tree without species missing data
dfs <- na.omit(df)
trees <- drop.tip(tree, c("Sorex_araneus", "Bos_indicus", "Equus_grevyi",
        "Ovis_aries", "Phoca_vitulina", "Ursus_maritimus", "Zalophus_californianus",
        "Pan_troglodytes"))

# tree with English names
tree_en <- tree
tree_en$tip.label <- df$English_name[
        match(tree$tip.label, df$Binomial_name_timetree)]
