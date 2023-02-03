#--------------------------------------
# Load data, reformat it and configure
# globally used variables
# 25 January 2023
#--------------------------------------

# read phenotypes
cb_area_length <- read.csv("../data/derived/csv/01_cb_area_length.csv",
  row.names = 1)
ctx_area_length <- read.csv("../data/derived/csv/02_ctx_area_length.csv",
  row.names = 1)
cb_width <- read.csv("../data/derived/csv/03_cb_width.csv", row.names = 1)
cb_period <- read.csv("../data/derived/csv/04_cb_period.csv", row.names = 1)
cb_thickness <- read.csv("../data/derived/csv/05_cb_thickness.csv", row.names = 1)

# load rosetta file to dataframe
rosetta <- read.csv("../data/external/rosetta.csv")

# merge all phenotypes in a single dataframe based on row names
df <- merge(
  data.frame("Log10Area.cb"=cb_area_length$Log10Area,
             "Log10Length.cb"=cb_area_length$Log10Length,
             row.names =rownames(cb_area_length)),
  data.frame("Log10Area.ctx"=ctx_area_length$Log10Area,
             "Log10Length.ctx"=ctx_area_length$Log10Length,
             row.names=rownames(ctx_area_length)),
  by = "row.names")
df <- merge(df, data.frame("Log10WidthMedian"=log10(cb_width$WidthMedian), row.names=rownames(cb_width)), by.x = "Row.names", by.y = 0)
df <- merge(df, data.frame("Log10PeriodMedian"=log10(cb_period$PeriodMedian), row.names=rownames(cb_period)), by.x = "Row.names", by.y = 0)
df <- merge(df, data.frame("Log10ThicknessMedian"=log10(cb_thickness$ThicknessMedian), row.names = rownames(cb_thickness)), by.x = "Row.names", by.y = 0, all.x = TRUE)

# merge rosetta with df based on row names
df <- merge(df, data.frame(
  "Binomial_name"=rosetta$Binomial_name,
  "Binomial_name_timetree"=rosetta$Binomial_name_timetree,
  "Microdraw_name"=rosetta$Microdraw_name,
  "English_name"=rosetta$English_name,
  "Provenance"=rosetta$Provenance,
  "Log10BodyWeight"=rosetta$LogBodyWeight,
  "Log10BrainWeight"=rosetta$LogBrainWeight
  ), by.x = "Row.names", by.y = "Microdraw_name", all.x=FALSE, all.y=TRUE)
rownames(df)<-df$Binomial_name_timetree
colnames(df)[1] <- "Microdraw_name"
df <- df[order(row.names(df)), ]

# exclude failed estimations
df["Sorex_araneus","Log10PeriodMedian"]<-NaN
df["Sorex_araneus","Log10WidthMedian"]<-NaN
df["Bos_indicus","Log10ThicknessMedian"]<-NaN
df["Equus_grevyi","Log10ThicknessMedian"]<-NaN
df["Ovis_aries","Log10ThicknessMedian"]<-NaN
df["Pan_troglodytes","Log10ThicknessMedian"]<-NaN
df["Phoca_vitulina","Log10ThicknessMedian"]<-NaN
df["Ursus_maritimus","Log10ThicknessMedian"]<-NaN
df["Zalophus_californianus","Log10ThicknessMedian"]<-NaN

write.csv(df, "../data/derived/csv/data.csv")

# # objective:
# dx <- read.csv("../data/derived/csv/data.csv", row.names = 1)
# rownames(dx) <- dx$Binomial_name_timetree
# dx <- dx[order(row.names(dx)), ]
# summary(comparedf(df, dx))


# phylogenetic tree
tree <- read.tree("../data/external/tree.nwk")
name.check(tree, df)

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
