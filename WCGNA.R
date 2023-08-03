library(readr)
library(WGCNA)
library(tidyverse)
library(phyloseq)
options(stringsAsFactors=FALSE)

# Import data
Disease_MS_CD_Expr <- as.data.frame(read_csv("Disease_MS_CD_Gene_Expression.csv", show_col_types = FALSE))
row.names(Disease_MS_CD_Expr) <- Disease_MS_CD_Expr$gene_name
Healthy_MS_CD_Expr<- as.data.frame(read_csv("Healthy_MS_CD_Gene_Expression.csv", show_col_types = FALSE))
row.names(Healthy_MS_CD_Expr) <- Healthy_MS_CD_Expr$gene_name

Disease_TST_MS_CD_Expr<- as.data.frame(read_csv("Disease_TST_MS_CD_Gene_Expression.csv", show_col_types = FALSE))
row.names(Disease_TST_MS_CD_Expr) <- Disease_TST_MS_CD_Expr$gene_name
Healthy_TST_MS_CD_Expr <- as.data.frame(read_csv("Healthy_TST_MS_CD_Gene_Expression.csv", show_col_types = FALSE))
row.names(Healthy_TST_MS_CD_Expr) <- Healthy_TST_MS_CD_Expr$gene_name

# Transpose data
Disease_MS_CD_datExpr <- t(Disease_MS_CD_Expr[,-1])
Healthy_MS_CD_datExpr <- t(Healthy_MS_CD_Expr[,-1])

Disease_TST_datExpr <- t(Disease_TST_MS_CD_Expr[,-1])
Healthy_TST_datExpr <- t(Healthy_TST_MS_CD_Expr[,-1])


# Choose the soft-thresholding power 
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(Disease_MS_CD_datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 7) 
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
#png("Figures/Disease_MS_CD_Scale_Independence.png")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", 
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
#dev.off()

# Mean connectivity as a function of the soft-thresholding power
#png("Figures/Disease_MS_CD_Mean_Connectivity.png")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

# Check scale-free topology
png("Figures/Disease_MS_CD_Scale_Free_Plot.png")
ADJ1 = abs(cor(Disease_MS_CD_datExpr,use="p"))^4
k = as.vector(apply(ADJ1,2,sum,na.rm=T))
sizeGrWindow(10,5) 
par(mfrow=c(1,2)) 
hist(k)
scaleFreePlot(k,main="Check scale free topology")
dev.off()

picked_power = 4
temp_cor <- cor       
cor <- WGCNA::cor
netwk <- blockwiseModules(as.matrix(Disease_MS_CD_datExpr),            
                          # == Adjacency Function ==
                          power = picked_power,               
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          minModuleSize = 5,
                          maxBlockSize = 100,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
#png("Figures/Disease_MS_CD_Cluster_Dendrogram.png")
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)
#dev.off()

# Get gene list from modules
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
module_df <- module_df[order(module_df$colors, decreasing=TRUE),]
write_delim(module_df,
            file = "Disease_MS_CD_gene_modules.txt",
            delim = "\t")

# Get module eigengenes
MEs0 <- moduleEigengenes(as.matrix(Disease_MS_CD_datExpr), mergedColors)$eigengenes
# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
# Add sample names
MEs0$treatment = row.names(MEs0)
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

# Plot module-treatment relationships in a heatmap
# png("Figures/Disease_MS_CD_Module_Trait.png")
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
# dev.off()

# Pull out modules of interest
modules_of_interest = c("red", "lightcyan", "black", "pink")
submod = module_df %>%
  subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id

# Get expression for genes from modules of interest
subexpr = subset(as.data.frame(Disease_MS_CD_datExpr), select = submod$gene_id)
write.csv(subexpr, "Disase_MS_CD_modules_of_interest.csv")
write_delim(subexpr,
            file = "Disease_MS_CD_modules_of_interest.txt",
            delim = "\t")