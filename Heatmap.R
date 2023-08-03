library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
options(stringsAsFactors=FALSE)

# Import data 

# WCGNA modules
Disease_MS_CD_Mod <- read.csv("Disease_MS_CD_gene_modules.txt", sep="\t")
Disease_TST_Mod <- read.csv("Disease_TST_MS_CD_gene_modules.txt", sep="\t")

# Expression Data
Disease_MS_CD_Expr <- as.data.frame(read_csv("Disease_MS_CD_Gene_Expression.csv", show_col_types = FALSE))
row.names(Disease_MS_CD_Expr) <- Disease_MS_CD_Expr$gene_name
Disease_MS_CD_datExpr <- t(Disease_MS_CD_Expr[,-1])

Healthy_MS_CD_Expr<- as.data.frame(read_csv("Healthy_MS_CD_Gene_Expression.csv", show_col_types = FALSE))
row.names(Healthy_MS_CD_Expr) <- Healthy_MS_CD_Expr$gene_name
Healthy_MS_CD_datExpr <- t(Healthy_MS_CD_Expr[,-1])

Disease_TST_MS_CD_Expr<- as.data.frame(read_csv("Disease_TST_MS_CD_Gene_Expression.csv", show_col_types = FALSE))
row.names(Disease_TST_MS_CD_Expr) <- Disease_TST_MS_CD_Expr$gene_name
Disease_TST_datExpr <- t(Disease_TST_MS_CD_Expr[,-1])

# Set Heatmap colors
col_fun = colorRamp2(c(0, 7, 15), c('green', 'black', 'red'))
col_fun(seq(-3,3))

# Disease MS_CD Expression 
#png("Figures/Disease_MS_CD_Gene_Expression_Heatmap.png")
Heatmap(as.matrix(Disease_MS_CD_datExpr), 
        name='mat', 
        col=col_fun, 
        column_title = 'MS/CD Genes', 
        row_title = 'Disease Samples')
#dev.off()

# Modules of interest for Disease_MS_CD
modules_of_interest = c("purple", "greenyellow", "brown", "pink", "turquoise")
submod = Disease_MS_CD_Mod %>%
  subset(colors %in% modules_of_interest)
subexpr_MS_CD = subset(as.data.frame(Disease_MS_CD_datExpr), select = submod$gene_id)
png("Figures/Disease_MS_CD_Modules_of_Interest_Heatmap.png")
Heatmap(as.matrix(subexpr_MS_CD), 
        name='mat', 
        col=col_fun,
        column_title = 'MS/CD Modules of Interest Genes',
        row_title = 'Disease Samples')
dev.off()

# Healthy MS_CD Expression 
#png("Figures/Healthy_MS_CD_Gene_Expression_Heatmap.png")
Heatmap(as.matrix(Healthy_MS_CD_datExpr), 
        name='mat', 
        col=col_fun, 
        column_title = 'MS/CD Genes', 
        row_title = 'Healthy Samples')
#dev.off()

# Modules of interest for Disease_TST_MS_CD
modules_of_interest = c("red", "lightcyan", "black", "pink", "purple", "blue")
submod = Disease_TST_Mod %>%
  subset(colors %in% modules_of_interest)
subexpr_TST = subset(as.data.frame(Disease_TST_datExpr), select = submod$gene_id)
# png("Figures/Disease_TST_MS_CD_Modules_of_Interest_Heatmap.png")
Heatmap(as.matrix(subexpr_TST), 
        name='mat', 
        col=col_fun,
        column_title = 'TST/MS/CD Genes', 
        row_title = 'Disease Samples')
# dev.off()