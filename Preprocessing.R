library(readr)
library(dplyr)

# Import data
Disease_MS.df <- read_csv("Disease_mRNA_LncRNA_Multiple_Sclerose.csv",show_col_types = FALSE)
Healthy_MS.df <- read_csv("Healthy_mRNA_LncRNA_Multiple_Sclerose_HC.csv",show_col_types = FALSE)

# MS/CD Gene list
MS_CD_genes <- read.delim("MS_CD_Genes.txt")

# Testosterone Gene list
TST_genes <- read.delim("Testosterone_Genes.txt")

TST_genes.list <- c('CYP19A1', 'FSHB', 'GNAS', 'HSD17B3', 'IGSIF1', 'LHB', 'MSH5', 'NDNF', 'NR5A1', 'SRY', 'TCF12')
TST_genes.list <- append(TST_genes.list, TST_genes$Symbol)
TST_genes.list <- TST_genes.list[!duplicated(TST_genes.list)]


# Filter Disease/Healthy by shared MS_CD genes 
Disease_MS_mod.df <- Disease_MS.df[Disease_MS.df$gene_name %in% MS_CD_genes$Symbol,]
Healthy_MS_mod.df <- Healthy_MS.df[Healthy_MS.df$gene_name %in% MS_CD_genes$Symbol,]

# Filter Disease/Healthy by Testosterone genes
Disease_TST_mod.df <- Disease_MS.df[Disease_MS.df$gene_name %in% TST_genes.list,]
Healthy_TST_mod.df <- Healthy_MS.df[Healthy_MS.df$gene_name %in% TST_genes.list,]

# Limit testosterone addition to 100 - ranked by significance in file 
Disease_TST_mod.df <- Disease_TST_mod.df[row.names(Disease_TST_mod.df) %in% 1:100, ]
Healthy_TST_mod.df <- Healthy_TST_mod.df[row.names(Healthy_TST_mod.df) %in% 1:100, ]

# Combine MS_CD and Testosterone data frames
Disease_MS_TST.df <- rbind(Disease_MS_mod.df, Disease_TST_mod.df)
Healthy_MS_TST.df <- rbind(Healthy_MS_mod.df, Healthy_TST_mod.df)

# Remove duplicates
Disease_MS_TST_dup.df <- distinct(Disease_MS_TST.df, gene_name, .keep_all = TRUE)
Healthy_MS_TST_dup.df <- distinct(Healthy_MS_TST.df, gene_name, .keep_all = TRUE)

# Drop columns
drops <- c('ID', 'gene_type', 'chr', 'start', 'end')
Disease_MS_mod.df <- Disease_MS_mod.df[, !(names(Disease_MS_mod.df) %in% drops)]
Healthy_MS_mod.df <- Healthy_MS_mod.df[, !(names(Healthy_MS_mod.df) %in% drops)]

Disease_MS_TST_mod.df <- Disease_MS_TST_dup.df[, !(names(Disease_MS_TST_dup.df) %in% drops)]
Healthy_MS_TST_mod.df <- Healthy_MS_TST_dup.df[, !(names(Healthy_MS_TST_dup.df) %in% drops)]

# Export files 
write.csv(Disease_MS_mod.df, file='Disease_MS_CD_Gene_Expression.csv', row.names=FALSE, na='')
write.csv(Healthy_MS_mod.df, file='Healthy_MS_CD_Gene_Expression.csv', row.names=FALSE, na='')
write.csv(Disease_MS_TST_mod.df, file='Disease_TST_MS_CD_Gene_Expression.csv', row.names=FALSE, na='')
write.csv(Healthy_MS_TST_mod.df, file='Healthy_TST_MS_CD_Gene_Expression.csv', row.names=FALSE, na='')