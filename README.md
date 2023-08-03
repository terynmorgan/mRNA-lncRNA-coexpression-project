# mRNA/lncRNA Co-Expression Networks in Multiple Sclerosis
Project done through INFO-B 646: Computational Systems Biology 

**Spring 2023** <br/>
**Programming Language:** R <br/>
**Description:** <br/>
This proejct aimed to detect if an association exists between Multiple Sclerosis (MS) and Celiac Disease (CD) by identifying co-expressed modules of genes associated with both disorders using ***Weighted Gene Correlation Network Analysis (WGCNA)***. Using the known association between MS and testosterone, this analysis also sought to analyze patterns of co-expression in the presence and absence of genes associated with testosterone using the ***Reactome Cytoscape Plugin***.  

**Data:** <br/>
mRNA and lncRNA gene expression data for MS was extracted from PltDB, a blood-platelets-based gene expression database. Extracted data included a disease and healthy gene expression dataset with 86 and 354 samples respectively, each with approximately 31,000 genes. 

The NCBI Gene Tool was used to obtain a list of associated genes with MS and CD, totaling 127 genes. The NCBI Gene Tool and Genetic Testing Registry tools were used to obtain genes associated with testosterone. This resulted in 1246 genes from the NCBI gene tool and 13 with clinical relevance from the Genetic Testing Registry, totaling 1259 genes.  

**Required Files:** </br>
<li> Disease_mRNA_LncRNA_Multiple_Sclerose.csv -> Csv file of diseased Multiple Sclerosis (MS) gene expression, 86 samples and 31k genes </li> 
<li> Healthy_mRNA_LncRNA_Multiple_Sclerose_HC.csv -> Csv file of healthy gene expression, 354 samples and 31k genes </li> 
<li> MS_CD_Genes.txt -> Txt file of genes shared by MS and Celiac Disease (CD) from NCBI's gene tool (127 genes) </li> 
<li> Testosterone_Genes.txt -> Txt file of genes associated with Testosterone (TST) from NCBI's gene tool (1246 genes) </li> 

Preprocessing.R -> R file to filter Disease/Healthy samples by genes in MS_CD_Genes.txt and Testosterone_Genes.txt
WCGNA.R -> R file for WCGNA analysis for the Disease/Healthy samples in MS/CD and MS/CD/TST gene expression subsets
Heatmap.R -> R file to create heatmaps from WCGNA modules of interest 

Healthy_MS_CD.cys -> Cytoscape file of co-expression Reactome network generated from healthy MS/CD expression data
Healthy_TST_MS_CD.cys -> Cytoscape file of co-expression Reactome network generated from healthy MS/CD/TST expression data
Disease_TST_MS_CD.cys ->Cytoscape file of co-expression Reactome network generated from disease MS/CD/TST expression data
Disease_MS_CD.cys -> Cytoscape file of co-expression Reactome network generated from disease MS/CD expression data

**Required packages:** </br>
WGCNA, readr, dplyr, tidyverse, phyloseq, ggplot2, RColorBrewer, circlize

**Output Files:** </br>
Disease_MS_CD_Gene_Expression.csv -> Csv file of subset of disease expression data by MS/CD genes 
Disease_TST_MS_CD_Gene_Expression.csv -> Csv file of subset of disease expression data by MS/CD/TST genes 
Healthy_MS_CD_Gene_Expression.csv -> Csv file of subset of healthy expression data by MS/CD genes 
Healthy_TST_MS_CD_Gene_Expression.csv -> Csv file of subset of healthy expression data by MS/CD/TST genes 

Disease_MS_CD_gene_modules.txt -> Txt file of genes and module colors from WCGNA analysis of disease MS/CD gene expression data
Disease_MS_CD_modules_of_interest.csv -> Csv file of expression data from genes in modules of interest from disease MS/CD expression WCGNA analysis

Disease_TST_MS_CD_gene_modules.txt -> Txt file of genes and module colors from WCGNA analysis of disease MS/CD/TST gene expression data
Disease_TST_MS_CD_modules_of_interest.csv ->  Csv file of expression data from genes in modules of interest from disease MS/CD/TST expression WCGNA analysis

Healthy_MS_CD_gene_modules.txt -> Txt file of genes and module colors from WCGNA analysis of healthy MS/CD gene expression data
Healthy_TST_MS_CD_gene_modules.txt -> Txt file of genes and module colors from WCGNA analysis of healthy MS/CD/TST gene expression data

Figures -> Folder of figures generated from analysis (READ_ME included in folder) 
