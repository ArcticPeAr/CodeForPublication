 
library(clusterProfiler)
library(org.Hs.eg.db) 
library(DOSE)
library(enrichplot)
library(readxl)
library(gprofiler2)

## Loading data with EntrezID, gene_names, FDR, P-value and log fold change
# The samples are C7-T8, C1-T2, and C1-T8
dataC1T2 <- read_excel("/home/petear/PublicationPaper/C1-T2.xlsx")
dataC1T8 <- read_excel("/home/petear/PublicationPaper/C1-T8.xlsx")
dataC7T8 <- read_excel("/home/petear/PublicationPaper/C7-T8.xlsx")

## Extracting the EntrezID from the data
# ENTREZID could be replaced with gene_names
genesC1T2 <- dataC1T2$ENTREZID
genesC1T8 <- dataC1T8$ENTREZID
genesC7T8 <- dataC7T8$ENTREZID

## Converting the EntrezID to numeric
genesC1T2 <- as.numeric(genesC1T2)
genesC1T8 <- as.numeric(genesC1T8)
genesC7T8 <- as.numeric(genesC7T8)



## Enrichment analysis 
# Biological Process
ego_BP_C1T2 <- enrichGO(gene = genesC1T2,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_BP_C1T8 <- enrichGO(gene = genesC1T8,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_BP_C7T8 <- enrichGO(gene = genesC7T8,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Cellular Component
ego_CC_C1T2 <- enrichGO(gene = genesC1T2,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "CC", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_CC_C1T8 <- enrichGO(gene = genesC1T8,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "CC", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_CC_C7T8 <- enrichGO(gene = genesC7T8,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "CC", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Molecular Function

ego_MF_C1T2 <- enrichGO(gene = genesC1T2,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "MF", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_MF_C1T8 <- enrichGO(gene = genesC1T8,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "MF", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_MF_C7T8 <- enrichGO(gene = genesC7T8,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "MF", # Biological Process. Change to "CC" for Cellular Component, or "MF" for Molecular Function
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

## Pathway analysis
# KEGG

ego_KEGG_C1T2 <- enrichKEGG(gene = genesC1T2,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_KEGG_C1T8 <- enrichKEGG(gene = genesC1T8,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_KEGG_C7T8 <- enrichKEGG(gene = genesC7T8,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)


# gProfileR
# Biological Process
gprofiler2_BP_C1T2 <- gprofiler2(genesC1T2, 
                                 organism = "hsapiens",
                                 domain = "biological_process",
                                 significant = TRUE)
                                 