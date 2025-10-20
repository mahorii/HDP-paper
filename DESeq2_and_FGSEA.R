#R version 4.4.1 (2024-06-14) -- "Race for Your Life"
#Copyright (C) 2024 The R Foundation for Statistical Computing
#Platform: aarch64-apple-darwin20

#DESeq2 analysis
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(4))
library("limma")

data <- read.csv("FeatureCounts_rawcounts_filt.CSV", header = TRUE, row.names = 1)

#DESeq2 data analysis
design <- read.csv("design.csv", row.names = 1, header=TRUE)
head(data)
head(design)
DESeq2CDS = DESeqDataSetFromMatrix(countData = data, colData=design, design=~batch + GA + sex + BMI_cat + disease)


DESeq2CDS_1 = DESeq(DESeq2CDS, parallel = TRUE)
resultsNames(DESeq2CDS_1)

# PE_vs_siPE
PE_vs_siPE = results(DESeq2CDS_1, contrast=c('disease', 'PE', 'siPE'))
PE_vs_siPE = PE_vs_siPE[order(PE_vs_siPE$padj),]
sum(PE_vs_siPE$padj < 0.05, na.rm=TRUE)
write.csv(PE_vs_siPE, file="PE_vs_siPE_account_for_GAwk_sex_BMI.csv")

#FGSEA
#Utility
library("dplyr")
library("tidyverse")
library("DT")
library("fgsea")
library("data.table")

#PE_vs_siPE_account_for_GA_sex_IUGR_BMI
res <- read.csv("PE_vs_siPE_account_for_GAwk_sex_BMI.csv")

#combine any symbols that are the same by averaging them
res2 <- res %>%
  dplyr::select(SYMBOL, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(SYMBOL) %>%
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)

pathways.hallmark <- gmtPathways("h.all.v7.0.symbols.gmt")

pathways.hallmark %>%
  head() %>%
  lapply(head)

#Run FGSEA
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=10000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>%
  DT::datatable()

#Write a table of results
fwrite(fgseaRes, file="PE_vs_siPE_account_for_GAwk_sex_BMI_FGSEA.txt")

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")+
  theme_minimal()