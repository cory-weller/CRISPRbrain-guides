#!/usr/bin/env Rscript
library(data.table)
library(foreach)
library(ggplot2)
library(clusterProfiler)
organism <- 'org.Hs.eg.db'
# BiocManager::install(organism, character.only = TRUE)         # Run if not installed
library(organism, character.only = TRUE)
# BiocManager::install("rWikiPathways")                         # Run if not installed
library(rWikiPathways)

# Import CRISPRbrain data
dat.2 <- fread('data/Tian_et_al_2019_2.csv')
dat.15 <- fread('data/Tian_et_al_2020_15.csv')

# Replace spaces in col names with underscore
setnames(dat.2, gsub(' ', '_', colnames(dat.2)))
setnames(dat.15, gsub(' ', '_', colnames(dat.15)))


get_entrez_ids <- function(symbol_list) {
    out <- select(org.Hs.eg.db, 
        keys = symbol_list,
        columns = c("ENTREZID", "SYMBOL"),
        keytype = "SYMBOL")
    return(out$ENTREZID)
}

dat.2.pwe <- enrichWP(get_entrez_ids(dat.2[Hit_Class=='Positive Hit']$Gene), organism = 'Homo sapiens')
dat.15.pwe <- enrichWP(get_entrez_ids(dat.15[Hit_Class=='Positive Hit']$Gene), organism = 'Homo sapiens')

dat.2.pwe <- dat.2.pwe@result
setDT(dat.2.pwe)
dat.15.pwe <- dat.15.pwe@result
setDT(dat.15.pwe)

fwrite(dat.2.pwe, file='pathway_enrichment_2_all.tsv', quote=F, row.names=F, col.names=T, sep='\t')
fwrite(dat.15.pwe, file='pathway_enrichment_15_all.tsv', quote=F, row.names=F, col.names=T, sep='\t')



#### guide RNA scoring #############################################################################

# Import and massage all data into single table
CML_day6 <- as.data.table(openxlsx::read.xlsx('data/1-s2.0-S0092867422005979-mmc1.xlsx', sheet='TabB_K562_day6_library'))
CML_day8 <- as.data.table(openxlsx::read.xlsx('data/1-s2.0-S0092867422005979-mmc1.xlsx', sheet='TabA_K562_day8_library'))
RPE_day7 <- as.data.table(openxlsx::read.xlsx('data/1-s2.0-S0092867422005979-mmc1.xlsx', sheet='TabC_RPE1_day7_library'))
CML_day6_summary <- as.data.table(openxlsx::read.xlsx('data/1-s2.0-S0092867422005979-mmc2.xlsx', sheet='TabB_K562_day6_summary_stat'))
CML_day8_summary <- as.data.table(openxlsx::read.xlsx('data/1-s2.0-S0092867422005979-mmc2.xlsx', sheet='TabA_K562_day8_summary_stat'))
RPE_day7_summary <- as.data.table(openxlsx::read.xlsx('data/1-s2.0-S0092867422005979-mmc2.xlsx', sheet='TabC_RPE1_summary_statistic'))
CML_6 <- merge(CML_day6, CML_day6_summary, by.x='unique.sgRNA.pair.ID', by.y='genetic.perturbation')
CML_8 <- merge(CML_day8, CML_day8_summary, by.x='unique.sgRNA.pair.ID', by.y='genetic.perturbation')
RPE_7 <- merge(RPE_day7, RPE_day7_summary, by.x='unique.sgRNA.pair.ID', by.y='genetic.perturbation')
CML_6[, celltype := 'CML']
CML_8[, celltype := 'CML']
RPE_7[, celltype := 'RPE']
CML_6[, day := 6]
CML_8[, day := 8]
RPE_7[, day := 7]
CML_6[, batch := 'CML_6']
CML_8[, batch := 'CML_8']
RPE_7[, batch := 'CML_7']
dat <- rbindlist(list(CML_6, CML_8, RPE_7))
rm(CML_day6, CML_day8, RPE_day7, CML_day6_summary, CML_day8_summary, RPE_day7_summary, CML_6, CML_8, RPE_7)
gc()


# exclude rows with duplicated guides
dat <- dat[`either.guide.duplicated?` == FALSE]


# exclude data with p > 0.05
setnames(dat, 'energy.test.(p-value)', 'p')
dat <- dat[p <= 0.05]


# Get cluster IDs and gene members
clusters <- as.data.table(openxlsx::read.xlsx('data/1-s2.0-S0092867422005979-mmc3.xlsx', sheet='perturbation clusters')[,c("X1","members")])
cluster_members <- foreach(i=0:(nrow(clusters)-1), .combine='rbind') %do% {
    cluster_n <- clusters[X1==i]$X1
    members <- CJ('cluster'=cluster_n, 'Gene'=unlist(strsplit(clusters[X1==i]$members, split=',')))
}


# Iterate over clusters 0:63 and assign cluster IDs to `dat`
for(i in 0:63) {
    genes_in_cluster <- cluster_members[cluster==i]$Gene
    dat[gene %in% genes_in_cluster, 'cluster' := i]
}


# Exclude guides not assigned to a cluster
dat <- dat[! is.na(cluster)]


# Rank guides by percent.knockdown
dat[, 'cluster_rank' := frank(percent.knockdown, ties.method='min'), by=list(cluster, batch)]


# Reorder and write output
fwrite(dat[order(batch, cluster, cluster_rank)], file='guides_merged.tsv', quote=F, row.names=F, col.names=T, sep='\t')