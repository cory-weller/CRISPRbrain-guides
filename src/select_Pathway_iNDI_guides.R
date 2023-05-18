#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(org.Hs.eg.db)

## GENES IN ENRICHED PATHWAYS - GUIDE SELECTION ######################

enrichment_file <- 'outputs/pathway_enrichment_15_all.tsv'

pathways  <- fread(enrichment_file)
pathways <- pathways[qvalue < 0.05 & Count >= 3]

# Convert to long form (one gene ID per row) by splitting GeneID string
pathways <- pathways[, j=list('ENTREZID'=strsplit(geneID, split='/')[[1]]), 
                                    by=list(ID, qvalue, Description)]


gene_symbols <- select(
                    org.Hs.eg.db, 
                    keys = pathways$ENTREZID,
                    columns = "SYMBOL",
                    keytype = "ENTREZID"
                )

gene_symbols <- as.data.table(gene_symbols)
gene_symbols <- unique(gene_symbols)

pathways <- merge(pathways, gene_symbols, by='ENTREZID')
setnames(pathways, 'ID', 'PathwayID')
setnames(pathways, 'Description', 'PathwayDescription')
setnames(pathways, 'qvalue', 'PathEnrichQval')

setcolorder(pathways, c('SYMBOL','ENTREZID','PathwayID','PathEnrichQval','PathwayDescription'))

# Some genes are present in multiple pathways
pathways[, list('EnrichedPathwayCount'=.N), by=list(SYMBOL,ENTREZID)][order(-EnrichedPathwayCount)]
#      SYMBOL ENTREZID EnrichedPathwayCount
#  1:    TP53     7157                   36
#  2:     JUN     3725                   31
#  3:   CASP3      836                   30
#  4:     BAX      581                   27
#  5:   MAPK8     5599                   25
#  6:     RB1     5925                   17
#  7:   EP300     2033                   11
#  8: GADD45G    10912                   10
#  9:     APC      324                   10
# 10:    ATF4      468                    9
# 11:  CDC25A      993                    9
# 12: EIF2AK3     9451                    7
# 13:   NANOG    79923                    5
# 14:  CDKN1C     1028                    4
# 15:    DFFA     1676                    4
# 16:  POU5F1     5460                    4
# 17:     SKI     6497                    3
# 18:    ABL2       27                    2
# 19:   ATG2B    55102                    2
# 20:  AMBRA1    55626                    2
# 21:   CCDC6     8030                    2
# 22:  UBE2E3    10477                    1
# 23:    FRS3    10817                    1
# 24:    CHD4     1108                    1
# 25:    AHCY      191                    1
# 26:    CYS1   192668                    1
# 27:     HTT     3064                    1
# 28:  METTL3    56339                    1
# 29:   TAOK1    57551                    1
# 30:    RORB     6096                    1
# 31:   AURKA     6790                    1
# 32:    KLF5      688                    1
# 33: MAP3K12     7786                    1
# 34:   KAT6A     7994                    1
# 35:   MAGI1     9223                    1
# 36:   SOCS6     9306                    1

# 36 Genes encompassing the enriched pathways
enriched_pathways_genes <- unique(pathways$SYMBOL)

# Import file with top guide per gene
top_2_file <- 'data/20200513_library_1_2_unbalanced.csv'
top_2_guides <- fread(top_2_file, header=TRUE)

# Rename col and exclude neg control guides
setnames(top_2_guides, 'gene', 'SYMBOL')
top_2_guides <- top_2_guides[SYMBOL != 'negative_control']

# Key by gene symbol for ease of reading order
setkey(top_2_guides, SYMBOL)

# Subset the 36 genes from enriched pathways
enriched_pathways_guides <- top_2_guides[SYMBOL %in% enriched_pathways_genes]
fwrite(enriched_pathways_guides,
        file='outputs/enriched_pathways_guides.csv',
        quote=F,
        row.names=F,
        col.names=T,
        sep=','
)

# all 36 of the genes from eriched pathways are present
# 3 have two options: ABL2, APC, POU5F1




## iNDI GENE GUIDE SELECTION #########################################

iNDI_genes_file <- 'data/iNDI_genes_extended.csv'
iNDI_genes <- fread(iNDI_genes_file)$gene
# 85 iNDI genes

# Subset to include only iNDI genes
iNDI_guides <- top_2_guides[SYMBOL %in% iNDI_genes]
fwrite(iNDI_guides,
        file='outputs/iNDI_guides.csv',
        quote=F,
        row.names=F,
        col.names=T,
        sep=','
)

# Find genes not present in the top 2 guides table
missing <- iNDI_genes[! iNDI_genes %in% unique(top_2_guides[SYMBOL %in% iNDI_genes]$SYMBOL)]
cat(paste0('Genes not covered by guides: ', missing, '\n'))

## Combined simplified table #########################################
subset_cols <- c('SYMBOL','sgID_A','protospacer_A','sgID_B','protospacer_B')

dat <- rbindlist(list(
    enriched_pathways_guides[, c(.SD, 'group'='Pathway_Enrichmant'), .SDcols=subset_cols],
    iNDI_guides[, c(.SD, 'group'='iNDI'), .SDcols=subset_cols]
))

fwrite(dat, file='outputs/pathways_iNDI_guides_merged.csv',
        quote=F,
        row.names=F,
        col.names=T,
        sep=','
)
