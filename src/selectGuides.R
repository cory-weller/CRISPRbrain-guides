#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(org.Hs.eg.db)

enrichment_file <- 'outputs/pathway_enrichment_15_all.tsv'
guides_file <- 'outputs/guides_all.tsv'

pathways  <- fread(enrichment_file)

signif_pathways <- pathways[qvalue < 0.05 & Count >= 3]

# Convert to long form (one gene ID per row) by splitting GeneID string
pathways.long <- signif_pathways[, j=list('ENTREZID'=strsplit(geneID, split='/')[[1]]), 
                                    by=list(ID, qvalue, Description)]


gene_symbols <- select(
                    org.Hs.eg.db, 
                    keys = pathways.long$ENTREZID,
                    columns = "SYMBOL",
                    keytype = "ENTREZID"
                )

gene_symbols <- as.data.table(gene_symbols)
gene_symbols <- unique(gene_symbols)

pathways <- merge(pathways.long, gene_symbols, by='ENTREZID')
setnames(pathways, 'ID', 'PathwayID')
setnames(pathways, 'Description', 'PathwayDescription')
setnames(pathways, 'qvalue', 'PathEnrichQval')

setcolorder(pathways, c('SYMBOL','ENTREZID','PathwayID','PathEnrichQval','PathwayDescription'))

# Some genes are present in multiple pathways
geneID_counts <- pathways[, list('EnrichedPathwayCount'=.N), by=SYMBOL][order(-EnrichedPathwayCount)]

#     geneID  N
#  1:   7157 36
#  2:   3725 31
#  3:    836 30
#  4:    581 27
#  5:   5599 25
#  6:   5925 17
#  7:   2033 11
#  8:  10912 10
#  9:    324 10
# 10:    993  9
# 11:    468  9
# 12:   9451  7
# 13:  79923  5
# 14:   1676  4
# 15:   1028  4
# 16:   5460  4
# 17:   6497  3
# 18:     27  2
# 19:  55626  2
# 20:  55102  2
# 21:   8030  2
# 22: 192668  1
# 23:  56339  1
# 24:    191  1
# 25:   6096  1
# 26:  10477  1
# 27:   7994  1
# 28:   9223  1
# 29:   1108  1
# 30:   7786  1
# 31:  57551  1
# 32:  10817  1
# 33:   6790  1
# 34:    688  1
# 35:   3064  1
# 36:   9306  1

# Do we want top guides for a given gene, or top guides for a given pathway?

guides <- fread('outputs/guides_all.tsv')

setnames(guides, 'gene', 'SYMBOL')
# guides <- guides[, .SD, .SDcols=c('SYMBOL', 