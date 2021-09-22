###########################
##ATAC-seq analysis ##
###########################
## Daria

library("DESeq2")
library("gplots")
library("ggplot2")
library("dplyr")
library(tidyverse)
library(IRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(biomaRt)
library(rtracklayer)
library(ggrepel)
library(viridis)
library(ChIPseeker)
library(DiffBind)

if (file.exists("/g/scb2/zaugg/bunina/LSD1/ATACseq/output/20160725_MomKid_ATACstats.csv")) {
  dir = "/g/scb2/zaugg"
} else {
  dir = "/Volumes/zaugg-1"
}
file.exists(file = file.path(dir, "bunina/LSD1/ATACseq/output", "20160725_MomKid_ATACstats.csv"))

#for DiffBind version 3.0.15:
### do loess correction, as the MA plot is very unbalanced!!
samples = read.csv(file.path(dir1, "bunina/LSD1/ATACseq/data", "sampleTableKidsDads_macs.csv"))
ATAC = dba(sampleSheet = samples)
ATAC = dba.count(ATAC) #done with old Diffbind version on the server, otherwise throws an error!
# saveRDS(ATAC, file = file.path(dir1, "bunina/LSD1/ATACseq/output/R.objects", "ATACkidsDadsCountsOnly.rds"))
ATACnloe = readRDS(file.path(dir1, "bunina/LSD1/ATACseq/output/R.objects", "ATACkidsDadsCountsOnly.rds"))
ATACnloe$config$AnalysisMethod <- DBA_EDGER
ATACnloe <- dba.normalize(ATACnloe, offsets=TRUE)
ATACnloe = dba.contrast(ATACnloe, design ="~Treatment+Condition") #account for batch effect
ATACnloe.a = dba.analyze(ATACnloe, bBlacklist = DBA_BLACKLIST_HG19)
saveRDS(ATACnloe.a, file = file.path(dir1, "bunina/LSD1/ATACseq/output/R.objects", "ATACloessNorm.rds"))
ATACnloe.a = readRDS(file.path(dir1, "bunina/LSD1/ATACseq/output/R.objects", "ATACloessNorm.rds"))
dba.show(ATACnloe.a, bContrasts = TRUE)

pdf(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "MAplotsLoessCorr.pdf"), width = 5, height = 4.5, useDingbats = FALSE)
dba.plotMA(ATACnloe.a, dotSize = 0.2) #method=edgeRGLM
dba.plotMA(ATACnloe.a, dotSize = 0.2, method = DBA_EDGER_BLOCK)
dev.off()

pdf(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "PCAplotsLoessCorr.pdf"), width = 5, height = 5, useDingbats = FALSE)
dba.plotPCA(ATACnloe.a, method = DBA_EDGER_BLOCK, dotSize = 1.5, label = DBA_ID)
dev.off()

dba.plotVolcano(ATACnloe.a, method = DBA_EDGER_BLOCK, dotSize = 0.2)
dba.plotVolcano(ATACnloe.a, dotSize = 0.2)

ATACnloe.r = dba.report(ATACnloe.a, method = DBA_EDGER_GLM) #DBA_EDGER_BLOCK
sum(ATACnloe.r$Fold > 0)
sum(ATACnloe.r$Fold < 0)
ATACnloe.r_all = dba.report(ATACnloe.a, th=1, DataType = DBA_DATA_FRAME) #doesn't work with GR! method = DBA_EDGER_BLOCK, 
ATACnloe.r_all_gr = makeGRangesFromDataFrame(ATACnloe.r_all[-62294,], keep.extra.columns = TRUE)



# diffTF_results_visualisation ----------------------------------------------------------

#new run of diffTF is here: /g/scb2/zaugg/bunina/LSD1/ATACseq/output/extension50
TFresults = read.table(file = file.path(dir1, "bunina/LSD1/ATACseq/output/extension50/LSD1_2019.summary.tsv"), 
                       header = T, sep = "\t")
hocoIDs = read.table(file = file.path(dir1, "bunina/LSD1/ATACseq/data", "HOCOTFID2ENSEMBL.txt"), header = T)
TFresults$ensembl = hocoIDs$ENSEMBL[match(TFresults$TF, hocoIDs$HOCOID)]
TFresults$type = ifelse(TFresults$weighted_meanDifference < 0, "UP_kids", "UP_dads")

# Filter out TFs not expressed in the fibroblasts:
TFresults$exprFib = ifelse(TFresults$ensembl %in% res_fibOrdered_df$ensembl, TRUE, FALSE)
table(TFresults$exprFib)
table(TFresults[,c(8,12)])
length(unique(TFresults$ensembl[TFresults$exprFib == TRUE]))
# write.csv(TFresults, file = file.path(dir1, "bunina/LSD1/ATACseq/output", "LSD1_2019.summaryEnsemblIDs.csv"), quote = F, row.names = F)
TFresults = read.csv(file.path(dir1, "bunina/LSD1/ATACseq/output", "LSD1_2019.summaryEnsemblIDs.csv"))

#stats of padj and cohenD:
TFresults$sig_padj = ifelse(TFresults$pvalueAdj < 0.05, TRUE, FALSE)
table(TFresults[,c("Cohend_factor", "sig_padj")])

TFresults_sig = TFresults[TFresults$Cohend_factor != 1 & TFresults$pvalueAdj < 0.05 & TFresults$exprFib == TRUE,]
write.csv(TFresults_sig, file = file.path(dir1, "bunina/LSD1/ATACseq/output", "LSD1_2019.summaryEnsemblIDs_signif.csv"), quote = F, row.names = F)
table(TFresults_sig$type)

TFresults_sig$TFname = egs_hg19_all$external_gene_name[match(TFresults_sig$ensembl, egs_hg19_all$ensembl_gene_id)]
TFresults$TFname = egs_hg19_all$external_gene_name[match(TFresults$ensembl, egs_hg19_all$ensembl_gene_id)]

library(clusterProfiler)
library(org.Hs.eg.db)
goObj_BP <- clusterProfiler::compareCluster(ensembl ~ type, data = TFresults_sig,
                                            fun = "enrichGO", universe = TFresults$ensembl,
                                            keyType = "ENSEMBL",
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP", readable = TRUE, 
                                            qvalueCutoff  = 0.1)
dotplot(goObj_BP, showCategory= 10)
goObj_BP_df = as.data.frame(goObj_BP@compareClusterResult)

goObj_BP_df = goObj_BP_df %>%
  dplyr::arrange(desc(Count)) #group_by(type) %>%

organDevel = goObj_BP_df$geneID[goObj_BP_df$ID == "GO:0048568"]
organDevel = unlist(sapply(strsplit(organDevel, "/"), "["))

AP1complex = c("JUN", "JUNB", "JUND", "FOS", "FOSL1", "FOSL2") #Fosl1/2 are Fra1/2!

TFresults_sig$GOterms = ifelse(TFresults_sig$TFname %in% AP1complex, "AP1complex",
                               ifelse(TFresults_sig$TFname %in% organDevel, "organDev", "other"))
table(TFresults_sig$GOterms)

### revert X-axis to match other plots direction LSD1mut vs LSD1wt:
TFresults_sig$LSD1mut_vs_WT = - TFresults_sig$weighted_meanDifference
TFresults$LSD1mut_vs_WT = - TFresults$weighted_meanDifference

TFresults$GOterms = ifelse(TFresults$TFname %in% AP1complex, "AP1complex",
                           ifelse(TFresults$TFname %in% organDevel, "organDev", "other"))

ageGenes_sig = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/output", "ageGenes_sig_pval.Rds"))
ageGenes300m = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/output", "ageGenes300m.Rds"))

TFresults$ageing = ifelse(TFresults$ensembl %in% c(ageGenes_sig$ensembl_gene_id, ageGenes300m$ensembl_gene_id), TRUE, FALSE)
TFresults_sig$ageing = ifelse(TFresults_sig$ensembl %in% c(ageGenes_sig$ensembl_gene_id, ageGenes300m$ensembl_gene_id), TRUE, FALSE)
table(TFresults_sig$ageing)
table(TFresults_sig[,c("ageing","type")])

pdf(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "diffTF2019_volcano_modJuly2021.pdf"), width = 7, height = 5, useDingbats = FALSE)
### remove low p-value TFs that have Cohen's D of 1 or are not expressed in fib and are not ageing:
TFresults %>%
  subset( (pvalueAdj > 0.05 & exprFib == TRUE) | (pvalueAdj < 0.05 & Cohend_factor != 1 & exprFib == TRUE) & ageing == FALSE) %>%
  ggplot() + theme_classic() + 
  geom_point(aes(x = LSD1mut_vs_WT, y = -log10(pvalueAdj), color = GOterms), alpha = 0.5) + 
  geom_text_repel(data = subset(TFresults_sig, GOterms %in% c("AP1complex", "organDev") & ageing == FALSE), aes(label = TF, x = LSD1mut_vs_WT, 
                                                                                                                y = -log10(pvalueAdj), color = GOterms), 
                  box.padding = 0.25, size = 2.5, segment.size = 0.1) +
  scale_color_manual(values = c("#5e3c99", "#e66101", "#404040")) +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2)
dev.off()


# Find nearest genes to peaks ---------------------------------------------

mart_hg19 = useMart(host='grch37.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ds_hg19 = useDataset('hsapiens_gene_ensembl', mart=mart_hg19)
egs_hg19 = getBM(attributes = c('ensembl_gene_id', 'chromosome_name','start_position', 'external_gene_name', 'end_position','strand'), mart=ds_hg19)

#take only prot coding and lncRNAs:
egs_hg19_prot = getBM(attributes = c('ensembl_gene_id', 'chromosome_name','start_position', 'external_gene_name', 'end_position','strand', 'entrezgene_id'), 
                      mart=mart_hg19, filters = "biotype", values = c("protein_coding", "lincRNA"))

egs_hg19$TSS = ifelse( egs_hg19$strand == "1", egs_hg19$start_position, egs_hg19$end_position )
egs_hg19$strand_m = ifelse(egs_hg19$strand == 1, "+", "-")
egs_hg19_tss = GRanges(seqnames = Rle( paste0('chr', egs_hg19$chromosome_name) ),
                       ranges = IRanges( start = egs_hg19$TSS,
                                         end = egs_hg19$TSS),
                       strand = Rle(egs_hg19$strand_m),
                       gene = egs_hg19$external_gene_name, ensembl = egs_hg19$ensembl_gene_id)

egs_hg19_tss = egs_hg19_tss[mcols(egs_hg19_tss)$ensembl %in% egs_hg19_prot$ensembl_gene_id]

#Find nearby protein-coding gene:
#(function taken from DiffBindSox2 script):
findNextProm = function(peakGR = sox2_diffP_s_downNP, promGR = egs_tss) {
  nearG = suppressWarnings(distanceToNearest(peakGR, promGR)) #GRanges function
  # length(unique(queryHits(nearG))) #all unique ones, good
  mcols(nearG)$ensembl = mcols(promGR)$ensembl[subjectHits(nearG)]
  mcols(nearG)$geneName = mcols(promGR)$geneName[subjectHits(nearG)]
  mcols(peakGR)$nearDist = NA
  mcols(peakGR)$nearDist[queryHits(nearG)] <- mcols(nearG)$distance
  mcols(peakGR)$nearGene = NA
  mcols(peakGR)$nearGene[queryHits(nearG)] <- mcols(nearG)$ensembl
  mcols(peakGR)$less50kb = ifelse(mcols(peakGR)$nearDist < 50000,"yes","no")
  mcols(peakGR)$less3kb = ifelse(mcols(peakGR)$nearDist < 3000,"yes","no")
  return(peakGR)
}

ATACnloe.r_all_gr
ATACnloe.r_al_prom = findNextProm(peakGR = ATACnloe.r_all_gr, promGR = egs_hg19_tss)

#add gene expression data:
dds_fib10 = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/output/Robjects",
                              "dds_fib10countsMinDadsKidsMyRNAseq.Rds"))
res_fib10 <- results(dds_fib10)
res_fibSig = subset(res_fib10, padj < 0.1 )
res_fibSig_df = as.data.frame(res_fibSig)
res_fibSig_df$geneType = ifelse(res_fibSig_df$log2FoldChange > 0, "up_kids", "down_kids")
res_fibSig_df$geneName = egs_hg19$external_gene_name[match(row.names(res_fibSig_df), egs_hg19$ensembl_gene_id)]
ageGenes300m = readRDS(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/output", "ageGenes300m.Rds"))
#do GO enrichments after removing ageing genes:
res_fibSig_df_noAge = subset(res_fibSig_df, !(row.names(res_fibSig_df) %in% ageGenes300m$ensembl_gene_id) )

ATACnloe.r_al_prom$geneLFC = res_fib10$log2FoldChange[match(ATACnloe.r_al_prom$nearGene, row.names(res_fib10))]
ATACnloe.r_al_prom$geneLFC_sig = ifelse(ATACnloe.r_al_prom$nearGene %in% row.names(res_fibSig_df_noAge), "sig", "other")
ATACnloe.r_al_prom$sig_both = ifelse(ATACnloe.r_al_prom$FDR < 0.05 & ATACnloe.r_al_prom$geneLFC_sig == "sig", "sig", "other")

# pdf(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "ATAC_RNAcorrScatter_new.pdf"), 
#     width = 4, height = 3.5, useDingbats = FALSE)
tiff(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "ATAC_RNAcorrScatter_new.tiff"), 
     width = 4, height = 3.5, units = "in", res = 500, compression = "lzw")
ggplot(data = as.data.frame(ATACnloe.r_al_prom), aes(x = Fold, y = geneLFC, color = sig_both)) + geom_point(size = 0.2, alpha = 0.5) + theme_classic() + 
  scale_color_manual(values = c("grey", "red")) + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + 
  xlab("ATACseq peak log2FC in LSD1mut vs. LSD1wt") + ylab("Peak-proximal RNA log2FC in LSD1mut vs. LSD1wt")
dev.off()

ATACnloe.r_al_prom_df = as.data.frame(ATACnloe.r_al_prom)
ATACnloe.r_al_prom_df = ATACnloe.r_al_prom_df[complete.cases(ATACnloe.r_al_prom_df[,c("Fold", "geneLFC")]),]
ATACnloe.r_al_prom_df$nearGene_name = egs_hg19$external_gene_name[match(ATACnloe.r_al_prom_df$nearGene, egs_hg19$ensembl_gene_id)]

#same for <3kb distance:
tiff(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "ATAC_RNAcorrScatter_new_less3kbTSS.tiff"), 
     width = 4, height = 3.5, units = "in", res = 500, compression = "lzw")
ggplot(data = subset(ATACnloe.r_al_prom_df, less3kb == "yes" ), aes(x = Fold, y = geneLFC, color = sig_both)) + geom_point(size = 0.2, alpha = 0.5) + theme_classic() + 
  scale_color_manual(values = c("grey", "red")) + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + 
  xlab("ATACseq peak log2FC in LSD1mut vs. LSD1wt") + ylab("Peak-proximal RNA log2FC in LSD1mut vs. LSD1wt") + 
  geom_text_repel(data = subset(ATACnloe.r_al_prom_df, less3kb == "yes" & sig_both == "sig" & Fold > 2.5 & geneLFC > 1.5), aes(label = nearGene_name), 
                  box.padding = 0.35, size = 3, segment.size = 0.2) + 
  geom_text_repel(data = subset(ATACnloe.r_al_prom_df, less3kb == "yes" & sig_both == "sig" & Fold < -1 & geneLFC < -1), aes(label = nearGene_name), 
                  box.padding = 0.35, size = 3, segment.size = 0.2)
dev.off()


cor(ATACnloe.r_al_prom_df$Fold, ATACnloe.r_al_prom_df$geneLFC)
cor(ATACnloe.r_al_prom_df$Fold[ATACnloe.r_al_prom_df$sig_both == "sig"], ATACnloe.r_al_prom_df$geneLFC[ATACnloe.r_al_prom_df$sig_both == "sig"])

ATACnloe.r_al_prom_df_3kb = subset(ATACnloe.r_al_prom_df, less3kb == "yes")
cor(ATACnloe.r_al_prom_df_3kb$Fold, ATACnloe.r_al_prom_df_3kb$geneLFC)
cor(ATACnloe.r_al_prom_df_3kb$Fold[ATACnloe.r_al_prom_df_3kb$sig_both == "sig"], ATACnloe.r_al_prom_df_3kb$geneLFC[ATACnloe.r_al_prom_df_3kb$sig_both == "sig"])
table(ATACnloe.r_al_prom_df_3kb$sig_both)
