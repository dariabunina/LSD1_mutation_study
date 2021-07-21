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

### ChromHMM overlaps ====
### compare differential peaks with all peaks in genomic features distribution:
# features are taken from the roadmap HMM in foreskin fibroblasts cell line E055 https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/indivModels/default_init/E055/n15/E055_15_coreMarks_dense.bed
# saved to ATACseq data folder:

# nhlf_hmm_data = import.bed(file.path(dir1, "bunina/LSD1/RNAseq/myRNAseqFibr/data", "wgEncodeBroadHmmNhlfHMM.bed"))
fib_hmm_data = rtracklayer::import.bed(file.path(dir1, "bunina/LSD1/ATACseq/data", "E055_15_coreMarks_dense.bed"))

#state names:
names_fib_hmm = read.table(file.path(dir1, "bunina/LSD1/ATACseq/data", "core15marks_HMM_names.txt"), sep = "\t")
names_fib_hmm_r = data.frame(ids = names_fib_hmm$V1[-1], numb = NA, names = NA, stringsAsFactors = FALSE)
names_fib_hmm_r$numb = sapply(strsplit(as.character(names_fib_hmm_r$ids), "_"), "[", 1)
names_fib_hmm_r$names = sapply(strsplit(as.character(names_fib_hmm_r$ids), "_"), "[", 2)

mcols(fib_hmm_data)$names_full = names_fib_hmm_r$names[match(mcols(fib_hmm_data)$name, names_fib_hmm_r$numb)]

# overlap peak midths only! to avoid multiple mappings to different features:
ATAC_all_n_mid = ATACnloe.r_all_gr
midths_atac = round(start(ATACnloe.r_all_gr) + width(ATACnloe.r_all_gr) / 2)
start(ATAC_all_n_mid) <- midths_atac
end(ATAC_all_n_mid) <- midths_atac

ol_atac_hmm = findOverlaps(ATAC_all_n_mid, fib_hmm_data)
mcols(ol_atac_hmm)$state = mcols(fib_hmm_data)$names_full[subjectHits(ol_atac_hmm)]

#now add the states to the original data granges - which has the same peak order as the midth GRanges:
mcols(ATACnloe.r_all_gr)$HMM_state = NA
mcols(ATACnloe.r_all_gr)$HMM_state[queryHits(ol_atac_hmm)] = mcols(ol_atac_hmm)$state[queryHits(ol_atac_hmm)]
mcols(ATACnloe.r_all_gr)$sigPeak = ifelse(mcols(ATACnloe.r_all_gr)$FDR < 0.05 & mcols(ATACnloe.r_all_gr)$Fold > 0, "UP_Kids", 
                                          ifelse(mcols(ATACnloe.r_all_gr)$FDR < 0.05 & mcols(ATACnloe.r_all_gr)$Fold < 0, "UP_Dads", "NonSignif"))

ATAC_all_n_bl_df = as.data.frame(ATACnloe.r_all_gr)
ATAC_hmm_stats = as.data.frame(table(ATAC_all_n_bl_df[,c(12:13)])) #16:17; 13:14

table(mcols(ATACnloe.r_all_gr)$sigPeak)
unique(mcols(ATACnloe.r_all_gr)$sigPeak)
total_peaks_df = data.frame(type = unique(mcols(ATACnloe.r_all_gr)$sigPeak), total = c(7364, 9061, 80041), stringsAsFactors = FALSE) #c(7048, 1207, 88119)

ATAC_hmm_stats$total = total_peaks_df$total[match(ATAC_hmm_stats$sigPeak, total_peaks_df$type)]

ATAC_hmm_stats$percent = 100*ATAC_hmm_stats$Freq / ATAC_hmm_stats$total
sum(ATAC_hmm_stats$percent[ATAC_hmm_stats$sigPeak == "NonSignif"]) #99.7%, good!
ATAC_hmm_stats = ATAC_hmm_stats[order(ATAC_hmm_stats$percent, decreasing = TRUE),]
ATAC_hmm_stats$HMM_state = factor(ATAC_hmm_stats$HMM_state, levels = unique(as.character(ATAC_hmm_stats$HMM_state)))

library(scales)
cbp1 <- c("#999999", scales::hue_pal()(2)[c(2,1)])

# pdf(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "ATACpeaksDistr_HMM_15states.pdf"), width = 4, height = 2.5, useDingbats = FALSE)
pdf(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "ATACpeaksDistr_HMM_15statesNew.pdf"), 
    width = 4, height = 2.5, useDingbats = FALSE)
ggplot(data = ATAC_hmm_stats, aes(x = HMM_state, y = percent)) + geom_bar(stat = "identity", aes(fill = sigPeak), position = "dodge", alpha = 0.9) +
  scale_fill_manual(values = cbp1) + theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#Fisher tests on selected states:
ATAC_hmm_stats$difference = ATAC_hmm_stats$total - ATAC_hmm_stats$Freq

hmm_fisher = list()
for (i in unique(ATAC_hmm_stats$HMM_state)) {
  for (j in c("UP_Dads", "UP_Kids")) {
    #order the rows properly to get correct oddsRatios directions:
    hmm_fisher[[paste(i, j, sep = "_")]] <- fisher.test(rbind(ATAC_hmm_stats[ATAC_hmm_stats$HMM_state == i & 
                                                                               ATAC_hmm_stats$sigPeak == j,c(3,6)], 
                                                              ATAC_hmm_stats[ATAC_hmm_stats$HMM_state == i & 
                                                                               ATAC_hmm_stats$sigPeak == "NonSignif",c(3,6)]))
  }
}
hmm_fisher_df = data.frame(id = names(hmm_fisher), pval = unlist(lapply(hmm_fisher, '[', 1)), 
                           oddR = unlist(lapply(hmm_fisher, '[', 3)), stringsAsFactors = FALSE)
hmm_fisher_df$p.adj = p.adjust(hmm_fisher_df$pval, method = "BH")
hmm_fisher_df$id = factor(hmm_fisher_df$id, levels = unique(hmm_fisher_df$id))

col_hmm = rep("black", nrow(hmm_fisher_df))
col_hmm[which(hmm_fisher_df$p.adj < 0.05)] = "red"

col_hmm1 = rep("#999999", nrow(hmm_fisher_df))
col_hmm1[which(hmm_fisher_df$p.adj < 0.05)] = "#E69F00"


pdf(file = file.path(dir1, "bunina/LSD1/ATACseq/output/plots", "ATACpeaksDistr_HMM_15states_pvalsNew.pdf"), 
    width = 4, height = 2.5, useDingbats = FALSE)
ggplot(data = hmm_fisher_df, aes(x = id, y = log2(oddR))) + geom_bar(stat = "identity", fill = col_hmm1) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(data = hmm_fisher_df, aes(x = id, y = log2(oddR))) + geom_point(color = col_hmm) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.2, alpha = 0.8) + 
  scale_y_continuous(position = "right")
dev.off()


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
